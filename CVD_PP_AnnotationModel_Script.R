#CVD-PP Variant annotation and model training script 
#Megan Ramaker#

#Load in required packages
library(pROC)
library(PRROC)
library(randomForest)
library(caret)
library(data.table)
library(plotly)
##Step 1: Variant Annotation
#Load in clinvar variants for 47 CVD genes 
Clinvar<-data.frame(fread("/Volumes/Sshare/Agenda/Ramaker/MonogenicDisease/CVD_Genes_Clinvar.txt",stringsAsFactors = F))
colnames(Clinvar)<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")

#Make new columns by extracting information from "INFO" column 
Clinvar$CLNSIG<-gsub("\\(|\\)|[0-9]|:", "", gsub(".*=","",gsub(";.*","", gsub(".*;CLNSIG.","",gsub("CLNSIGINCL=.*","",Clinvar$INFO)))))
Clinvar$CLNSIG_Simple<-NA
Clinvar$CLNSIG_Simple[grep("association", Clinvar$CLNSIG)]<-"Association"
Clinvar$CLNSIG_Simple[grep("not_provided", Clinvar$CLNSIG)]<-"Not_provided"
Clinvar$CLNSIG_Simple[grep("onfers_sensitivity", Clinvar$CLNSIG)]<-"Confers_sensitivity"
Clinvar$CLNSIG_Simple[grep("risk_factor", Clinvar$CLNSIG)]<-"Risk_factor"
Clinvar$CLNSIG_Simple[grep("Pathogenic", Clinvar$CLNSIG)]<-"Pathogenic"
Clinvar$CLNSIG_Simple[grep("Likely_pathogenic", Clinvar$CLNSIG)]<-"Likely_pathogenic"
Clinvar$CLNSIG_Simple[grep("Likely_benign", Clinvar$CLNSIG)]<-"Likely_benign"
Clinvar$CLNSIG_Simple[grep("Benign", Clinvar$CLNSIG)]<-"Benign"
Clinvar$CLNSIG_Simple[grep("Uncertain_significance", Clinvar$CLNSIG)]<-"Uncertain_significance"
Clinvar$CLNSIG_Simple[grep("Conflicting", Clinvar$CLNSIG)]<-"Conflicting_interpretations"
Clinvar$CLNSIG_Simple[grep("CLNSIG=Conflicting_interpretations_of_pathogenicity", Clinvar$INFO)]<-"Conflicting_interpretations"
Clinvar$RSID<-gsub(";.*","", gsub("ALLELEID.*","", gsub(".*;RS.","",Clinvar$INFO)))
Clinvar$RSID<-paste0("rs", Clinvar$RSID)
Clinvar$RSID[which(Clinvar$RSID=="rs")]<-NA
Clinvar$GENE<-gsub(":.*","",gsub(".*;GENEINFO.","",Clinvar$INFO))
Clinvar$LNREVSTAT<-gsub(".*criteria_provided,_","", gsub(";.*","",gsub(".*LNREVSTAT=","",Clinvar$INFO)))
unique(Clinvar$LNREVSTAT)
unique(Clinvar$CLNSIG)
table(Clinvar$CLNSIG_Simple)
unique(Clinvar$RSID)
table(Clinvar$CLNSIG_Simple[grep("conflicting_interpretations", Clinvar$LNREVSTAT)])
Clinvar$CLNSIG_Simple[which(is.na(Clinvar$CLNSIG_Simple))]<-"Uncertain_significance"
Clinvar$CLNSIG_Simple[which(Clinvar$CLNSIG_Simple=="Association"|Clinvar$CLNSIG_Simple=="Not_provided"|Clinvar$CLNSIG_Simple=="Confers_sensitivity"|Clinvar$CLNSIG_Simple=="Risk_factor")]<-"Uncertain_significance"
#table(Clinvar$CLNSIG_Simple)

#Read in CADD Scores and VEP results generated from the "simple" file created above
CADD<-data.frame(fread(file = "/Volumes/Sshare/Agenda/Ramaker/MonogenicDisease/CVD_Clinvar_GRCh38_CADD.tsv", header = T, stringsAsFactors = F))
VEP<-data.frame(fread(file = "/Volumes/Sshare/Agenda/Ramaker/MonogenicDisease/CVD_Clinvar_GRCh38_VEP.txt", header = T, stringsAsFactors = F))

#Make "Key" for each data set to merge
VEP$CHROM<-gsub(":.*$", "", VEP$Location)
VEP$START<-gsub("-.*","", gsub(".*:", "", VEP$Location))
VEP$STOP<-gsub("-.*$","", gsub(".*:", "", VEP$Location))
VEP$Key<-paste0(VEP$CHROM,"_", VEP$START, "_",VEP$Allele)  
CADD$Key<-paste0(CADD$X.Chrom,"_", CADD$Pos, "_", CADD$Alt)
Clinvar$Key<-paste0(Clinvar$CHROM, "_", Clinvar$POS, "_", Clinvar$ALT)
write.table(Clinvar, file="/Volumes/Sshare/Agenda/Ramaker/MonogenicDisease/ML_Clinvar_Dataset.txt", sep = "\t", quote = F, col.names = T, row.names = F)
#Merge data based on keys - remove duplicate entries
CADD<-CADD[which(CADD$Key %in% Clinvar$Key),]
CADD<-CADD[which(!duplicated(CADD$Key)),]
Clinvar<-Clinvar[which(Clinvar$Key %in% CADD$Key),]
Clinvar<-Clinvar[which(!duplicated(Clinvar$Key)),]
rownames(CADD)<-CADD$Key
rownames(Clinvar)<-Clinvar$Key
Clinvar<-Clinvar[order(rownames(Clinvar)),]
CADD<-CADD[order(rownames(CADD)),]
identical(rownames(Clinvar), rownames(CADD))

#Merge Clinvar and CADD
Merged_Annot<-merge(Clinvar,CADD)

#Make sure all genes have variants retained
unique(Merged_Annot$GeneName)
length(unique(Merged_Annot$GENE)) #48

#Merge VEP with Merged annot
VEP<-VEP[which(!duplicated(VEP$Key)),]
VEP<-VEP[which(VEP$Key%in%Merged_Annot$Key),]
Merged_Annot<-Merged_Annot[which(Merged_Annot$Key %in% VEP$Key),]
rownames(Merged_Annot)<-Merged_Annot$Key
rownames(VEP)<-VEP$Key
length(which(rownames(VEP)%in%rownames(Merged_Annot)))
VEP<-VEP[order(rownames(VEP)),]
Merged_Annot<-Merged_Annot[order(rownames(Merged_Annot)),]
identical(rownames(Merged_Annot), rownames(VEP))
Merged_Annot<-cbind(Merged_Annot,VEP)

#Remove what is not needed
rm(CADD,VEP,Clinvar)

#Make sure still have variants retained for 47 genes
length(unique(Merged_Annot$GENE)) 

#Gtex expression data subset in another script called CVD-PP gtex expression - merge with this script
#Load in Gtex expression for clinvar monogenic CVD variants - Key is column 12 and expression is column 8
#This has been lifted over to the new genome build GRCh38
Merged_Annot_Hits<-data.frame(fread("/Volumes/Sshare/Agenda/Ramaker/MonogenicDisease/Merged_Annot_GTF_liftover", sep = "\t", stringsAsFactors = F))
#Could have just written out the merged annot file in bed format and skipped next 4 lines but oh well
Merged_Annot$Expression<-0
for (i in 1:nrow(Merged_Annot_Hits)) {
  Merged_Annot$Expression[which(Merged_Annot$Key==Merged_Annot_Hits$V12[i])]<-Merged_Annot_Hits$V8[i]
}

#Calculate the distance to a nearest pathogenic variant 
#Seperate by pathogenicity
Pathogenic<-Merged_Annot[which(Merged_Annot$CLNSIG_Simple=="Pathogenic"|Merged_Annot$CLNSIG_Simple=="Likely_pathogenic"),] #n = 6639
Non_Pathogenic<-Merged_Annot[-which(Merged_Annot$CLNSIG_Simple=="Pathogenic"|Merged_Annot$CLNSIG_Simple=="Likely_pathogenic"),] #n = 53778
#Find closest pathogenic variant to each non path variant and record distance
Non_Pathogenic$DistToPathVar<-NA
for(i in 1:nrow(Non_Pathogenic)){
  tmp<-Pathogenic[which(Pathogenic$CHROM==Non_Pathogenic$CHROM[i]),]
  Non_Pathogenic$DistToPathVar[i]<-(abs(tmp$POS-Non_Pathogenic$POS[i]))[which.min(abs(tmp$POS-Non_Pathogenic$POS[i]))]
}
length(which(Non_Pathogenic$DistToPathVar==0))/nrow(Non_Pathogenic)

#Find next closest pathogenic variant to each pathogenic variant (i.e not itself)
Pathogenic$DistToPathVar<-NA
for(i in 1:nrow(Pathogenic)){
  print(i)
  tmp<-Pathogenic[which(Pathogenic$CHROM==Pathogenic$CHROM[i]&Pathogenic$POS!=Pathogenic$POS[i]),]
  if(nrow(tmp)==0){
    Pathogenic$DistToPathVar[i]<-0
  }else{
    Pathogenic$DistToPathVar[i]<-(abs(tmp$POS-Pathogenic$POS[i]))[which.min(abs(tmp$POS-Pathogenic$POS[i]))]
  }
}

#Combine pathogenic and non pathogenic variants back together
Merged_Annot<-data.frame(rbind(Pathogenic,Non_Pathogenic))
rm(tmp, Merged_Annot_Hits, Non_Pathogenic, Pathogenic, i)
table(Merged_Annot$CLNSIG_Simple)

#Annotations to use as predictors
#CADD Phred, GnomadAF, MaxEnScanDIFF (splicing score - use absolute value bc it's the difference between ref and alt), Distance to path var (manual calc) and exon expression from gtex
# Merged_Annot$PHRED
# Merged_Annot$gnomAD_AF
# abs(Merged_Annot$MaxEntScan_diff)
#Need to adjust the distributions and remove - in variables
#Sift - -> 1 #not using this
#PolyPheno - -> 0 #not using this
#Gnomad - -> 0 & multiple by 50000
#MaxEnt Diff - -> 0


#Reannotate genes bc should be 47
# bed<-Merged_Annot[,c("CHROM","POS","POS","Key")]
# bed$POS<-bed$POS-1
# bed$CHROM<-paste0("chr", bed$CHROM)
# #write.table(bed, file = "/Volumes/Sshare/Agenda/Ramaker/MonogenicDisease/Merged_Annot.bed", sep = "\t", quote = F, col.names = F, row.names = F)
# system("bedtools intersect -a /Volumes/Sshare/Agenda/Ramaker/MonogenicDisease/Merged_Annot.bed -b /Volumes/Sshare/Agenda/Ramaker/MonogenicDisease/CVD_Genes.bed -wa -wb > Merged_Annot_CVD_Genes_Intersect.bed")
bed<-data.frame(fread("/Volumes/Sshare/Agenda/Ramaker/MonogenicDisease/Merged_Annot_CVD_Genes_Intersect.bed", sep = "\t", stringsAsFactors = F))
rownames(bed)<-bed$V4
idx<-intersect(rownames(bed), rownames(Merged_Annot))
bed<-bed[idx,]
Merged_Annot<-Merged_Annot[idx,]
identical(rownames(bed), rownames(Merged_Annot))
Merged_Annot$GENE<-bed$V8
length(unique(Merged_Annot$GENE))
table(Merged_Annot$GENE)
Merged_Annot<-Merged_Annot[order(Merged_Annot$CHROM, Merged_Annot$POS),]
rm(bed, tmp, idx)
#Save image now- because this is before transforming values for random forest
#save.image("/Volumes/Sshare/Agenda/Ramaker/MonogenicDisease/Clinvar_RandomForest_Dataset.RData")

##Step 2: Model training
#Start here and load Clinvar_RandomForest_Dataset.RData file 
#load("/Volumes/Sshare/Agenda/Ramaker/MonogenicDisease/Clinvar_RandomForest_Dataset.RData")
table(Merged_Annot$CLNSIG_Simple)
#Need to adjust predictors to work in random forest
Merged_Annot$gnomAD_AF[which(Merged_Annot$gnomAD_AF=="-")]<-0
Merged_Annot$gnomAD_AF<-as.numeric(Merged_Annot$gnomAD_AF)
Merged_Annot$gnomAD_AF<-Merged_Annot$gnomAD_AF*50000 #scale AFs for RF
Merged_Annot$MaxEntScan_diff[which(Merged_Annot$MaxEntScan_diff=="-")]<-0
Merged_Annot$MaxEntScan_diff<-as.numeric(Merged_Annot$MaxEntScan_diff)
Merged_Annot$MaxEntScan_diff<-abs(Merged_Annot$MaxEntScan_diff)

#Create a clinical significance variable for random forest - should be Benign or Pathogenic
#Pathogenic and likely pathogenic -> Pathogenic
#Everything else is Benign
Merged_Annot$CLNSIG_RF<-"Benign"
Merged_Annot$CLNSIG_RF[grep("pathogenic", Merged_Annot$CLNSIG_Simple, ignore.case = T)]<-"Pathogenic"
table(Merged_Annot$CLNSIG_RF) #Benign = 53792, Pathogenic = 6625
tmp<-data.frame(table(Merged_Annot$CLNSIG_Simple))
table(Merged_Annot$CLNSIG_Simple, Merged_Annot$CLNSIG_RF) 

Merged_Annot$CLNSIG_Simple2<-"Benign"
Merged_Annot$CLNSIG_Simple2[grep("pathogenic", Merged_Annot$CLNSIG_Simple, ignore.case = T)]<-"Pathogenic"
Merged_Annot$CLNSIG_Simple2[grep("Uncertain", Merged_Annot$CLNSIG_Simple, ignore.case = T)]<-"VUS"
table(Merged_Annot$CLNSIG_Simple2)

#Benign (4173) Conflicting_interpretations (6223) Likely_benign (15833)          
#Likely_pathogenic (3846) Pathogenic (2793) Uncertain_significance (27549)

#Load in variants for which Benign/pathogenic status has changed from 2014 to 2022 (each year) (n = 663)
load("/Volumes/Sshare/Agenda/Ramaker/MonogenicDisease/clinvar_changes_2014-2022_annotated.RData")

#Keep a "Master" annotation dataframe before removing variants that have changed
Merged_Annot_Master<-Merged_Annot

#Remove variants that changes from benign -> path or Path -> benign
Merged_Annot<-Merged_Annot[-which(rownames(Merged_Annot)%in%rownames(changes)),]
dim(Merged_Annot) #59,753 variants remain


#Subset to variables for model
Merged_Annot_RF<-Merged_Annot[,c("CLNSIG_RF","PHRED","gnomAD_AF","MaxEntScan_diff","DistToPathVar",
                                 "LNREVSTAT", "Expression", "CLNSIG_Simple")]
dim(Merged_Annot_RF)#59753 variants

#Define test set - 33,782 of VUS and Conflicting interp. variants NOT in old clinvar with changes 
Discovery_Merged_Annot<-Merged_Annot[which(Merged_Annot$CLNSIG_Simple=="Uncertain_significance"|Merged_Annot$CLNSIG_Simple=="Conflicting_interpretations"),]
Discovery_Merged_Annot_RF<-Merged_Annot_RF[which(Merged_Annot_RF$CLNSIG_Simple=="Uncertain_significance"|Merged_Annot_RF$CLNSIG_Simple=="Conflicting_interpretations"),]
dim(Discovery_Merged_Annot_RF) #33,782  variants
table(Discovery_Merged_Annot_RF$CLNSIG_Simple) #Conflicting Int. = 6,223, VUS = 27,559

#Define train set - 25,971 Variants : Benign (19,974), Likely Ben. (15,833), Likely Path. (3,706), Pathogenic (2,681)
Merged_Annot_RF<-Merged_Annot_RF[-which(Merged_Annot_RF$CLNSIG_Simple=="Uncertain_significance"|Merged_Annot_RF$CLNSIG_Simple=="Conflicting_interpretations"),]
dim(Merged_Annot_RF) #25971 Variants
table(Merged_Annot_RF$CLNSIG_RF)#19,974 Benign, 5,997 Pathogenic
table(Merged_Annot_RF$CLNSIG_Simple)#Benign (4167), Likely Ben. (15807), Likely Path. (3414), Pathogenic (2583)

summary(Merged_Annot_RF$gnomAD_AF)
#Save image for easier analysis downstream
#save.image(file = "/Volumes/Sshare/Agenda/Ramaker/MonogenicDisease/MonogenicDisease_MachineLearning_Input_06282023.RData")
load(file = "/Volumes/Sshare/Agenda/Ramaker/MonogenicDisease/MonogenicDisease_MachineLearning_Input_06282023.RData")

#Make random forest output variable numeric
Merged_Annot_RF$CLNSIG_RF[which(Merged_Annot_RF$CLNSIG_RF=="Benign")]<-0
Merged_Annot_RF$CLNSIG_RF[which(Merged_Annot_RF$CLNSIG_RF=="Pathogenic")]<-1
Merged_Annot_RF$CLNSIG_RF<-as.numeric(Merged_Annot_RF$CLNSIG_RF)
table(Merged_Annot_RF$CLNSIG_RF)
table(Merged_Annot_RF$CLNSIG_Simple)
class(Merged_Annot_RF$CLNSIG_RF)

#Define number of trees to grow and random seeds 
trees<-c(25,50,100,250,500,750)
seeds<-c(5,23,604,799,50501)

#Test out
#trees<-c(10,20)
#seeds<-c(5,23)
#Random forest loop using 6 fold cross validation evenly split between pathogenic and benign variants
RF_trees<-list()
Pred_trees<-list()
ROC_AUC_trees<-list()
MSE_trees<-list()
Importance_trees<-list()
OptCutoff_trees<-list()
for(tree in trees){
  RF_seeds<-list()
  Pred_seeds<-list()
  ROC_AUC_seeds<-list()
  MSE_seeds<-list()
  Importance_seeds<-list()
  OptCutoff_seeds<-list()
  for (seed in seeds) {
    RF_k<-list()
    Pred_k<-list()
    Importance_k<-list()
    ROC_AUC_k<-list()
    MSE_k<-list()
    OptCutoff_k<-list()
    #Randomly shuffle the data
    path<-Merged_Annot_RF[which(Merged_Annot_RF$CLNSIG_RF==1),]
    set.seed(seed)
    path<-path[sample(nrow(path)),]
    path_fold<-cut(seq(1,nrow(path)), breaks = 6, labels = F)
    ben<-Merged_Annot_RF[which(Merged_Annot_RF$CLNSIG_RF!=1),]
    set.seed(seed)
    ben<-ben[sample(nrow(ben)),]
    ben_fold<-cut(seq(1,nrow(ben)), breaks = 6, labels = F)
    #Perform 10 fold cross validation
    for(k in 1:6){
      #Segement your data by fold using the which() function 
      print(k)
      path_testIndexes <- which(path_fold==k,arr.ind=TRUE)
      ben_testIndexes <- which(ben_fold==k,arr.ind=TRUE)
      testPath <- path[path_testIndexes,]
      testBen <- ben[ben_testIndexes,]
      trainPath <- path[-path_testIndexes,]
      trainBen <- ben[-ben_testIndexes,]
      train<-rbind(trainPath, trainBen)
      test<-rbind(testPath,testBen)
      print(tree)
      set.seed(seed)
      RF<-randomForest(x = train[,c(2,3,4,5,7)], y = train[,1], ntree = tree, importance = T)
      Pred<-predict(RF, newdata = test[,c(2,3,4,5,7)])
      Truth<-test[,1]
      ROC_AUC<-roc(as.numeric(Truth),as.numeric(Pred))$auc[1]
      OptCutoff<-Pred[which.min(unlist(lapply(Pred,function(x) sum(Pred[which(Truth==0)]>(x-0.0001))+sum(Pred[which(Truth==1)]<(x-0.0001)))))]
      MSE<-mean((Truth-Pred)^2)
      RF_k[[paste0(k)]]<-RF
      Pred_k[[paste0(k)]]<-Pred
      ROC_AUC_k[[paste0(k)]]<-ROC_AUC
      MSE_k[[paste0(k)]]<-MSE
      Importance_k[[paste0(k)]]<-RF$importance
      OptCutoff_k[[paste0(k)]]<-OptCutoff
    }
    RF_seeds[[paste0(seed)]]<-RF_k
    Pred_seeds[[paste0(seed)]]<-Pred_k
    ROC_AUC_seeds[[paste0(seed)]]<-ROC_AUC_k
    MSE_seeds[[paste0(seed)]]<-MSE_k
    Importance_seeds[[paste0(seed)]]<-Importance_k
    OptCutoff_seeds[[paste0(seed)]]<-OptCutoff_k
  }
  RF_trees[[paste0(tree)]]<-RF_seeds
  Pred_trees[[paste0(tree)]]<-Pred_seeds
  ROC_AUC_trees[[paste0(tree)]]<-ROC_AUC_seeds
  MSE_trees[[paste0(tree)]]<-MSE_seeds
  Importance_trees[[paste0(tree)]]<-Importance_seeds
  OptCutoff_trees[[paste0(tree)]]<-OptCutoff_seeds
}

#Step 3: Model interrogation 

#Top tree should have the maximum median ROC AUC across all seeds and splits - 750
Top_Tree<-trees[which.max(unlist(lapply(1:6, function(x) median(unlist(ROC_AUC_trees[[x]])))))] 

boxplot(unlist(ROC_AUC_trees$`25`), unlist(ROC_AUC_trees$`50`), unlist(ROC_AUC_trees$`100`),
        unlist(ROC_AUC_trees$`250`), unlist(ROC_AUC_trees$`500`), unlist(ROC_AUC_trees$`750`),
        ylab = "ROC AUC", xlab = "Number of Trees", names = names(ROC_AUC_trees))

#N tree = 750, median AUC = 0.9968
#Top model should use the top number of trees and have the median ROC AUC amongst all models using that n tree
Top_Model<-names(sort(unlist(ROC_AUC_trees[[paste0(Top_Tree)]]))[16]) #AUC = 0.9968
Top_Model # "5.5"
Top_RF<-RF_trees$`750`$`5`$`5`

#Define optimal naive bayes classifier cut off
OptCutoff<-OptCutoff_trees$`750`$`5`$`5` #0.55678945

#Pull scores from all models with top n trees - 5
clinvar_predictions<-data.frame(unlist(Pred_trees[[paste0(Top_Tree)]]))
#Remove names from random forest list from rownames and save as key column
clinvar_predictions<-cbind(gsub("^.*\\.","",rownames(clinvar_predictions)),clinvar_predictions)
rownames(clinvar_predictions)<-NULL
colnames(clinvar_predictions)<-c("Variant","Prediction")

#Make a table of predictions and pull most voted prediction across seeds 
tmp<-data.frame(table(clinvar_predictions$Variant))
tmp$score<-0
for(i in 1:nrow(tmp)){
  print(i)
  tmp$score[i]<-median(clinvar_predictions$Prediction[which(clinvar_predictions$Variant==tmp$Var1[i])])
}
tmp$Prediction<-0
tmp$Prediction[which(tmp$score>OptCutoff)]<-1
table(tmp$Prediction)
clinvar_predictions<-tmp[!duplicated(tmp$Var1),]
colnames(clinvar_predictions)<-c("Variant","NumVotes","MedianScore","Prediction")
rownames(clinvar_predictions)<-clinvar_predictions$Variant
rm(tmp)
#Combine predictions with Merged Annot (just because..)
idx<-intersect(rownames(Merged_Annot),rownames(clinvar_predictions))
Merged_Annot<-Merged_Annot[idx,]
clinvar_predictions<-clinvar_predictions[idx,]
identical(rownames(Merged_Annot), rownames(clinvar_predictions))
Merged_Annot$Pred<-clinvar_predictions$MedianScore
Merged_Annot$Prediction<-clinvar_predictions$Prediction


#Pull importance scores from all models using top n trees and calculate the mean importance
importance_scores<-list()
for(i in 1:5){
  for (j in 1:6) {
    importance_scores[[paste0(i, "_", j)]]<-RF_trees[[paste0(Top_Tree)]][[i]][[j]]$importance
  }
}
Importance_Scores_DF<-data.frame(do.call(cbind, importance_scores))
MeanImportance<-data.frame(do.call(cbind, list(rowMeans(Importance_Scores_DF[,grep("IncMSE", colnames(Importance_Scores_DF))]),
                                               rowMeans(Importance_Scores_DF[,grep("IncNodePurity", colnames(Importance_Scores_DF))]))))

colnames(MeanImportance)<-c("IncMSE","IncNodePurity")

#Step 4: Model Validation
#Predict variants for which clinical significance changed from 2014, i.e the "changes" data - 241 variants
changes$Pred<-predict(Top_RF, newdata = changes[,c("PHRED","gnomAD_AF","MaxEntScan_diff","DistToPathVar","Expression")])
changes$Prediction<-0
changes$Prediction[which(changes$Pred>OptCutoff)]<-1
table(changes$CLNSIG_Simple,changes$Prediction)
changes$CLNSIG_RF[which(changes$CLNSIG_RF=="Benign")]<-0
changes$CLNSIG_RF[which(changes$CLNSIG_RF=="Pathogenic")]<-1
table(changes$CLNSIG_RF)
changes$CLNSIG_Simple2<-"VUS"
changes$CLNSIG_Simple2[grep("pathogenic", changes$CLNSIG_Simple, ignore.case = T)]<-"Pathogenic"
changes$CLNSIG_Simple2[grep("benign", changes$CLNSIG_Simple, ignore.case = T)]<-"Benign"
changes$Previous_CLNSIG_Simple2<-"VUS"
changes$Previous_CLNSIG_Simple2[grep("pathogenic", changes$Previous_CLNSIG_Simple, ignore.case = T)]<-"Pathogenic"
changes$Previous_CLNSIG_Simple2[grep("benign", changes$Previous_CLNSIG_Simple, ignore.case = T)]<-"Benign"

table(changes$CLNSIG_Simple2, changes$Previous_CLNSIG_Simple2)

#First test performance on all changes variants
pr.curve(changes$CLNSIG_RF,changes$Pred)
pr.curve(changes$CLNSIG_RF,changes$Prediction)

#Annotate CADD binary classifier
changes$CADD_Validation<-NA
changes$CADD_Validation[which(changes$PHRED>=15)]<-1
changes$CADD_Validation[which(changes$PHRED<15)]<-0

changes$SIFT_Validation<-NA
changes$SIFT_Validation[which(changes$SIFTcat=="deleterious")]<-1
changes$SIFT_Validation[which(changes$SIFTcat=="tolerated")]<-0
changes$SIFTval #good
length(which(is.na(changes$SIFTval)))

changes$PolyPhen_Validation<-NA
changes$PolyPhen_Validation[grep("damaging",changes$PolyPhenCat)]<-1
changes$PolyPhen_Validation[grep("benign",changes$PolyPhenCat)]<-0
changes$PolyPhenVal #good
length(which(is.na(changes$PolyPhenVal)))

changes$MetaSVM_pred_Validation<-NA
changes$MetaSVM_pred_Validation[which(changes$MetaSVM_pred=="D")]<-1
changes$MetaSVM_pred_Validation[which(changes$MetaSVM_pred=="T")]<-0
changes$MetaSVM_score<-as.numeric(changes$MetaSVM_score) #good but a lot of NAs
length(which(is.na(changes$MetaSVM_score)))

changes$MetaLR_pred_Validation<-NA
changes$MetaLR_pred_Validation[which(changes$MetaLR_pred=="D")]<-1
changes$MetaLR_pred_Validation[which(changes$MetaLR_pred=="T")]<-0
changes$MetaLR_score<-as.numeric(changes$MetaLR_score) #good but a lot of NAs
length(which(is.na(changes$MetaLR_score)))

changes$MutationTaster_pred_Validation<-NA
changes$MutationTaster_pred_Validation[grep("N|P", changes$MutationTaster_pred)]<-0
changes$MutationTaster_pred_Validation[grep("A|D", changes$MutationTaster_pred)]<-1
changes$MutationTaster_converted_rankscore<-as.numeric(changes$MutationTaster_converted_rankscore) #good but NAs
length(which(is.na(changes$MutationTaster_converted_rankscore)))

changes$MutationAssessor_pred_Validation<-NA
changes$MutationAssessor_pred_Validation[grep("L|N", changes$MutationAssessor_pred)]<-0
changes$MutationAssessor_pred_Validation[grep("H|M", changes$MutationAssessor_pred)]<-1
changes$MutationAssessor_rankscore<-as.numeric(changes$MutationAssessor_rankscore) #good but NAs
length(which(is.na(changes$MutationAssessor_rankscore)))

changes$PROVEAN_Validation<-NA
changes$PROVEAN_Validation[grep("N", changes$PROVEAN_pred)]<-0
changes$PROVEAN_Validation[grep("D", changes$PROVEAN_pred)]<-1
changes$PROVEAN_converted_rankscore<-as.numeric(changes$PROVEAN_converted_rankscore) #good but NAs
length(which(is.na(changes$PROVEAN_converted_rankscore)))

changes$fathmm.XF_coding_pred_Validation<-NA
changes$fathmm.XF_coding_pred_Validation[grep("N", changes$fathmm.XF_coding_pred)]<-0
changes$fathmm.XF_coding_pred_Validation[grep("D", changes$fathmm.XF_coding_pred)]<-1
changes$fathmm.XF_coding_score<-as.numeric(changes$fathmm.XF_coding_score) #good but NAs
length(which(is.na(changes$fathmm.XF_coding_score)))

changes$BayesDel_addAF_pred_Validation<-NA
changes$BayesDel_addAF_pred_Validation[grep("T", changes$BayesDel_addAF_pred)]<-0
changes$BayesDel_addAF_pred_Validation[grep("D", changes$BayesDel_addAF_pred)]<-1
changes$BayesDel_addAF_score<-as.numeric(changes$BayesDel_addAF_score) #good but NAs
length(which(is.na(changes$BayesDel_addAF_score)))

changes$REVEL_rankscore<-as.numeric(changes$REVEL_rankscore)
changes$REVEL_rankscore_Validation<-NA
changes$REVEL_rankscore_Validation[which(changes$REVEL_rankscore>=0.75)]<-1
changes$REVEL_rankscore_Validation[which(changes$REVEL_rankscore<0.75)]<-0

#Add REVEL and CardioBoost
#Make a table for cardio boost submission
changes$cardioboost_var<-paste0(changes$GENE,":c.",changes$cDNApos,changes$REF,">",changes$ALT)
#write.table(paste0(changes$GENE,":c.",changes$cDNApos,changes$REF,">",changes$ALT),
#            file = "/Volumes/Sshare/Agenda/Ramaker/MonogenicDisease/changes_cardioboost_input.txt", sep = "\t", quote = F, row.names = F, col.names = F)

cardioboost_cm<-data.frame(fread("/Volumes/Sshare/Agenda/Ramaker/MonogenicDisease/cardiacboost_cm_predict.csv", sep = ",", stringsAsFactors = F))
table(cardioboost_cm$Classification.with.90..confidence.level, useNA = "always")
cardioboost_cm<-cardioboost_cm[-which(is.na(cardioboost_cm$Probablity.of.pathogenicity)),]

cardioboost_arm<-data.frame(fread("/Volumes/Sshare/Agenda/Ramaker/MonogenicDisease/cardiacboost_arm_predict.csv", sep = ",", stringsAsFactors = F))
table(cardioboost_arm$Classification.with.90..confidence.level, useNA = "always")
cardioboost_arm<-cardioboost_arm[-which(is.na(cardioboost_arm$Probablity.of.pathogenicity)),]

identical(cardioboost_arm$Upload_variant, changes$cardioboost_var)
changes[,paste0(colnames(cardioboost_arm)[5:6],"_arm")]<-NA
changes[,paste0(colnames(cardioboost_cm)[5:6],"_cm")]<-NA
length(which(changes$cardioboost_var%in%cardioboost_cm$Upload_variant))
for(i in 1:nrow(changes)){
  if(changes$cardioboost_var[i]%in%cardioboost_cm$Upload_variant){
    changes[i,paste0(colnames(cardioboost_cm)[5:6],"_cm")]<-cardioboost_cm[which(cardioboost_cm$Upload_variant==changes$cardioboost_var[i]),5:6]
  }
}
changes$cardioboost_cm_validation<-NA
changes$cardioboost_cm_validation[grep("Path",changes$Classification.with.90..confidence.level_cm)]<-1
changes$cardioboost_cm_validation[grep("Benign",changes$Classification.with.90..confidence.level_cm)]<-0
changes$cardioboost_cm_validation[grep("VUS",changes$Classification.with.90..confidence.level_cm)]<-0

for(i in 1:nrow(changes)){
  if(changes$cardioboost_var[i]%in%cardioboost_arm$Upload_variant){
    changes[i,paste0(colnames(cardioboost_arm)[5:6],"_arm")]<-cardioboost_arm[which(cardioboost_arm$Upload_variant==changes$cardioboost_var[i]),5:6]
  }
}

changes$cardioboost_arm_validation<-NA
changes$cardioboost_arm_validation[grep("Path",changes$Classification.with.90..confidence.level_arm)]<-1
changes$cardioboost_arm_validation[grep("Benign",changes$Classification.with.90..confidence.level_arm)]<-0
changes$cardioboost_arm_validation[grep("VUS",changes$Classification.with.90..confidence.level_arm)]<-0

##CM - 27 variants (2 benign, 25 pathogenic)
#Linear predictions
#CM - 0.9622597, CVD-PP - 0.9622597
pr.curve(as.numeric(changes$CLNSIG_RF[-which(is.na(changes$Probablity.of.pathogenicity_cm))]),as.numeric(changes$Probablity.of.pathogenicity_cm[-which(is.na(changes$Probablity.of.pathogenicity_cm))]))
pr.curve(as.numeric(changes$CLNSIG_RF[-which(is.na(changes$Probablity.of.pathogenicity_cm))]),as.numeric(changes$Pred[-which(is.na(changes$Probablity.of.pathogenicity_cm))]))
#Binary predictions
#CM - 0.7252452, CVD-PP - 0.9622597
pr.curve(as.numeric(changes$CLNSIG_RF[-which(is.na(changes$cardioboost_cm_validation))]),as.numeric(changes$cardioboost_cm_validation[-which(is.na(changes$cardioboost_cm_validation))]))
pr.curve(as.numeric(changes$CLNSIG_RF[-which(is.na(changes$cardioboost_cm_validation))]),as.numeric(changes$Pred[-which(is.na(changes$cardioboost_cm_validation))]))

#36 variants (1 benign, 35 pathogenic)
#Linear predictions
#ARM - 0.9591049, CVD-PP - 0.9860138 
pr.curve(as.numeric(changes$CLNSIG_RF[-which(is.na(changes$Probablity.of.pathogenicity_arm))]),as.numeric(changes$Probablity.of.pathogenicity_arm[-which(is.na(changes$Probablity.of.pathogenicity_arm))]))
pr.curve(as.numeric(changes$CLNSIG_RF[-which(is.na(changes$Probablity.of.pathogenicity_arm))]),as.numeric(changes$Pred[-which(is.na(changes$Probablity.of.pathogenicity_arm))]))
#Binary predictions
#ARM - 0.6231188, CVD-PP - 0.9860138
pr.curve(as.numeric(changes$CLNSIG_RF[-which(is.na(changes$cardioboost_arm_validation))]),as.numeric(changes$cardioboost_arm_validation[-which(is.na(changes$cardioboost_arm_validation))]))
pr.curve(as.numeric(changes$CLNSIG_RF[-which(is.na(changes$cardioboost_arm_validation))]),as.numeric(changes$Pred[-which(is.na(changes$cardioboost_arm_validation))]))

#Assess accuracy for linear outcome
changes.linear<-changes[,c("CLNSIG_RF","Pred","PHRED","SIFTval","PolyPhenVal",
                           "MetaSVM_score","MetaLR_score","MutationTaster_converted_rankscore",
                           "MutationAssessor_rankscore","PROVEAN_converted_rankscore",
                           "fathmm.XF_coding_score","BayesDel_addAF_score","REVEL_rankscore")]
changes.linear[changes.linear==""]<-NA
changes.linear<-changes.linear[-which(rowSums(is.na(changes.linear))>0),]

Mod_AUC_Linear.pr<-pr.curve(as.numeric(changes.linear$CLNSIG_RF),as.numeric(changes.linear$Pred), curve = T) #0.9724753 if using pred, 0.5222 if using prediction
Mod_AUC_Linear.pr #0.993205
CADD_AUC_Linear.pr<-pr.curve(as.numeric(changes.linear$CLNSIG_RF),as.numeric(changes.linear$PHRED))#0.3
CADD_AUC_Linear.pr #0.3081436 
SIFT_AUC_Linear.pr<-pr.curve(as.numeric(changes.linear$CLNSIG_RF),as.numeric(changes.linear$SIFTval))#0.3
SIFT_AUC_Linear.pr #0.9902358 
Polyphen_AUC_Linear.pr<-pr.curve(as.numeric(changes.linear$CLNSIG_RF),as.numeric(changes.linear$PolyPhenVal))#0.3
Polyphen_AUC_Linear.pr #0.8247592
MetaSVM_AUC_Linear.pr<-pr.curve(as.numeric(changes.linear$CLNSIG_RF),as.numeric(changes.linear$MetaSVM_score))#0.3
MetaSVM_AUC_Linear.pr #0.4391515
MetalR_AUC_Linear.pr<-pr.curve(as.numeric(changes.linear$CLNSIG_RF),as.numeric(changes.linear$MetaLR_score))#0.3
MetalR_AUC_Linear.pr #0.9932088
MutationAssessor_AUC_Linear.pr<-pr.curve(as.numeric(changes.linear$CLNSIG_RF),as.numeric(changes.linear$MutationAssessor_rankscore))#0.3
MutationAssessor_AUC_Linear.pr #0.993205
MutationTaster_AUC_Linear.pr<-pr.curve(as.numeric(changes.linear$CLNSIG_RF),as.numeric(changes.linear$MutationTaster_converted_rankscore))#0.3
MutationTaster_AUC_Linear.pr #0.993205 
PROVEAN_AUC_Linear.pr<-pr.curve(as.numeric(changes.linear$CLNSIG_RF),as.numeric(changes.linear$PROVEAN_converted_rankscore))#0.3
PROVEAN_AUC_Linear.pr #0.993205 
FATHMM_AUC_Linear.pr<-pr.curve(as.numeric(changes.linear$CLNSIG_RF),as.numeric(changes.linear$fathmm.XF_coding_score))#0.3
FATHMM_AUC_Linear.pr #0.993205 
BayesDel_AUC_Linear.pr<-pr.curve(as.numeric(changes.linear$CLNSIG_RF),as.numeric(changes.linear$BayesDel_addAF_score))#0.3
BayesDel_AUC_Linear.pr #0.9932667 
REVEL_AUC_Linear.pr<-pr.curve(as.numeric(changes.linear$CLNSIG_RF),as.numeric(changes.linear$REVEL_rankscore))#0.3
REVEL_AUC_Linear.pr #0.993205

#Assess accuracy for binary outcome
changes.binary<-changes[rownames(changes.linear),]
colnames(changes.binary)

Mod_AUC_Binary.pr<-pr.curve(as.numeric(changes.binary$CLNSIG_RF),as.numeric(changes.binary$Prediction), curve = T) #0.9724753 if using pred, 0.5222 if using prediction
Mod_AUC_Binary.pr
CADD_AUC_Binary.pr<-pr.curve(as.numeric(changes.binary$CLNSIG_RF),as.numeric(changes.binary$CADD_Validation))#0.5
CADD_AUC_Binary.pr
PolyPhen_AUC_Binary.pr<-pr.curve(as.numeric(changes.binary$CLNSIG_RF),as.numeric(changes.binary$PolyPhen_Validation)) #0.5228503
PolyPhen_AUC_Binary.pr
MetaSVM_AUC_Binary.pr<-pr.curve(as.numeric(changes.binary$CLNSIG_RF),as.numeric(changes.binary$MetaSVM_pred_Validation)) #0.5086307
MetaSVM_AUC_Binary.pr
MetaLR_AUC_Binary.pr<-pr.curve(as.numeric(changes.binary$CLNSIG_RF),as.numeric(changes.binary$MetaLR_pred_Validation)) #0.5081147
MetaLR_AUC_Binary.pr
MutationAssesor_AUC_Binary.pr<-pr.curve(as.numeric(changes.binary$CLNSIG_RF),as.numeric(changes.binary$MutationAssessor_pred_Validation)) #0.5386667 
MutationAssesor_AUC_Binary.pr
MutationTaster_AUC_Binary.pr<-pr.curve(as.numeric(changes.binary$CLNSIG_RF),as.numeric(changes.binary$MutationTaster_pred_Validation)) #0.5022791 
MutationTaster_AUC_Binary.pr
PROVEAN_AUC_Binary.pr<-pr.curve(as.numeric(changes.binary$CLNSIG_RF), as.numeric(changes.binary$PROVEAN_Validation))#0.5224268
PROVEAN_AUC_Binary.pr
Fathmm_AUC_Binary.pr<-pr.curve(as.numeric(changes.binary$CLNSIG_RF), as.numeric(changes.binary$fathmm.XF_coding_pred_Validation))#0.5249126
Fathmm_AUC_Binary.pr
BayesDel_AUC_Binary.pr<-pr.curve(as.numeric(changes.binary$CLNSIG_RF), as.numeric(changes.binary$BayesDel_addAF_pred_Validation))#0.5050417 
BayesDel_AUC_Binary.pr
SIFT_AUC_Binary.pr<-pr.curve(as.numeric(changes.binary$CLNSIG_RF),as.numeric(changes.binary$SIFT_Validation) ) #0.5165861
SIFT_AUC_Binary.pr
REVEL_AUC_Binary.pr<-pr.curve(as.numeric(changes.binary$CLNSIG_RF),as.numeric(changes.binary$REVEL_rankscore_Validation) ) #0.505
REVEL_AUC_Binary.pr

performance_comp<-data.frame(cbind(c(Mod_AUC_Linear.pr$auc.integral, CADD_AUC_Linear.pr$auc.integral,Polyphen_AUC_Linear.pr$auc.integral,MetaSVM_AUC_Linear.pr$auc.integral,
                                     MetalR_AUC_Linear.pr$auc.integral,MutationAssessor_AUC_Linear.pr$auc.integral,MutationTaster_AUC_Linear.pr$auc.integral,PROVEAN_AUC_Linear.pr$auc.integral,
                                     FATHMM_AUC_Linear.pr$auc.integral,BayesDel_AUC_Linear.pr$auc.integral,SIFT_AUC_Linear.pr$auc.integral, REVEL_AUC_Linear.pr$auc.integral),
                                   c(Mod_AUC_Binary.pr$auc.integral, CADD_AUC_Binary.pr$auc.integral,PolyPhen_AUC_Binary.pr$auc.integral,MetaSVM_AUC_Binary.pr$auc.integral,
                                     MetaLR_AUC_Binary.pr$auc.integral,MutationAssesor_AUC_Binary.pr$auc.integral,MutationTaster_AUC_Binary.pr$auc.integral,PROVEAN_AUC_Binary.pr$auc.integral,
                                     Fathmm_AUC_Binary.pr$auc.integral,BayesDel_AUC_Binary.pr$auc.integral,SIFT_AUC_Binary.pr$auc.integral, REVEL_AUC_Binary.pr$auc.integral)))
rownames(performance_comp)<-c("CVD-PP","CADD","PolyPhen","MetaSVM","MetalR","MutationAssessor","MutationTaster","PROVEAN","FATHMM","BayesDel","SIFT", "REVEL")

colnames(performance_comp)<-c("Linear","Binary")

write.table(performance_comp,
            file = "/Volumes/Sshare/Agenda/Ramaker/MonogenicDisease/MachineLearning_Manuscript/PR_AUC_Comparisons.txt",
            sep = "\t", col.names = T, row.names = T, quote = F)

changes<-changes[rownames(changes.linear),]
changes$CLNSIG_Simple2<-"VUS"
changes$CLNSIG_Simple2[grep("pathogenic", changes$CLNSIG_Simple, ignore.case = T)]<-"Pathogenic"
changes$CLNSIG_Simple2[grep("benign", changes$CLNSIG_Simple, ignore.case = T)]<-"Benign"
changes$Previous_CLNSIG_Simple2<-"VUS"
changes$Previous_CLNSIG_Simple2[grep("pathogenic", changes$Previous_CLNSIG_Simple, ignore.case = T)]<-"Pathogenic"
changes$Previous_CLNSIG_Simple2[grep("benign", changes$Previous_CLNSIG_Simple, ignore.case = T)]<-"Benign"


#Step 5:Predict VUS & Conflicting interpretation variants, i.e the "Discovery_Merged_Annot_RF" data
Discovery_Merged_Annot_RF$CLNSIG_RF[which(Discovery_Merged_Annot_RF$CLNSIG_RF=="Benign")]<-0
Discovery_Merged_Annot_RF$CLNSIG_RF[which(Discovery_Merged_Annot_RF$CLNSIG_RF=="Pathogenic")]<-1
table(Discovery_Merged_Annot_RF$CLNSIG_RF)
Discovery_Merged_Annot_RF$Pred<-predict(Top_RF, newdata = Discovery_Merged_Annot_RF[,c("PHRED","gnomAD_AF","MaxEntScan_diff","DistToPathVar","Expression")])
Discovery_Merged_Annot_RF$Prediction<-0
Discovery_Merged_Annot_RF$Prediction[which(Discovery_Merged_Annot_RF$Pred>OptCutoff)]<-1
table(Discovery_Merged_Annot_RF$Prediction)
summary(Discovery_Merged_Annot_RF$Pred)


rm(Truth,trees,tree,seeds,seed,Pred,path_testIndexes,path_fold,MSE,k,j,idx,i,ben_testIndexes,ben_fold,
   trainPath,trainBen,train,testPath,testBen,test,SIFT_AUC_Binary.pr,SIFT_AUC_Linear.pr,ROC_AUC_seeds,
   ROC_AUC_k,RF_seeds,RF_k,RF,REVEL_AUC_Binary.pr,REVEL_AUC_Linear.pr,PROVEAN_AUC_Binary.pr,PROVEAN_AUC_Linear.pr,
   Pred_seeds,Pred_k,Polyphen_AUC_Linear.pr,PolyPhen_AUC_Binary.pr,path,OptCutoff_k,OptCutoff_seeds,MutationAssesor_AUC_Binary.pr,
   MutationAssessor_AUC_Linear.pr,MetaLR_AUC_Binary.pr,MetalR_AUC_Linear.pr,MetaSVM_AUC_Binary.pr,MetaSVM_AUC_Linear.pr,
   MutationTaster_AUC_Binary.pr,MutationTaster_AUC_Linear.pr,MSE_seeds,MSE_k,Mod_AUC_Binary.pr,Mod_AUC_Linear.pr,Importance_seeds,
   Importance_k,importance_scores,Fathmm_AUC_Binary.pr,FATHMM_AUC_Linear.pr,changes.binary,changes.linear,
   cardioboost_arm,cardioboost_cm,CADD_AUC_Binary.pr,CADD_AUC_Linear.pr,ben,BayesDel_AUC_Binary.pr,BayesDel_AUC_Linear.pr,)

#save.image("/Volumes/Sshare/Agenda/Ramaker/MonogenicDisease/Clinvar_RandomForest_Model_Regression_06282023_simpleoutput.RData")
load("/Volumes/Sshare/Agenda/Ramaker/MonogenicDisease/Clinvar_RandomForest_Model_Regression_06282023_simpleoutput.RData")

#Figures for manuscript (need to fix these):

#Bar plot of number of variants per gene, colored by clinical significance
#with a broken y axis because of TTN - taken from http://sickel.net/blogg/?p=688
#Redo with 3 colors and groups 
Gene_DF2<-data.frame(table(Merged_Annot_Master$GENE,Merged_Annot_Master$CLNSIG_Simple2))
Gene_DF2<-Gene_DF2[order(Gene_DF2$Var1, Gene_DF2$Var2),]
tmp<-data.frame(Gene_DF2[1:3,3],row.names = Gene_DF2$Var2[1:3])
colnames(tmp)[1]<-as.character(Gene_DF2$Var1[1])
for(i in seq(4, 141, by = 3)){
  tmp<-cbind(tmp,Gene_DF2[i:(i+2),3])
  colnames(tmp)[ncol(tmp)]<-as.character(Gene_DF2$Var1[i])
}
tmp<-tmp[,order(colSums(tmp), decreasing = T)]
cols<-c("#8ECAE6","#9B2226","#FFB703")
#Broken bar plot
lower=c(0,5000)
upper=c(8500,18000)
y_outer=21

lowspan=c(0,11)
topspan=c(lowspan[2]+1,21)

ylabel="Number of Variants"
xlabel="Gene"
legendtext=colnames(tmp)

cnvrt.coords <-function(x,y=NULL){
  # Stolen from the teachingDemos library, simplified for this use case
  xy <- xy.coords(x,y, recycle=TRUE)
  cusr <- par('usr')
  cplt <- par('plt')	
  plt <- list()
  plt$x <- (xy$x-cusr[1])/(cusr[2]-cusr[1])
  plt$y <- (xy$y-cusr[3])/(cusr[4]-cusr[3])
  fig <- list()
  fig$x <- plt$x*(cplt[2]-cplt[1])+cplt[1]
  fig$y <- plt$y*(cplt[4]-cplt[3])+cplt[3]
  return( list(fig=fig) )
}

subplot <- function(fun, x, y=NULL){
  # Stolen from the teachingDemos library, simplified for this use case
  old.par <- par(no.readonly=TRUE)
  on.exit(par(old.par))
  xy <- xy.coords(x,y)
  xy <- cnvrt.coords(xy)$fig
  par(plt=c(xy$x,xy$y), new=TRUE)
  fun
  tmp.par <- par(no.readonly=TRUE)
  return(invisible(tmp.par))
}

#Make Figure 2A
tiff("/Volumes/Sshare/Agenda/Ramaker/MonogenicDisease/MachineLearning_Manuscript/GenePlot_broken.tiff", width = 800, height = 700)
plot(c(0,1),c(0,y_outer),type='n',axes=FALSE,ylab=ylabel,xlab='',lwd=7)

subplot(barplot(as.matrix(tmp),col=cols,ylim=lower,xpd=FALSE,las=3),x=c(0,1),y=lowspan)

subplot(barplot(
  as.matrix(tmp),
  col=cols,
  ylim=upper,
  xpd=FALSE,
  names.arg=vector(mode="character",length=ncol(tmp))), 
  x=c(0,1),
  y=topspan)

lowertop=lowspan[2]+0.1     # Where to end the lower axis
breakheight=0.5   # Height of the break
upperbot=lowertop+breakheight # Where to start the upper axes
markerheight=0.4 # Heightdifference for the break markers
markerwidth=.04  # With of the break markers
abline(h = 0, col = "black")
lines(c(0,0),c(1,lowertop))
lines(c(markerwidth/-2,markerwidth/2),c(lowertop-        
                                          markerheight/2,lowertop+markerheight/2))
lines(c(0,0),c(upperbot,14))
lines(c(markerwidth/-2,markerwidth/2),c(upperbot-    
                                          markerheight/2,upperbot+markerheight/2))

dev.off()

#Make Figure 2B
#histogram of the proportion of clinvar variants in each gene that are VUS
percent_VUS<-data.frame(tmp[3,]/colSums(tmp))
tiff("/Volumes/Sshare/Agenda/Ramaker/MonogenicDisease/MachineLearning_Manuscript/ProportionVUS_Genes.tiff", width = 800, height = 700)
barplot(as.matrix(percent_VUS)[order(percent_VUS[1,], decreasing = T)],
        ylim = c(0,0.6), names.arg = colnames(percent_VUS)[order(percent_VUS[1,], decreasing = T)],
        las = 2, col = "black", ylab = "Proportion VUS")
dev.off()

#Make Figure 3 
#Boxplot of changed (validation) variant CVD-PP prediction scores by 2022 clinvar classification
tiff("/Volumes/Sshare/Agenda/Ramaker/MonogenicDisease/MachineLearning_Manuscript/CVDPP_Validation_Score_Boxplot.tiff", width = 580)
boxplot(changes$Pred~changes$CLNSIG_Simple,notch = T , outline = F, 
        ylab = "CVD-PP Score", xlab = "ClinVar 2022 Classification")
stripchart(as.numeric(changes$Pred)~changes$CLNSIG_Simple,
           method = "jitter",
           pch = 19,
           col = c("black","black", "red","red"),
           vertical = TRUE,
           add = TRUE)
dev.off()

#Make Figure 4
#Histogram of VUS scores with naive bayes cut off indicated in red dashed line
tiff("/Volumes/Sshare/Agenda/Ramaker/MonogenicDisease/MachineLearning_Manuscript/VUS_CVDPP_PredictionScores_Histogram.tiff", width = 580)
hist(Discovery_Merged_Annot_RF$Pred, xlab = "RF Prediction", 
     main = "Histogram of VUS CVD-PP prediction scores", breaks = 100, col = "black")
abline(v = OptCutoff, col = "red", lwd = 4, lty = 3)
dev.off()
