#Script for calculating mean gtex exon expression in Cardiac tissue for use in CVD-PP
#Megan Ramaker
#10/4/2023
library(data.table)
gtex_samples<-data.frame(fread("/Volumes/Sshare/Agenda/Ramaker/MonogenicDisease/GTEx_Data_V6_Annotations_SampleAttributesDS.txt", sep = "\t", stringsAsFactors = F))
gtex_samples<-gtex_samples[which(gtex_samples$SMTS=="Heart"),]

gtex<-data.frame(fread("/Volumes/Sshare/Agenda/Ramaker/MonogenicDisease/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_exon_reads.txt", sep = "\t",stringsAsFactors = F, select = c("Id", gtex_samples$SAMPID)))
rownames(gtex)<-gtex[,"Id"]
gtex<-gtex[,-1]
gtex<-as.matrix(gtex)
#Divide each exon by (sum of reads for that sample/1000000) to get counts per million
for (col in 1:ncol(gtex)) {
  gtex[,col]<-(gtex[,col]/(sum(gtex[,col])/1000000))
}
summary(rowMeans(gtex))
hist(rowMeans(gtex))
rownames(gtex)[which.max(rowMeans(gtex))]

#Remove what isn't expressed at all
#gtex<-gtex[-which(rowSums(gtex)==0),]
summary(rowMeans(gtex[-which(rowSums(gtex)<=0.06),]))

gtex<-gtex[-which(rowSums(gtex == 0) >= 206),]
summary(rowMeans(gtex))

gtf<-data.frame(fread("~/Documents/CALERIE/Analysis_Files/Homo_sapiens.GRCh37.87.gtf", stringsAsFactors = F))
gtf<-gtf[which(gtf$V3=="exon"),]
#gtex ids are ens gene id . gene version _ exon num
gtf$Id<-gsub("\\\".*$", "",gsub(".*gene_id \\\"", "", gtf$V9))
gtf$version<-gsub("\\\".*$", "",gsub(".*gene_version \\\"", "", gtf$V9))
gtf$exon<-gsub("\\\".*$", "",gsub(".*exon_number \\\"", "", gtf$V9))
gtf<-gtf[order(gtf$V2),]
gtf$key<-paste0(gtf$Id,".",gtf$version,"_",gtf$exon)
gtf<-gtf[which(!duplicated(gtf$key)),]
rownames(gtf)<-gtf$key
idx<-intersect(rownames(gtf),rownames(gtex))
gtex<-gtex[idx,]
gtf<-gtf[idx,]
identical(rownames(gtf),rownames(gtex))

#next include median exon exp for each variant in merged annot
gtex_medians<-lapply(1:nrow(gtex), function(x) median(as.numeric(gtex[x,])))
gm<-data.frame(unlist(gtex_medians))
rownames(gm)<-rownames(gtex)
gtex_medians<-gm
rm(gm)
colnames(gtex_medians)[1]<-"median_cpm"
identical(rownames(gtf),rownames(gtex_medians))
gtf<-cbind(gtf,gtex_medians[,1])
colnames(gtf)[14]<-"gtex_median_cpm"

#Create bed file for genome build liftover
#gtex_gtf_bed<-gtf[,c(1,4,5,13)]
#gtex_gtf_bed$X..genebuild.last.updated.2013.09<-paste0("chr",gtex_gtf_bed$X..genebuild.last.updated.2013.09)
#write.table(gtex_gtf_bed, file = "/Volumes/Sshare/Agenda/Ramaker/MonogenicDisease/gtex_gtf_GRCh37.bed", sep = "\t", col.names = F, row.names = F, quote = F)
liftover<-data.frame(fread(file = "/Volumes/Sshare/Agenda/Ramaker/MonogenicDisease/hglft_genome_3b6cc_e70900.bed", sep = "\t", stringsAsFactors = F))
rownames(liftover)<-liftover$V4

idx<-intersect(rownames(gtf), rownames(liftover))
gtf<-gtf[idx,]
liftover<-liftover[idx,]
identical(rownames(liftover),rownames(gtf))
gtf<-cbind(gtf,liftover[,1:3])
colnames(gtf)[c(1,4,5,15:17)]<-c("GRCh37Chr","GRCh37start","GRCh37stop","GRCh38Chr","GRCh38start","GRCh38stop")
gtf$GRCh38Chr<-gsub("chr","",gtf$GRCh38Chr)
gtf$GRCh38Chr<-gsub("_.*$","",gtf$GRCh38Chr)
gtf<-gtf[which(gtf$GRCh37Chr==gtf$GRCh38Chr),]
bed<-gtf[,c(15:17,9,10:12,14)]
bed$GRCh38Chr<-paste0("chr",bed$GRCh38Chr)
write.table(bed, file = "/Volumes/Sshare/Agenda/Ramaker/MonogenicDisease/GTEx_Cardiac_MedianExp_Liftover.bed", sep = "\t", col.names = F,row.names = F, quote = F)
write.table(gtf, file = "/Volumes/Sshare/Agenda/Ramaker/MonogenicDisease/GTEx_Cardiac_MedianExp_Liftover.txt", sep = "\t", col.names = T,row.names = F, quote = F)



tmp<-Merged_Annot[,c("CHROM","POS","POS","Key")]
tmp$CHROM<-paste0("chr", tmp$CHROM)
tmp$POS<-tmp$POS-1
write.table(tmp, file = "/Volumes/Sshare/Agenda/Ramaker/MonogenicDisease/Merged_Annot.bed",col.names = F, row.names = F, quote = F, sep = "\t")

tmp<-gtf[,c(1,4,5,2,3,6:14)]
# tmp<-tmp[which(tmp$V7=="+"),]
# colnames(tmp)<-1:ncol(tmp)
# tmp1<-gtf[,c(1,5,4,2,3,6:14)]
# tmp1<-tmp1[which(tmp1$V7=="-"),]
# colnames(tmp1)<-1:ncol(tmp1)
# tmp2<-rbind(tmp,tmp1)

tmp$X..genebuild.last.updated.2013.09<-paste0("chr", tmp$X..genebuild.last.updated.2013.09)
write.table(tmp, file = "/Volumes/Sshare/Agenda/Ramaker/MonogenicDisease/GTF.bed",col.names = F, row.names = F, quote = F, sep = "\t")