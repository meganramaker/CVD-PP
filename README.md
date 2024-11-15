# CVD-PP

A shiny app for accessing CVD-PP ClinVar VUS scores and predictions.

The following packages are required to run the CVD-PP shiny application and can be installed and loaded using these lines of code:

if (!require("DT")) install.packages('DT')

if (!require("htmltools")) install.packages('htmltools')

if (!require("data.table")) install.packages('data.table')

library(DT,htmltools,data.table)

To run this app use the following line of code: 

shiny::runGitHub(username = 'meganramaker',repo = 'CVD-PP',ref = 'main')
