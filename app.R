library(shiny)
library(DT)
library(data.table)
VUS<-readRDS("Clinvar_VUS_RF_Predictions.rds")
ui <- basicPage(
  h2("CVD-PP VUS Table"),
  DT::dataTableOutput("mytable")
)

server <- function(input, output) {
  output$mytable = DT::renderDataTable({
    VUS
  })
}

shinyApp(ui, server)
