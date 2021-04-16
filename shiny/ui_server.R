

# Packages
## shinypackages
library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(shinyBS)
library(shinyjs)

## Array Statistics pkgs
library(limma)
library(affy)
library(annotate)
library(mouse4302.db)  # load chip-specific annotation

# Data frame combination pkg
library(dplyr)

## Visualisation pkgs
library(ggplot2)
library(affyPLM)  # Chip pseudo images
library(scatterplot3d)
library(pheatmap)  # Heat map


#load the example data
load("/home/ninomoriaty/0Exam/FGT-ICA/shiny/Affy.RData")


# panels
## Quality Control
qc <- tabItem("qc",
              
              fluidRow(
                width = 12,
                box(
                  width = 12,
                  status = "info",
                  solidHeader = TRUE,
                  title = div(strong("Quality Control")),
                  helpText("Tips: This panel is used to show the quality control plots. You can select the plots that interest you and present the figure in customised size."),
                  br()
                )
              ),
              
              sidebarLayout(
                
                sidebarPanel(
                  radioButtons("qcPlotOption", "Plot Options:", 
                               choices=c("Histogram(Curve)","Boxplot","MA plot")), 
                  sliderInput("heightx",
                              "Figure height (px):",
                              min = 0,
                              max = 900,
                              value = 400),
                  sliderInput("widthx",
                              "Figure width (px):",
                              min = 0,
                              max = 1600,
                              value = 400)
                  ),
                
                # Show a plot of the generated distribution
                mainPanel(
                  uiOutput("qcPlot.ui")
                )
              )
)

# Differential Expression
de <- tabItem("de",
              
              # Figure Section
              fluidRow(
                width = 12,
                box(
                  width = 12,
                  status = "info",
                  solidHeader = TRUE,
                  title = div(strong("Differential Expression")),
                  h3("Figures Section"),
                  helpText("Tips: This panel is used to show the cluster dendrogram, PCA figure. You can select the plots that interest you and present the figure in customised size."),
                  br()
                )
              ),
              
              sidebarLayout(
                
                sidebarPanel(
                  radioButtons("dePlotOption", "Plot Options:", 
                               choices=c("Cluster Dendrogram", "PCA", "Venn diagram", "Heatmap")), 
                  sliderInput("heightdex",
                              "Figure height (px):",
                              min = 0,
                              max = 900,
                              value = 400),
                  sliderInput("widthdex",
                              "Figure width (px):",
                              min = 0,
                              max = 1600,
                              value = 400)
                ),
                
                # Show a plot of the generated distribution
                mainPanel(
                  uiOutput("dePlot.ui")
                )
              ),
              
              # Table section
              fluidRow(
                width = 12,
                box(
                  width = 12,
                  status = "info",
                  solidHeader = TRUE,
                  h3("Table Section"),
                  helpText("Tips: This panel is used to show the data table for differential genes from Limma and fold change tables. You can search elements and sort columns through the data table."),
                  br()
                )
              ),
              
              sidebarLayout(
                
                sidebarPanel(
                  radioButtons("deTableOption", "Table Options:", 
                               choices=c("Limma", "Fold change")), 

                ),
                
                # Show a plot of the generated distribution
                mainPanel(
                  
                    uiOutput("DeStatTable.ui")
                  
                )
              )
              
              
              
              
              
)

# Functional Enrichment
fe <- tabItem("fe",
              
              fluidRow(
                box(
                  width = 12,
                  status = "info",
                  solidHeader = TRUE,
                  title = div(strong("Functional Enrichment")),
                  helpText("This table shows the summary of the functional enrichment results and differential expression genes."),
                  br()
                )
              ),
              
              sidebarLayout(
                
                sidebarPanel(

                  
                ),
                
                
                mainPanel(
                  
                  dataTableOutput("feTable")
                  
                )
              )

)

# ui.R The ui interface arrangement
ui <- fluidPage(
  # Title panel
  # Application title
  titlePanel("AffyWorkflow R shiny Prototype"),
  
  mainPanel(
    tabsetPanel(
      tabPanel("Quality Control", qc),
      tabPanel("Differential Expression", de),
      tabPanel("Functional Enrichment", fe)
    ))
  
)

# server.R The server script to response
server <- function(input, output,session) {
  
  # Quality Control
  output$qcPlot.ui <- renderUI({
    plotOutput("qcPlot", width = input$widthx, height = input$heightx)
  })
  
  output$qcPlot <- renderPlot({
    # Choose the plot option
    if (input$qcPlotOption == "Histogram(Curve)") {
      hist(mydata, xlab = 'Log intensity', ylab = 'Density')
      legend("topright", rownames(dfPhenoData), col = rownames(dfPhenoData))
    } else if (input$qcPlotOption == "Boxplot") {
      boxplot(values, col = lsColours, xaxt = "n")  
      text(seq_along(mydata), par("usr")[3] - 0.5, labels = rownames(dfPhenoData), 
           srt = 20, adj = 1, xpd = TRUE, cex = 0.8)
    } else if (input$qcPlotOption =="MA plot") {
      mva.pairs(values)
    }
  }, execOnResize = F)
  
  # Differential Expression
  output$dePlot.ui <- renderUI({
    plotOutput("dePlot", width = input$widthdex, height = input$heightdex)
  })
  
  output$DeStatTable.ui <- renderUI({
    dataTableOutput('DeStatTable')
  })
  
  output$dePlot <- renderPlot({
    
    # Choose the plot option
    if (input$dePlotOption == "Cluster Dendrogram") {
      plot(hc)
    } else if (input$dePlotOption == "PCA") {
      s3d<-scatterplot3d(pca$x[,1:3], pch=19, 
                         color = rainbow(length(dfPhenoData$sample)))
      s3d.coords <- s3d$xyz.convert(pca$x[,1:3])
      text(s3d.coords$x, s3d.coords$y, labels = colnames(values), 
           pos = 3,offset = 0.5, cex = 0.8)
    } else if (input$dePlotOption == "Venn diagram") {
      vennDiagram(clas)
    } else if (input$dePlotOption == "Heatmap"){
      pheatmap(plotme,scale="row")
    }
    
  }, execOnResize = F)
  
  output$DeStatTable <- renderDataTable({
   if (input$deTableOption == "Fold change"){
     all.data
   } else if (input$deTableOption == "Limma") {
     myresults_LFC
   } 
  },
    options = list(
      pageLength = 5,
      lengthMenu = c(5, 10, 15, 20)
    )
  )
  
  
  # Functional Enrichment
  output$feTable <- renderDataTable({
    dfFinal
  }, 
  options = list(
    pageLength = 5,
    lengthMenu = c(5, 10, 15, 20)))
  
  observeEvent(input$refresh, {
    session$invalidate
  })
  
  
}

  
# Run the application 
shinyApp(ui = ui, server = server)
