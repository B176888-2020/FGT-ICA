

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
              
              fluidRow(
                box(
                  width = 12,
                  status = "info",
                  solidHeader = TRUE,
                  title = div(strong("Differential Expression")),
                  p(".",
                    style = "font-size:18px; font-weight:500;line-height:40px;"),
                  br()
                )
              ),
              
              sidebarLayout(
                sidebarPanel(
                  sliderInput("font_row",
                              "Font size row:",
                              min = 6,
                              max = 14,
                              value = 10),
                  sliderInput("font_col",
                              "Font size col:",
                              min = 6,
                              max = 14,
                              value = 10)
                ),
                
                # Show a plot of the generated distribution
                mainPanel(
                  plotOutput("dePlot", width = "100%")
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
                  p(".",
                    style = "font-size:18px; font-weight:500;line-height:40px;"),
                  br()
                )
              ),
              
              sidebarLayout(
                
                sidebarPanel(
                  sliderInput("font_row",
                              "Font size row:",
                              min = 6,
                              max = 14,
                              value = 10),
                  sliderInput("font_col",
                              "Font size col:",
                              min = 6,
                              max = 14,
                              value = 10)
                ),
                
                # Show a plot of the generated distribution
                mainPanel(
                  plotOutput("fePlot", width = "100%")
                )
              )
)

# The ui interface arrangement
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

# Define server logic required to draw a histogram
server <- function(input, output,session) {
  
  # UI Part
  output$qcPlot.ui <- renderUI({
    plotOutput("qcPlot", width = input$widthx, height = input$heightx)
  })
  
  # Plot Part
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
  
  output$dePlot <- renderPlot({
    
    hist(mydata, xlab = 'Log intensity', ylab = 'Density')
    
  }, execOnResize = F)
  
  output$fePlot <- renderPlot({
    
    hist(mydata, xlab = 'Log intensity', ylab = 'Density')
    
  }, execOnResize = F)
  
  observeEvent(input$refresh, {
    session$invalidate
  })
  
  
}

  
# Run the application 
shinyApp(ui = ui, server = server)
