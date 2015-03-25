library(shiny)
load('dssTC.Rdata')
library(DESeq2)
library(ggplot2)
library(DT)

shinyUI(fluidPage(

  # Application title
  singleton(tags$script(type="text/javascript", src="addons.js")),
  
  fluidRow(
      column(
          4, 
          titlePanel("DEview"),
          conditionalPanel(
              condition = 'input.plot',
              plotOutput("distPlot")
          ),
          selectInput(
              inputId = 'type', 
              label = 'Condition', 
              choices = c('TS', colnames(colData(ddsTC))[1:2])
          ),
          conditionalPanel(
              condition = 'input.type != "TS"',
              selectInput('p1', 'var 1', levels(colData(ddsTC)[['strain']])),
              selectInput('p2', 'var 2', levels(colData(ddsTC)[['strain']])),
              
              checkboxInput('filter', 'Filter'),
              conditionalPanel(
                  condition = 'input.filter',
                  selectInput('which', 'Which value filter on', colnames(colData(ddsTC))[1:2]),
                  selectInput('what', 'What to filter on', levels(colData(ddsTC)[['strain']]))
              )
          ),
          actionButton(inputId = 'apply', label = 'Apply filters and conditions')
      ),
      column(
          8,
          tags$br(),
          DT::dataTableOutput("data"),
          downloadButton('downloadData', label = "Get result table as CSV", class = NULL),
          downloadButton('downloadFigure', label = "Get figure as PDF", class = NULL),
          downloadButton('downloadDataFlt', label = "Get filtered results as CSV", class = NULL),
          tags$br(),
          tags$strong('Numeric column definitions:'),
          verbatimTextOutput('info'),
          div( 
              class='', id='debug', tags$hr(),
              'Debug console: ', tags$br(), tags$textarea(id='debug_cmd', rows=4, style='width:88%'),
              actionButton('debug_submit', 'Submit'), verbatimTextOutput("debug_out")
          ),
          DT::dataTableOutput("data2")

      )
      
  )

  # Sidebar with a slider input for number of bins
    
   

))
