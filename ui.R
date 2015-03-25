library(shiny)
load('dssTC.Rdata')
library(DESeq2)
library(ggplot2)
library(DT)

shinyUI(fluidPage(
    
    # Application title
    singleton(tags$script(type="text/javascript", src="addons.js")),
    singleton(tags$link(type="text/css", rel="stylesheet", href="style.css")),
    
    fluidRow(
        column(
            4, 
            titlePanel("DEview"),
            conditionalPanel(
                condition = 'input.plot',
                plotOutput("distPlot")
            ),
            
            selectInput(
                inputId = 'test', 
                label = 'Statistical test', 
                choices = list('Wald (pairwise comparison)'="Wald", "Likelihood-ratio test (model comparison)"='LRT')
            ),
            
            
            conditionalPanel(
                condition = 'input.test == "Wald"',
                selectInput(
                    inputId = 'type', 
                    label = 'Condition', 
                    choices = c(colnames(colData(ddsTC))[1:2]),
                    selected = 'stage'
                ),
                selectInput('p1', 'var 1', levels(colData(ddsTC)[['stage']]), selected=head(levels(colData(ddsTC)[['stage']]), 1)),
                selectInput('p2', 'var 2', levels(colData(ddsTC)[['stage']]), selected=tail(levels(colData(ddsTC)[['stage']]), 1))
            ),
            
            conditionalPanel(
                condition = 'input.test == "LRT"',
                textInput('m1', 'Model formula', '~ stage'),
                textInput('m0', 'Reduced formula to compare against', '~ 1')
                
            ),
            
            checkboxInput('filter', 'Filter condition (recommended to filter unused condition(s) for Wald test)', value = TRUE),
            conditionalPanel(
                condition = 'input.filter',
                conditionalPanel(
                    condition = 'input.test == "LRT"', 
                    selectInput('which', 'Which value filter on', colnames(colData(ddsTC))[1:2], selected = 'strain')
                ),
                radioButtons('what', 'Use following [strain]:', levels(colData(ddsTC)[['strain']]), selected = 'N2', inline=FALSE)
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
