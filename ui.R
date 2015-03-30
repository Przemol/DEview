library(shiny)
load('SE.Rdata')
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
            tabsetPanel(type = 'tabs', position = 'above', tabPanel('Statistical test',
            conditionalPanel(condition = 'input.advstat', selectInput(
                inputId = 'test',
                label = 'Statistical test', 
                choices = list('Wald (pairwise comparison)'="Wald", "Likelihood-ratio test (model comparison)"='LRT', "Default results from DEseq file"='asis')
            )), 
            conditionalPanel(
                condition = 'input.test == "Wald"',
                div(class='row', div(class='col-md-3',
                    radioButtons(
                        inputId = 'type', 
                        label = 'Compare', 
                        choices = colnames(colData(SE))[1:2],
                        selected = colnames(colData(SE))[2]
                    )
                ), div(class='col-md-3',
                       radioButtons('p1', paste0('[',colnames(colData(SE))[2],']'), levels(colData(SE)[[2]]), selected=head(levels(colData(SE)[[2]]), 1))
                ), div(class='col-md-3',
                       radioButtons('p2', 'versus', levels(colData(SE)[[2]]), selected=tail(levels(colData(SE)[[2]]), 1))
                ))
            ),
            
            conditionalPanel(
                condition = 'input.test == "LRT"',
                textInput('m1', 'Model formula', paste0('~ ', colnames(colData(SE))[2] )),
                textInput('m0', 'Reduced formula to compare against', '~ 1')
                
            ),
            conditionalPanel(condition = 'input.test == "LRT"',
                checkboxInput('filter', 'Filter condition (recommended to filter unused condition(s) for Wald test)', value = TRUE)
            ),
            conditionalPanel(
                condition = 'input.filter',
                conditionalPanel(
                    condition = 'input.test == "LRT"', 
                    selectInput('which', 'Which value filter on', colnames(colData(SE))[1:2], selected = colnames(colData(SE))[1])
                ),
                radioButtons('what', paste0('Use following [',colnames(colData(SE))[1],']'), levels(colData(SE)[[1]]), selected = 'N2', inline=FALSE)
            ),
            
            actionButton(inputId = 'apply', label = 'Apply settings', class='btn-success')
            ), tabPanel(
                'Outputs/plot optins',
                checkboxGroupInput('add', 'Add to CSV', list('Raw counts'='R', 'Normalized counts'='NR', 'Robust RPKM'='RPKM'), inline = TRUE),
                
                downloadButton('downloadData', label = "Get result table as CSV", class = NULL),
                downloadButton('downloadDataFlt', label = "Get filtered results as CSV", class = NULL),
                
                
                tags$hr(),
                
                radioButtons('plotValues', 'Plot: ', list('All data points'='a', 'Filtered'='f', 'Select'='s'), inline=TRUE),
                conditionalPanel(condition = 'input.plotValues == "s"', div(class='row', 
                    div(class='col-md-4',
                        checkboxGroupInput('plotValues_f1', colnames(colData(SE))[1], levels(colData(SE)[[1]]), selected=levels(colData(SE)[[1]]))
                    ), div(class='col-md-4',
                        checkboxGroupInput('plotValues_f2', colnames(colData(SE))[2], levels(colData(SE)[[2]]), selected=levels(colData(SE)[[2]]))
                ))),
                
                radioButtons('plotType', 'Plot type: ', list('Norm counts'='N', 'RPKM'='fpkm', 'RPM'='fpm'), inline=TRUE),
                radioButtons('plotScale', 'Plot scale: ', list('log10'='log10', 'log2'='log2', 'linear'='N'), inline=TRUE),
                downloadButton('downloadFigure', label = "Get figure as PDF", class = NULL),
                
                tags$hr(),
                'Advanced',
                checkboxInput('debug', 'Debug console'),
                checkboxInput('advstat', 'Enable advanced stats options')
            
            ))
        ),
        column(
            8,
            tags$br(),
            DT::dataTableOutput("data"),

            tags$strong('Numeric column definitions:'),
            verbatimTextOutput('info'),
            conditionalPanel(
                condition = 'input.debug',
                div( 
                    class='', id='debug', tags$hr(),
                    'Debug console: ', tags$br(), tags$textarea(id='debug_cmd', rows=4, style='width:88%'),
                    actionButton('debug_submit', 'Submit'), verbatimTextOutput("debug_out")
                )
                #DT::dataTableOutput("data2")
            )
            
        )
        
    )
    
    # Sidebar with a slider input for number of bins
    
    
    
))
