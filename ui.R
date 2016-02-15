library(shiny)
library(DESeq2)
library(ggplot2)
library(DT)

shinyUI(navbarPage(
    'DEview', id = 'nvpage',


    tabPanel('Settings',
             singleton(tags$script(type="text/javascript", src="addons.js")),
             singleton(tags$link(type="text/css", rel="stylesheet", href="style.css")),
             
    # Application title

######## Settings ###############


        tabsetPanel(type = 'tabs', position = 'above', tabPanel('Statistical test',
            
            
            div(class='row', div(class='col-md-4',
                   selectInput('SEdata', 'Dataset', choices = dir('data'), width = '400px' )
            ), class='row', div(class='col-md-3',
                   shiny::fileInput('newfile', 'Add dataset')
            ), div(class='col-md-1',
                   tags$br(), shiny::downloadButton('dlfile', 'Download', class='btn btn-success btn-sm')
            ), div(class='col-md-1',
                   tags$br(), shiny::actionButton('rmfile', 'Remove', icon = icon('remove'), class='btn btn-danger btn-sm', position='middle')
            )),
            tags$hr(),
            
            conditionalPanel(condition = 'input.advstat', selectInput(
                inputId = 'test',
                label = 'Statistical test', 
                choices = list('Wald (pairwise comparison)'="Wald", "Likelihood-ratio test (model comparison)"='LRT')
            )), 
            conditionalPanel(
                condition = 'input.test == "Wald"',
                div(class='row', div(class='col-md-2',
                    radioButtons(
                        inputId = 'type', 
                        label = 'Compare', 
                        choices = ""
                    )
                ), div(class='col-md-2',
                       radioButtons('p1', paste0('[',']'), "")
                ), div(class='col-md-2',
                       radioButtons('p2', 'versus', "")
                ))
            ),
            conditionalPanel(
                condition = 'input.test == "LRT"',
                textInput('m1', 'Model formula', paste0('~ ', "strain" )),
                textInput('m0', 'Reduced formula to compare against', '~ 1')
                
            ),
            conditionalPanel(condition = 'input.test == "LRT"',
                checkboxInput('filter', 'Filter condition (recommended to filter unused condition(s) for Wald test)', value = TRUE)
            ),
            conditionalPanel(
                condition = 'input.filter',
                conditionalPanel(
                    condition = 'input.test == "LRT"', 
                    selectInput('which', 'Which value filter on', "")
                ),
                radioButtons('what', paste0('Use following [',']'), "", inline=FALSE)
            ),
            actionButton(inputId = 'apply', label = 'Apply settings', class='btn-success')
            
            ), tabPanel(
                'Outputs/plot optins',
                checkboxGroupInput('add', 'Add to CSV', list('Raw counts'='R', 'Normalized counts'='NR', 'RPKM'='RPKM'), inline = TRUE),
                
                downloadButton('downloadData', label = "Get result table as CSV", class = NULL),
                downloadButton('downloadDataFlt', label = "Get filtered results as CSV", class = NULL),
                
                
                tags$hr(),
                
                radioButtons('plotValues', 'Plot: ', list('All data points'='a', 'Filtered'='f', 'Select'='s'), inline=TRUE),
                conditionalPanel(condition = 'input.plotValues == "s"', div(class='row', 
                    div(class='col-md-4',
                        checkboxGroupInput('plotValues_f1', "", "")
                    ), div(class='col-md-4',
                        checkboxGroupInput('plotValues_f2', "", "")
                ))),
                
                radioButtons('plotType', 'Plot type: ', list('Norm counts'='N', 'RPKM'='fpkm', 'RPM'='fpm'), inline=TRUE),
                radioButtons('plotScale', 'Plot scale: ', list('log10'='log10', 'log2'='log2', 'linear'='N'), inline=TRUE),
                downloadButton('downloadFigure', label = "Get figure as PDF", class = NULL),
                
                tags$hr(),
                'Advanced',
                checkboxInput('debug', 'Debug console'),
                checkboxInput('advstat', 'Enable advanced stats options'),
                textInput('annosoeurce', 'Annotation link', 'http://www.wormbase.org/species/c_elegans/gene/')
            
            ), tabPanel(
                "Data and design",
                DT::dataTableOutput("design"),
                checkboxInput('desall', 'Show all')
                #,
                
                #selectInput('dataset', 'Dataset', dir('data'))
            )
        )
            ######## EIS ###############
    ), tabPanel('Results',
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
            
    ),
    tabPanel('Plot',
        plotOutput("distPlot")
    )
        
    
    
    # Sidebar with a slider input for number of bins
    
    
    
))
