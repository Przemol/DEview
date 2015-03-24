library(shiny)
load('dssTC.Rdata')
library(DESeq2)
library(ggplot2)
library(DT)


shinyServer(function(input, output, session) {

    rv <- reactiveValues(info='start', dds=NULL, plot=NULL, table=NULL)
    
    output$distPlot <- renderPlot({
        if(is.null(input$plot)) return()
        gene <- input$plot
        data <- plotCounts(rv$dds, gene, intgroup=c("stage","strain"), returnData=TRUE)
        info <- sapply(mcols(rv$dds[gene,])[,1:2], as.character)
        #bm <- round( mcols(rv$dds[gene,])[,3], 2) 
        g <- ggplot(data, aes(x=stage, y=count, color=strain, group=strain)) + 
            geom_point() + scale_y_log10() +
            stat_smooth(se=FALSE,method="loess", size = 1.3) +
            ggtitle(paste0(
                info[1],  ' (', gene, ')\n ', info[2] #, '\n FDR=', res[gene,]$padj, '; BM=', bm
            ))
        rv$plot <- g
        g
        
        
    })
    
    observe({
        if(input$type == 'TS') return()
        cc <- rev(levels(colData(ddsTC)[[input$type]]))
        updateSelectInput(session, 'p1', choices = cc )
        updateSelectInput(session, 'p2', choices = cc, selected = cc[2] )
        updateSelectInput(session, 'which', choices = colnames(colData(ddsTC))[1:2][colnames(colData(ddsTC))[1:2] != input$type] )
        updateSelectInput(session, 'what', choices = levels(colData(ddsTC)[[input$which]]) )
        
    }, priority = 10)

    
    getResultTable <- reactive({
   
            progress <- shiny::Progress$new(session, min=0, max=3)
            on.exit(progress$close())
            
            progress$set(message = 'Calculating statistical tests', value = 1,
                         detail = 'This may take a while...')
            
            if(input$type == 'TS') {
                res <- results(ddsTC)
                res$log2FoldChange <- res$lfcSE <- ''
                rv$info <- paste0(
                    'Log2 Fold Change and its standard error available for paired tests only!\n', 
                    paste0(elementMetadata(res)[-(2:3),2], collapse='\n')
                )
                rv$dds <- ddsTC
            } else {
                if(input$filter) {
                    dds <- ddsTC[,colData(ddsTC)[[input$which]] == input$what]
                    design(dds)  <- as.formula(paste('~', input$type))
                    dds <- DESeq(dds)
                    rv$dds <- dds
                } else {
                    dds <- ddsTC
                    rv$dds <- ddsTC
                }
                progress$set(message = 'Calculating tests results', value = 2)
                res <- results(dds, contrast=c(input$type, input$p1, input$p2), test="Wald")
                rv$info <- paste0(elementMetadata(res)[,2], collapse='\n')
            }
            
            progress$set(message = 'Building result table', value = 3,
                         detail = 'This may take a while...')
            
            tab <- cbind(
                ID=rownames(res), 
                gene=mcols(ddsTC)$geneName,
                info=mcols(ddsTC)$desc,
                as.data.frame(res), 
                Plot=''
            )
            return(tab)
   
    })
    
    output$info <- renderText({
        input$apply
        rv$table <- isolate( getResultTable() )
        return(rv$info)
    })

    output$data = DT::renderDataTable({
        action = session$registerDataObj('iris', rv$table, shiny:::dataTablesJSON)
        sketch = htmltools::withTags(table(
            tableHeader(rv$table),
            tableFooter(rep('', ncol(rv$table)))
        ))
        callback = JS("
            $('table tfoot th').slice(0,2).each( function () {
                var title = $('table thead th').eq( $(this).index() ).text();
                var width = $('table thead th').eq( $(this).index() ).width()+25;
                $(this).html( '<input type=\"text\" placeholder=\"'+title+'\" style=\"width:'+width+'px;\" />' );
            } );

            $('table tfoot th').slice(3,9).each( function () {
                var title = $('table thead th').eq( $(this).index() ).text();
                var width = 50;
                $(this).html( '<input class=\"min\" type=\"text\" placeholder=\"'+'min'+'\" style=\"width:'+width+'px;\" /><br />' +
                              '<input class=\"max\" type=\"text\" placeholder=\"'+'max'+'\" style=\"width:'+width+'px;\" />' );
            } );
            $('table tfoot th').css('text-align', 'right');
            $('table tfoot th').css('padding', 5);

            table.columns().eq( 0 ).each( function ( colIdx ) {
                $( 'input', table.column( colIdx ).footer() ).on( 'keyup change', function () {

                    if(this.className == 'min') {
                        var flt = $(this).val() + ',' + $(this).siblings('.max').val();
                        table.column( colIdx ).search( flt ).draw();
                    } else if(this.className == 'max') {
                        var flt = $(this).siblings('.min').val() + ',' + $(this).val();
                        table.column( colIdx ).search( flt ).draw();
                    } else {
                        table.column( colIdx ).search( this.value ).draw();
                    }
                } );
            } );
        ")
        
        btn <- JS('["copy", "csv", "xls", "pdf", "print",
                {
                    "sExtends": "copy",
                    "sButtonText": "Copy columns",
                    "mColumns": [ 0, 1, 4 ],
                      oSelectorOpts: {
        filter: "applied"
    }
                },
                {
                    "sExtends": "ajax",
                    "sButtonText": "Visible columns",
                    "mColumns": "visible"
                },
{"sExtends": "text", "sButtonText": "Select filtered", "fnClick": function ( node, conf ) {
                 this.fnSelectAll( true )
                 alert("Total selections: "+ this.fnGetSelectedData().length +" rows!")
               } }, "select_none"
            ]')
        
        datatable(
            rv$table, 
            server = TRUE, 
            rownames = FALSE,
            extensions = 'TableTools',
            container = sketch,
            callback = callback,
            options = list(
                processing = TRUE,
                order=JS('[[ 7, "asc" ]]'),
                pageLength = 10,
                lengthMenu=JS('[[10, 25, 50, 100, 1000, -1], [10, 25, 50, 100, 1000, "All"]]'),
                columns = JS(readLines('colDef.js')),
                dom = 'T<"clear">lfrtip',
                tableTools = list(aButtons=btn, sRowSelect="multi", sSwfPath = copySWF(dest='www', pdf = TRUE)),
                ajax = list(
                    url = action, 
                    type = 'POST', 
                    data = JS(
                        'function(d) {',
                        'd.search.caseInsensitive = false;',
                        'd.escape = true;',
                        '}'
                    )
                )
            )
        )  %>% formatRound(4:5, 2)    
        
    })

#     output$data <- renderDataTable(
#         datatable(rv$table)
# #         ,options = list(
# #             pageLength = 10,
# #             order=I('[[ 7, "asc" ]]'),
# #             lengthMenu=I('[[10, 25, 50, 100, -1], [10, 25, 50, 100, "All"]]'),
# #             columns = I(readLines('colDef.js'))
# #         )
#     )

    
    output$downloadData <- downloadHandler(
        filename = function() {
            paste('DEview_data_', gsub(' ', '_', Sys.time()), '.csv', sep='')
        },
        content = function(con) {
            write.csv(rv$table[,-ncol(rv$table)], con, row.names = FALSE)
        }
    )
    
    output$downloadFigure <- downloadHandler(
        filename = function() {
            paste('DEview_figure_', gsub(' ', '_', Sys.time()), '.pdf', sep='')
        },
        content = function(con) {
            if(is.null(rv$plot)) stop('Plot something first (green button in result table).')
            pdf(con, height=6, width=12)
            print(rv$plot)
            dev.off()
            #ggsave(con, rv$plot, height=6, width=12)
        }
    )
    
    output$downloadDataFlt = downloadHandler(
        filename = function() {
            paste('DEview_data_', gsub(' ', '_', Sys.time()), '.csv', sep='')
        },
        content = function(file) {
            s = input$data_rows_all
            message(length(s))
            message(length(input$data_rows_current))
            write.csv(rv$table[s, -ncol(rv$table), drop = FALSE], file, row.names = FALSE)
        }
    )
    

})
