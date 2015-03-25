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
                res$log2FoldChange <- res$lfcSE <- NA
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
            
            tab <- data.frame(
                ID=rownames(res), 
                gene=as.character(mcols(ddsTC)$geneName),
                info=as.character(mcols(ddsTC)$desc),
                as.data.frame(res), 
                Plot=NA,
                stringsAsFactors = FALSE
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
    
        ")
        
        btn <- JS('[
            "copy", "csv", "xls", "pdf", "print",
            {"sExtends": "text", "sButtonText": "Select all visible", "fnClick": function ( node, conf ) {
                 this.fnSelectAll( true );
                 alert("Total selections: "+ this.fnGetSelectedData().length +" rows!");
               } }, "select_none"
            ]')
        
        datatable(
            rv$table, 
            server = TRUE, 
            rownames = FALSE,
            extensions = c('TableTools'),
            container = sketch,
            callback = JS(readLines('callback.js')),,
            options = list(
                processing = TRUE,
                serverSide= TRUE,
                
                deferRender = TRUE,
                scrollY = 385,
                
                order=JS('[[ 7, "asc" ]]'),
                pageLength = 10,
                lengthMenu=JS('[[10, 25, 50, 100, 1000, -1], [10, 25, 50, 100, 1000, "All"]]'),
                columns = JS(readLines('colDef.js')),
                dom = 'T<"clear">lfrtip',
                tableTools = list(aButtons=btn, sRowSelect="os", sSwfPath = copySWF(dest='www', pdf = TRUE)),
                #fnServerParams= JS("function(params) {Shiny.shinyapp.sendInput({sel:this.DataTable().ajax.params()});}"),
                ajax = list(
                    url = action, 
                    type = 'POST'
                )
            )
        )  %>% formatRound(4:7, 2)    
        
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

    getFLT <- function(data, q) {
        n <- nrow(data)
        ci <- q$search[["caseInsensitive"]] == TRUE
        i <- seq_len(n)
        if (q$search[["value"]] != "") {
            i0 <- apply(data, 2, function(x) {
                shiny:::grep2(q$search[["value"]], as.character(x), fixed = q$search[["regex"]] == 
                          FALSE, ignore.case = ci)
            })
            i <- intersect(i, unique(unlist(i0)))
        }
        if (length(i)) 
            for (j in seq_len(length(q$columns))) {
                col <- q$columns[[j]]
                if (col[["searchable"]] != TRUE) 
                    next
                if ((k <- col[["search"]][["value"]]) == "") 
                    next
                j <- as.integer(j)
                dj <- data[, j]
                r <- shiny:::commaToRange(k)
                ij <- if (length(r) == 2 && is.numeric(dj)) {
                    which(dj >= r[1] & dj <= r[2])
                } else {
                    shiny:::grep2(k, as.character(dj), fixed = col[["search"]][["regex"]] == 
                              FALSE, ignore.case = ci)
                }
                i <- intersect(ij, i)
                if (length(i) == 0) 
                    break
            }
        if (length(i) != n) 
            data <- data[i, , drop = FALSE]

        
        for (ord in q$order) {
            k <- ord[["column"]]
            d <- ord[["dir"]]
            
            data <- data[order(data[,as.integer(k) + 1], decreasing = (d!="asc")), , drop = FALSE]
            
        }

        return(data)
    }
    
    output$downloadDataFlt = downloadHandler(
        filename = function() {
            paste('DEview_data_', gsub(' ', '_', Sys.time()), '.csv', sep='')
        },
        content = function(file) {
#             s = input$selected
#             message(length(s))
#             message(length(input$data_rows_current))
            out <- getFLT(rv$table, input$sel)
            write.csv(out[, -ncol(rv$table), drop = FALSE], file, row.names = FALSE)
        }
    )
# 
# output$data2 <- DT::renderDataTable({    
#     action = dataTableAjax(session, rv$table, rownames = FALSE)
#     datatable(
#         rv$table, 
#         #server = TRUE, 
#         rownames = FALSE,
# 
#         options = list(
#             #ajax = list(url = action),
#             columns = JS(readLines('colDef.js'))
#         )
#     )
# })

output$debug_out <- renderPrint({
    if(input$debug_submit==0) return()
    isolate( eval(parse(text=input$debug_cmd)) )
})
    

})
