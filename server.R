library(shiny)
library(DESeq2)
library(ggplot2)
library(scales)
library(DT)

options("shiny.maxRequestSize" = -1)

shinyServer(function(input, output, session) {

    rv <- reactiveValues(info='start', dds=NULL, plot=NULL, table=NULL, SE=NULL)
    if(!file.exists("cache")) if(!dir.create('cache')) stop('FS is not writable.') else message('cache dir created') else message('cache OK')
    
    observe({
        message('loading: ', input$SEdata)
        rv$SE <- get(load(file.path('data', input$SEdata)))
        updateRadioButtons(session, 'type', choices = colnames(colData(rv$SE))[1:2], selected = colnames(colData(rv$SE))[2])
    })
    
    output$distPlot <- renderPlot({
        if(is.null(input$plot)) return()
        gene <- input$plot
        
        
        if(input$plotValues=='a') {
            message('plotsetsetup')
            rv$plotset <- DESeqDataSet(rv$SE, ~ 1)
        } else if (input$plotValues=='f') {
            rv$plotset <- rv$dds
        } else {
            message('plotsetsetup')
            inc <- (colData(rv$SE)[[1]] %in% input$plotValues_f1) & (colData(rv$SE)[[2]] %in% input$plotValues_f2)
            rv$plotset <- DESeqDataSet( rv$SE[,inc], ~1 )
        }
        
        plotset <- rv$plotset
        intgroups <- colnames(colData(rv$SE))
        
        
        data <- plotCounts(plotset, gene, intgroup=intgroups, returnData=TRUE)
        if(input$plotType == 'fpm') data$count <- fpm(plotset)[gene,]
        if(input$plotType == 'fpkm') data$count <- fpkm(plotset)[gene,]
        
        info <- sapply(mcols(plotset[gene,])[,1:2], as.character)
        #bm <- round( mcols(plotset[gene,])[,3], 2) 
        secondGroup <- intgroups[intgroups!='stage']
        
        g <- ggplot(data, aes_string( x="stage", y="count", color=secondGroup, group=secondGroup )) + 
            geom_point() + guides(color=guide_legend(title=secondGroup)) +
            stat_smooth(se=FALSE,method="loess", size = 1.3) +
            ggtitle(paste0(
                info[1],  ' (', gene, ')\n ', info[2] #, '\n FDR=', res[gene,]$padj, '; BM=', bm
            ))
        
        if(input$plotScale == 'N') g <- g + scale_y_continuous(breaks=pretty_breaks(n=10))
        
        if(input$plotScale == 'log10') g <- g + scale_y_continuous(
            trans = 'log10',
            breaks = trans_breaks('log10', function(x) 10^x, n=10),
            #labels = trans_format('log10', math_format(10^.x))
            labels = round
        )

        
        if(input$plotScale == 'log2') g <- g + scale_y_continuous(
            trans = log2_trans(), 
            breaks = trans_breaks('log2', function(x) 2^x, n=10)
            #labels = trans_format('log2', math_format(2^.x))
        )
        
        rv$plot <- g
        g
        
        
    })
    
    sw <- observe({
        #if(is.null(rv$table)) return()
        
        if(input$test == "Wald") {
            cc <- levels(colData(rv$SE)[[input$type]])
            updateRadioButtons(session, 'p1', label=paste0('[', input$type,']'), choices = cc, selected = head(cc, 1) )
            updateRadioButtons(session, 'p2', choices = cc, selected = tail(cc, 1) )
        }
        if(input$test == "Wald") {
            wh <- colnames(colData(rv$SE))[1:2][colnames(colData(rv$SE))[1:2] != input$type]
            updateSelectInput(session, 'which', choices = wh)
        } else {
            wh <- colnames(colData(rv$SE))[1:2]
            updateSelectInput(session, 'which', choices = wh, selected=input$which)
        }
        
        
        ft <- levels(colData(rv$SE)[[input$which]])
        info <- paste0("Use following [",input$which,"]:")
        updateRadioButtons(session, 'what', label=info, choices = ft, selected=head(ft, 1), inline=FALSE )
        
    })

    
    getResultTable <- reactive({
   
            progress <- shiny::Progress$new(session, min=0, max=3)
            on.exit(progress$close())
            
            progress$set(message = 'Initiating...', value = 0)
            
            if(input$test == 'asis') {
                progress$set(message = 'Calculating TS', value = 1)
                load('ddsTC.Rdata')
                res <- results(ddsTC)
                res$log2FoldChange <- res$lfcSE <- NA
                rv$info <- paste0(
                    'Log2 Fold Change and its standard error available for paired tests only!\n', 
                    paste0(elementMetadata(res)[-(2:3),2], collapse='\n')
                )
                rv$dds <- dds <- ddsTC
                
            } else if(input$test == 'LRT') {
                
                
                fltstring <- if(input$filter) paste0(input$which, '_', input$what) else 'AllData'
                
                fn <- file.path('cache', paste0(
                    sub('\\.', '_', input$SEdata),
                    '_cache_LTR_', fltstring, '_', 
                    paste(as.character(as.formula(input$m1)), collapse=''), '_',
                    paste(as.character(as.formula(input$m0)), collapse=''), '.rda'
                ))
                
                if(file.exists(fn)) {
                    progress$set(message = 'Loading from cache', value = 1, detail=""); message('Loading from cache: ', fn)
                    rv$dds <- dds <- get(load(fn))
                } else {
                    progress$set(message = 'Calculating statistical tests (LRT)', detail = 'This may take a while...', value = 1)
                    if(input$filter) {
                        dds <-DESeqDataSet(rv$SE[,colData(rv$SE)[[input$which]] == input$what], design = as.formula(input$m1) )
                    } else {
                        dds <-DESeqDataSet(rv$SE, design = as.formula(input$m1) )
                    }
                    design(dds) <- as.formula(input$m1)
                    rv$dds <- dds <- DESeq(dds, test="LRT", reduced = as.formula(input$m0) )
                    save(dds, file=fn)
                }
                
                res <- results(dds)
                res$log2FoldChange <- res$lfcSE <- NA
                rv$info <- paste0(
                    'Log2 Fold Change and its standard error available for paired tests only!\n', 
                    paste0(elementMetadata(res)[-(2:3),2], collapse='\n')
                )
                
            } else {
                
                fltstring <- if(input$filter) paste0(input$which, '_', input$what) else 'AllData'
                fn <- file.path('cache', paste0(
                    sub('\\.', '_', input$SEdata), 
                    '_cache_Wald_', fltstring, '_', input$type, '.rda')
                )
                
                if(file.exists(fn)) {
                    progress$set(message = 'Loading from cache', value = 1, detail=''); message('Loading from cache: ', fn)
                    rv$dds <- get(load(fn))
                } else {
                    progress$set(message = 'Calculating statistical tests', value = 1, detail = 'This may take a while...')
                    if(input$filter) {
                        dds <- DESeqDataSet(
                            rv$SE[,colData(rv$SE)[[input$which]] == input$what], 
                            design = as.formula(paste('~', input$type))
                        )
                    } else {
                        dds <- DESeqDataSet(rv$SE, design = as.formula(paste('~', input$type))) 
                    }
                    rv$dds <- dds <- DESeq(dds)
                    save(dds, file=fn)
                }
                
                progress$set(message = 'Calculating Wald tests results', value = 2, detail=paste('Contrast: ', input$type, input$p1, input$p2))
                res <- results(dds, contrast=c(input$type, input$p1, input$p2), test="Wald", addMLE=TRUE)
                rv$info <- paste0(elementMetadata(res)[,2], collapse='\n')
            }
            
            progress$set(message = 'Building result table', value = 3, detail = '')
            
            tab <- data.frame(
                ID=rownames(res), 
                gene=if(length(mcols(dds)$geneName)) as.character(mcols(dds)$geneName) else "",
                
                seqID=if(length(mcols(dds)$seqID)) as.character(mcols(dds)$seqID) else "",
                info=if(length(mcols(dds)$desc)) as.character(mcols(dds)$desc)else "",
                as.data.frame(res), 
                Plot=NA,
                stringsAsFactors = FALSE
            )
            return(tab)
   
    })
    
    output$info <- renderText({
        if( input$apply == 0 ) return('Press "Apply settings" to see the results table.')
        rv$table <- isolate( getResultTable() )
        return(rv$info)
    })

    output$data = DT::renderDataTable({
        if(is.null(rv$table)) return()
        
        action = session$registerDataObj('iris', rv$table, shiny:::dataTablesJSON)
        sketch = htmltools::withTags(table(
            tableHeader(rv$table),
            tableFooter(rep('', ncol(rv$table)))
        ))
        callback = JS("
    
        ")
        
        btn <- JS('[{"sExtends": "div", "sButtonClass": "button-label"},
            "copy", "csv", "xls", "pdf", "print",
            {"sExtends": "text", "sButtonText": "Sel. all", "fnClick": function ( node, conf ) {
                 this.fnSelectAll( true );
                 alert("Total selections: "+ this.fnGetSelectedData().length +" rows!");
               } }, "select_none"
            ]')
        
        datatable(
            rv$table, 
            #server = TRUE, 
            rownames = FALSE,
            extensions = c('TableTools', 'ColReorder', 'ColVis'),
            container = sketch,
            callback = JS(readLines('callback.js')),,
            options = list(
                processing = TRUE,
                serverSide = TRUE,
                
                deferRender = TRUE,
                scrollY = 385,
                
                colReorder = list(realtime = TRUE),
                
                order = JS('[[ 9, "asc" ]]'),
                pageLength = 10,
                lengthMenu = JS('[[10, 25, 50, 100, 1000, -1], [10, 25, 50, 100, 1000, "All"]]'),
                columns = JS(readLines('colDef.js')),
                dom = 'RCT<"clear">lfrtip',
                tableTools = list(aButtons=btn, sRowSelect="os", sSwfPath = copySWF(dest='www', pdf = TRUE)),
                #fnServerParams= JS("function(params) {Shiny.shinyapp.sendInput({sel:this.DataTable().ajax.params()});}"),
                ajax = list(
                    url = action, 
                    type = 'POST'
                )
            )
        )  %>% formatRound(5:9, 2)
        
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
            out <- rv$table[,-ncol(rv$table)]
            if('R' %in% input$add) out <- cbind(out, counts(rv$dds, normalized=FALSE))
            if('NR' %in% input$add) out <- cbind(out, counts(rv$dds, normalized=TRUE))
            #if('RPM' %in% input$add) out <- cbind(out, fpm(rv$dds))
            if('RPKM' %in% input$add) out <- cbind(out, fpkm(rv$dds, robust = FALSE))
            write.csv(out, con, row.names = FALSE)
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
        content = function(con) {
#             s = input$selected
#             message(length(s))
#             message(length(input$data_rows_current))
            # counts(rv$dds, normalized=FALSE)[ rownames(getFLT(rv$table, input$sel)), ]
            out <- getFLT(rv$table, input$sel)
            out <- out[, -ncol(rv$table), drop = FALSE]
            if('R' %in% input$add) out <- cbind(out, counts(rv$dds, normalized=FALSE)[rownames(out),,drop = FALSE] )
            if('NR' %in% input$add) out <- cbind(out, counts(rv$dds, normalized=TRUE)[rownames(out),,drop = FALSE] )
            if('RPM' %in% input$add) out <- cbind(out, fpm(rv$dds)[rownames(out),,drop = FALSE] )
            if('RPKM' %in% input$add) out <- cbind(out, fpkm(rv$dds)[rownames(out),,drop = FALSE] )
            write.csv(out, con, row.names = FALSE)
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

output$design <- DT::renderDataTable({
    #tmp <- if(input$desall) SE else SE[,colData(SE)[[input$which]] == input$what]
    des <- as.data.frame(colData(rv$SE))
    rownames(des) <- 1:nrow(des)
    if(!input$desall) des <- des[des[[input$which]] == input$what,]

    datatable(des, rownames = TRUE)
})

output$debug_out <- renderPrint({
    if(input$debug_submit==0) return()
    isolate( eval(parse(text=input$debug_cmd)) )
})
   
observe({
    file.copy(input$newfile$datapath, file.path('data', input$newfile$name))
    updateSelectInput(session, inputId = 'SEdata', label = 'Dataset', choices = dir('data') )
}) 

output$dlfile <- downloadHandler(input$SEdata, content = function(x) { 
    file.copy(file.path('data', input$SEdata), x) 
})

observe({
    if(input$rmfile == 0) return()
    isolate({
        file.remove(file.path('data', input$SEdata))
        file.remove(file.path('cache', dir('cache', pattern=paste0(sub('\\.', '_', input$SEdata), '_cache'))))
        updateSelectInput(session, inputId = 'SEdata', label = 'Dataset', choices = dir('data') )
    })
}) 

})
