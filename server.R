#   ScanGEO: parallel mining of high-throughput gene expression data
#   Copyright (C) <2017>  <Katja Koeppen>

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
#    General Public License for more details <http://www.gnu.org/licenses/>.

shinyServer(function(input, output, session){
               

        datasetInput <- eventReactive(input$find, {
                query <- paste("select title, gds, pubmed_id, type, platform_organism,
                                update_date from gds where platform_organism like",
                                "'%",input$organism, "%'", "and description like", "'%",
                                toupper(input$text),"%'", sep="")
                dbGetQuery(con, query)
                })
        
        
        output$Dim <- renderText(dim(datasetInput())[1])
        
        
        observeEvent(input$find, {
        
        if (input$text > 0){
                if (dim(datasetInput())[1] == 0) {
                        output$message <- renderText(paste(dim(datasetInput())[1],
                                          "studies found \n searching for \"", isolate(input$text), "\"", 
                                          "\n in", isolate(input$organism),
                                          "\n Try modifying your search term"))        
                } else {
                        output$message <- renderText(paste(dim(datasetInput())[1],
                                          "studies found \n searching for \"", isolate(input$text), "\"",  
                                          "\n in", isolate(input$organism)))
                        }
        } else {
        output$message <- renderText(paste(dim(datasetInput())[1],
                          "studies found", "in", isolate(input$organism)))
        }
        
        output$status <- renderUI({
                        tagList(
                                verbatimTextOutput('message'),
                                hr(),
                                downloadButton('downloadTable', "Download summary table")
                        )
                })
        
        output$table <- renderDataTable({
                datasetInput()
        }, options = list(searching = FALSE))
        
        output$ScanButton <- renderUI({
                actionButton('scan', "ScanGEO", icon("spinner"),
                style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
        })
        })
        
   
        
        output$downloadTable <- downloadHandler(
                filename = function() {
                        paste(input$organism, input$text, "Table", '.csv', sep='')
                },
                content = function(file) {
                        write.csv(datasetInput(), file, row.names = FALSE)
                }
        )
       
        genesList <- reactive({
                Org[[input$organism]]
        })
        
        KEGGList <- reactive({
                names(MasterListAll[[input$organism]])
        })
        
        observe({
        updateSelectizeInput(session, 'gene', server = TRUE, choices = genesList())
        updateSelectizeInput(session, 'KEGG', server = TRUE, choices = KEGGList())
        })
        
        
        custom_genes <- reactive({
                if (is.null(input$genefile))
                        return(NULL)
                read.csv(input$genefile$datapath, header = FALSE, nrows = 200,
                         stringsAsFactors = FALSE)[, 1]
        })
        
        
        output$customfile <- renderUI({
                fileInput('genefile', 'Upload CSV file (limit = 200 gene symbols)', accept='.csv')
        })
        
        observeEvent(input$genefile, {
        unmapped <- reactive({
                setdiff(custom_genes(), unlist(Org[input$organism]))
                })
        output$unmapped <- renderText({unmapped()})
        
        if (length(unmapped()) > 0){
                output$error <- renderUI({
                        tagList(
                        h5("The following gene symbols were not found in the selected organism:"),
                        verbatimTextOutput('unmapped')
                        )
                })
        } else {
                output$error <- renderUI({})      
             }
        })
        
        observeEvent({length(input$gene) != 0 | length(input$KEGG) != 0 | length(input$genefile) != 0}, {
                Genes <- unique(c(as.character(unlist(MasterListAll[[input$organism]][input$KEGG])), as.character(input$gene), custom_genes()))
                if (length(Genes) != 0) {
                output$ScanTime <- renderText(paste("Estimated scan time:", 
                round(predict(RF, newdata = data.frame(Studies = dim(datasetInput())[1], Geneinput = length(Genes)))/60, 1), "minutes"))
                output$time <- renderUI({
                        verbatimTextOutput('ScanTime')
                        })
                }
        })
        
        
        observeEvent(input$scan, {
                
                Genes <- unique(c(as.character(unlist(MasterListAll[[input$organism]][input$KEGG])), as.character(input$gene), custom_genes()))
                
                if (length(Genes) == 0) {
                        output$done <- renderText({'No genes selected'})
                        output$ui <- renderUI({
                                verbatimTextOutput('done')
                        })
                } else {       
                        Session.ID <- gsub(" ", "", Sys.time(), fixed = TRUE)
                        Session.ID <- gsub("-", "", Session.ID, fixed = TRUE)
                        Session.ID <- gsub(":", "", Session.ID, fixed = TRUE)
                        
                        setwd(First.wd)
                        Session.path <- paste(First.wd, "/results/", Session.ID, sep ="")
                        dir.create(Session.path)
                        setwd(Session.path)
                        
                        write.csv(datasetInput(), file = paste("02_GDS_summary_", input$organism, '.csv', sep=''), 
                                  row.names = FALSE)
                        file.copy(from = paste(First.wd, "/01_README.pdf", sep = ""), to = Session.path)
                        
                        withProgress(message = 'Scanning GEO data base', value = 0, {
                                ScanGeo(Genes, datasetInput()$gds, input$alpha)})
                        if (length(list.files(pattern = '\\.pdf$')) > 0) {
                                
                                Composite.Name <- paste(Session.ID, '.zip', sep='')
                                SystemCall <- paste("find . \\( -name '*.pdf' -or -name '*.csv' \\) -print | zip ", 
                                                    Composite.Name, " -@", sep ="")
                                system(SystemCall)
                                
                                output$done <- renderText({'Scan complete!'})
                                output$ui <- renderUI({
                                        tagList(
                                        verbatimTextOutput('done'),
                                        downloadButton('download', 'Download results'),
                                        hr(),
                                        actionButton('reset', "Reset", icon("refresh", class = "fa-spin"),
                                                     style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                                        )
                                })
                        } else {
                                output$done <- renderText({'No significant genes!'})
                                output$ui <- renderUI({
                                tagList(
                                verbatimTextOutput('done'),
                                hr(),
                                actionButton('reset', "Reset", icon("refresh", class = "fa-spin"),
                                             style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                                )
                                })
                        }
               }
           
               

        output$download <- downloadHandler(
                filename = Composite.Name,
                content <- function(file) {
                        file.copy(from = paste(Session.path, "/", Composite.Name, sep = ""), to = file)
                },
                contentType = "application/zip")
        })
        
        
        observeEvent(input$reset, {
                updateSelectInput(session, 'organism', label = "Organism:", choices = Genus)
                updateTextInput(session, 'text', label = "Additional search term:", value = "")
                updateSelectizeInput(session, 'gene', server = TRUE, choices = genesList())
                updateRadioButtons(session, 'alpha', label = "Significance level alpha",
                                   choices = list("0.05" = 0.05, "0.01" = 0.01, "0.001" = 0.001),
                                   selected = 0.05)
                output$table <- renderDataTable({})
                output$ScanButton <- renderUI({})
                output$time <- renderUI({})
                output$ui <- renderUI({})
                output$status <- renderUI({})
                output$error <- renderUI({})
                output$customfile <- renderUI({
                        fileInput('genefile', 'Upload CSV file (limit = 200 gene symbols)', accept='.csv')
                })
        
              
        })

        
}) 
