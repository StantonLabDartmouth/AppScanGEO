shinyUI(
        fluidPage(
                tags$img(style="height:150px; width:150px", src="Scan.png", align = "right"),
                titlePanel("ScanGEO - paralell mining of high-throughput gene expression data"),
                
                
                fluidRow(
                        column(4,
                                wellPanel(
                                selectInput('organism', label = "Organism:",
                                            choices = Genus),

                                textInput('text', label = "Additional search term:"),
                                tags$small("Select an organism, enter one optional
                                           search term and push the button below to
                                           find relevant GEO data sets."),
                                hr(),

                                actionButton('find', "Find GEO data sets", icon("question-circle"), 
                                             style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                                hr(),
                                
                                uiOutput("status")
                                ),


                                wellPanel(

                                selectizeInput('gene', "Enter gene symbols", multiple = TRUE, choices = NULL),
                                hr(),
                                selectizeInput('KEGG', "Select KEGG pathway", multiple = FALSE, choices = NULL),


                                radioButtons('alpha', label = "Significance level alpha",
                                             choices = list("0.05" = 0.05, "0.01" = 0.01, "0.001" = 0.001),
                                             selected = 0.05),

                                tags$small("Enter text to select gene symbols and/or select a KEGG pathway,
                                           chose a significance level and press 'ScanGEO'.
                                           If significant genes are found, pdf plots can be downloaded
                                           with the 'Download results' button once the scan is complete."),
                                hr(),
                                
                                uiOutput("time"),
                                
                                hr(),

                                uiOutput("ScanButton"),
                                
                                hr(),
                                
                                uiOutput("ui")
                                
                                )),
                               
                               
                        column(8,
                                mainPanel(
                                       tabsetPanel(
                                               tabPanel("Table",
                                                        dataTableOutput(outputId="table")),
                                                tabPanel("Documentation",
                                                        tags$iframe(style="height:600px; width:800px; scrolling=yes",
                                                                    src="index.html"),
                                                        h4("Click on slides to start, use arrow buttons to advance or go back")
                                                )

                                       )
                              ))
                        )))
