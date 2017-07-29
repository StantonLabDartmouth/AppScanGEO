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

#    If this code is helpful to your research, please cite our related publication:
#    https://academic.oup.com/bioinformatics/article-abstract/doi/10.1093/bioinformatics/btx452/3965322/ScanGEO-parallel-mining-of-highthroughput-gene

shinyUI(
        fluidPage(
                tags$img(style="height:150px; width:150px", src="Scan.png", align = "right"),
                titlePanel("ScanGEO - parallel mining of high-throughput gene expression data"),
                
                fluidRow(
                        column(4,
                               br(),
                               br(),
                               br(),
                               br(),
                               tabsetPanel(
                                tabPanel("Select Studies",        
                                wellPanel(
                                selectInput('organism', label = "Organism:", choices = Genus),

                                textInput('text', label = "Additional search term:"),
                                tags$small("Select an organism, enter one optional
                                           search term and push the button below to
                                           find relevant GEO data sets."),
                                hr(),

                                actionButton('find', "Find GEO data sets", icon("question-circle"), 
                                             style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                                hr(),
                                
                                uiOutput("status")
                                )),
                               
                                tabPanel("KEGG Pathway",
                                wellPanel(
                                selectizeInput('KEGG', "Select KEGG pathway", multiple = FALSE, choices = NULL)
                                )),
                                
                                tabPanel("Custom Genes",
                                wellPanel(
                                selectizeInput('gene', "Enter gene symbols", multiple = TRUE, choices = NULL),
                                hr(),
                                textInput('wildcard', "Wildcard search (e.g. all genes starting with 'MIR' or 'LINC')", value = ""),
                                hr(),
                                uiOutput('wild_UI'),
                                hr(),
                                uiOutput('customfile'),
                                hr(),
                                uiOutput('error')
                                )),
                                
                                tabPanel("Scan",
                                wellPanel(
                                radioButtons('alpha', label = "Significance level alpha", selected = 0.05,
                                             choices = list("0.05" = 0.05, "0.01" = 0.01, "0.001" = 0.001)),

                                tags$small("Chose a significance level and press 'ScanGEO'.
                                           If significant genes are found, pdf plots and csv files can be downloaded
                                           with the 'Download results' button once the scan is complete."),
                                hr(),
                                                        
                                uiOutput("time"),
                                                        
                                hr(),
                                                        
                                uiOutput("ScanButton"),
                                                        
                                hr(),
                                                        
                                uiOutput("ui")
                                ))
                        
                        )),
                        
                        column(8,       
                                mainPanel(
                                        tabsetPanel(
                                        tabPanel("Table",
                                                dataTableOutput(outputId="table")),
                                        tabPanel("Significant Genes",
                                                dataTableOutput(outputId="sig_genes")),
                                        tabPanel("Significant Studies",
                                                dataTableOutput(outputId="sig_studies")),
                                        tabPanel("Documentation",
                                                tags$iframe(style="height:600px; width:800px; scrolling=yes", src="index.html"),
                                                h4("Click on slides to start, use arrow buttons to advance or go back")
                                                )))
                                       
                              
                        )
                )
        )
)
