#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
options(shiny.maxRequestSize=70*1024^2)
options(stringsAsFactors = F)

library(shiny)
library(shinythemes)
library(reshape2)
library(edgeR)
library(DESeq2)
library(limma)
library(pheatmap)
library(ggplot2)
library(DT)
library(yulab.utils)
library(rvcheck)
library(KEGG.db)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(igraph)
# library(openxlsx)
source("./script/DEG_heatmap_volcanol.R")

# Define UI for application 
ui <- fluidPage(
  theme = shinytheme('readable'),
  #shinythemes::themeSelector()
  #theme = shinytheme('sandstone'),
  navbarPage("ceRNA", 
             tabPanel("Step by Step",
                      fluidRow(
                        column(2,
                               wellPanel(
                                 h4(strong("Input your expression file")),
                                 selectInput("file1",label= "Choose an example or your own data", 
                                             choices = c("Example_array" = "Example1",
                                                         "Example_RNAseq" = "Example1_RNAseq",
                                                         "Your own data" = "load_my_own1")),
                                 conditionalPanel("input.file1 == 'Example1'",
                                                  downloadButton('downloadEx1', 'Download example')),
                                 conditionalPanel("input.file1 == 'Example1_RNAseq'",
                                                  downloadButton('downloadEx1_RNAseq', 'Download example')),
                                 conditionalPanel("input.file1 == 'load_my_own1'",
                                                  fileInput('loadfile1', 'Choose xlsx File', 
                                                            accept=c('.xlsx','text/csv', 'text/comma-separated-values,text/plain', '.csv', '.txt')))
                               ),
                               wellPanel(
                                 h4(strong("Input your sample information file")),
                                 selectInput("file2",label= "Choose an example or your own data", 
                                             choices = c("Example_array"="Example2",
                                                         "Example_RNAseq"="Example2_RNAseq",
                                                         "Your own data" = "load_my_own2")),
                                 conditionalPanel("input.file2 == 'Example2'",
                                                  downloadButton('downloadEx2', 'Download example')),
                                 conditionalPanel("input.file2 == 'Example2_RNAseq'",
                                                  downloadButton('downloadEx2_RNAseq', 'Download example')),
                                 conditionalPanel("input.file2 == 'load_my_own2'",
                                                  fileInput('loadfile2', 'Choose xlsx File', 
                                                            accept=c('.xlsx','text/csv', 'text/comma-separated-values,text/plain', '.csv', '.txt')))
                               ),
                               wellPanel(
                                 h4(strong("Input your gene ID information file")),
                                 selectInput("file3",label= "Choose an example or your own data", 
                                             choices = c("Example_array"="Example3", 
                                                         "Example_RNAseq"="Example3_RNAseq",
                                                         "Your own data" = "load_my_own3")),
                                 conditionalPanel("input.file3 == 'Example3'",
                                                  downloadButton('downloadEx3', 'Download example')),
                                 conditionalPanel("input.file3 == 'Example3_RNAseq'",
                                                  downloadButton('downloadEx3_RNAseq', 'Download example')),
                                 conditionalPanel("input.file3 == 'load_my_own3'",
                                                  fileInput('loadfile3', 'Choose xlsx File', 
                                                            accept=c('.xlsx','text/csv', 'text/comma-separated-values,text/plain', '.csv', '.txt')))
                               ),
   
                               conditionalPanel("input.cPanels1 == 2",
                                                wellPanel(
                                                  h4(strong("Heatmap")),
                                                  selectInput("select_heatmap_filetype", 
                                                              label= "The file type:",
                                                              choices=c("array", "RNAseq"),
                                                              selected = "array",
                                                              multiple = FALSE),
                                                  sliderInput("select_heatmap_pvalue", 
                                                              label = "The cutoff of pvalue:",
                                                              min = 0, max = 1, value = 0.05, step = 0.01),
                                                  sliderInput("select_heatmap_fc", 
                                                              label = "The cutoff of fold change:",
                                                              min = 0, max = 5, value = 1, step = 0.1),
                                                  selectInput("select_heatmap_Scale", 
                                                              label= "scale",
                                                              choices=c("row", "column", "none"),
                                                              selected = "row",
                                                              multiple = FALSE),
                                                  checkboxInput("select_heatmap_Clusterrows", 
                                                                label= "Clustered by rows",
                                                                TRUE),
                                                  checkboxInput("select_heatmap_Clustercols", 
                                                                label= "Clustered by columns",
                                                                TRUE),
                                                  checkboxInput("select_Annotation_legend", 
                                                                label= "Show annatation legend",
                                                                TRUE),
                                                  checkboxInput("select_Showrownames", 
                                                                label= "Show row names",
                                                                FALSE)
                                                )),
                               
                               conditionalPanel("input.cPanels1 == 2",
                                                h4(strong("Download")),
                                                wellPanel(
                                                  textInput("fnameHeatmap", "filename", value = "Heatmap"),
                                                  downloadButton('DownloadHeatmap', 'Download Heatmap'),
                                                  br(),
                                                  br(),
                                                  downloadButton('DownloadDEGheatmaptable', 'Download DEG table')
                                                )),
                               
                               
                               conditionalPanel("input.cPanels1 == 3",
                                                wellPanel(
                                                  h4(strong("Volcano plot")),
                                                  selectInput("select_volcano_filetype", 
                                                              label= "The file type:",
                                                              choices=c("array", "RNAseq"),
                                                              selected = "array",
                                                              multiple = FALSE),
                                                  sliderInput("select_volcano_pvalue", 
                                                              label = "The cutoff of pvalue:",
                                                              min = 0, max = 1, value = 0.05, step = 0.01),
                                                  sliderInput("select_volcano_fc", 
                                                              label = "The cutoff of fold change:",
                                                              min = 0, max = 5, value = 1.5, step = 0.1),
                                                  textInput("select_volcano_Title", "Title name", value = "PD vs NC")
                                                )),
                               
                               conditionalPanel("input.cPanels1 == 3",
                                                h4(strong("Download")),
                                                wellPanel(
                                                  textInput("fnamevolcano", "filename", value = "Volcano"),
                                                  downloadButton('Downloadvolcano', 'Download volcano plot'),
                                                  br(),
                                                  br(),
                                                  downloadButton('DownloadDEGvolcanotable', 'Download DEG table')
                                                )),

                               conditionalPanel("input.cPanels1 == 4",
                                                wellPanel(
                                                  h4(strong("Enrichment")),
                                                  selectInput("select_enrichment_filetype", 
                                                              label= "The file type:",
                                                              choices=c("array", "RNAseq"),
                                                              selected = "array",
                                                              multiple = FALSE),
                                                  selectInput("select_enrichment_genesymbol", 
                                                              label= "The up or down regulated differential gene:",
                                                              choices=c("up", "down"),
                                                              selected = "up",
                                                              multiple = FALSE),
                                                  selectInput("select_enrichment_methods", 
                                                              label= "The methods of enrichment:",
                                                              choices=c("KEGG","GO","Reactome"),
                                                              selected = "KEGG",
                                                              multiple = FALSE),
                                                  sliderInput("select_enrichment_pvalue", 
                                                              label = "The pvalue of differential genes:",
                                                              min = 0, max = 1, value = 0.05, step = 0.01),
                                                  sliderInput("select_enrichment_fc", 
                                                              label = "The fold change of differential genes:",
                                                              min = 0, max = 5, value = 1, step = 0.1),
                                                  sliderInput("select_enrichment_pvaluecutoff", 
                                                              label = "The pvalue of enrichment:",
                                                              min = 0, max = 1, value = 1, step = 0.01),
                                                  sliderInput("select_enrichment_qvaluecutoff", 
                                                              label = "The qvalue of enrichment:",
                                                              min = 0, max = 1, value = 1, step = 0.01)
                                                )),
                               
                               conditionalPanel("input.cPanels1 == 4",
                                                h4(strong("Download")),
                                                wellPanel(
                                                  textInput("fnameenrichment", "filename", value = "Enrichment"),
                                                  downloadButton('Downloadenrichment', 'Download enrichment plot'),
                                                  br(),
                                                  br(),
                                                  downloadButton('Downloadenrichmenttable', 'Download enrichment table')
                                                )),
                               
                               conditionalPanel("input.cPanels1 == 5",
                                                wellPanel(
                                                  h4(strong("Classification")),
                                                  selectInput("select_classification_filetype", 
                                                              label= "The file type:",
                                                              choices=c("array", "RNAseq"),
                                                              selected = "array",
                                                              multiple = FALSE),
                                                  selectInput("select_classification_type", 
                                                              label= "Show the type of  differential genes:",
                                                              choices=c("All", "mRNA","miRNA","lncRNA"),
                                                              selected = "mRNA",
                                                              multiple = FALSE),
                                                  sliderInput("select_classification_pvalue", 
                                                              label = "The pvalue of differential genes:",
                                                              min = 0, max = 1, value = 0.05, step = 0.01),
                                                  sliderInput("select_classification_fc", 
                                                              label = "The fold change of differential genes:",
                                                              min = 0, max = 5, value = 1.5, step = 0.1)
                                                )),
                               
                               conditionalPanel("input.cPanels1 == 5",
                                                h4(strong("Download")),
                                                wellPanel(
                                                  textInput("fnameclassification", "filename", value = "classification"),
                                                  downloadButton('Downloadclassificationtable', 'Download classification table')
                                                )),
                               
                               conditionalPanel("input.cPanels1 == 6",
                                                wellPanel(
                                                  h4(strong("DEG Correlation")),
                                                  selectInput("select_correlation_filetype", 
                                                              label= "The file type:",
                                                              choices=c("array", "RNAseq"),
                                                              selected = "array",
                                                              multiple = FALSE),
                                                  selectInput("select_correlation_type", 
                                                              label= "Show the type of correlation:",
                                                              choices=c("miRNA-mRNA","miRNA-lncRNA","all"),
                                                              selected = "miRNA-mRNA",
                                                              multiple = FALSE),
                                                  selectInput("select_calculate_cor_methods", 
                                                              label= "The method of correlation:",
                                                              choices=c("pearson","spearman"),
                                                              selected = "pearson",
                                                              multiple = FALSE),
                                                  sliderInput("select_calculate_cor_DEG_p", 
                                                              label = "The pvalue of differential genes:",
                                                              min = 0, max = 1, value = 0.05, step = 0.01),
                                                  sliderInput("select_calculate_cor_DEG_fc", 
                                                              label = "The fold change of differential genes:",
                                                              min = 0, max = 5, value = 1, step = 0.1),
                                                  sliderInput("select_calculate_cor_pvalue", 
                                                              label = "The pvalue of correlation:",
                                                              min = 0, max = 1, value = 0.05, step = 0.01),
                                                  sliderInput("select_calculate_cor_corvalue", 
                                                              label = "The value of correlation:",
                                                              min = -1, max = 1, value = -0.4, step = 0.01),
                                                  # enrichment
                                                  selectInput("select_cor_enrichment_methods", 
                                                              label= "The methods of enrichment:",
                                                              choices=c("KEGG","GO","Reactome"),
                                                              selected = "KEGG",
                                                              multiple = FALSE),
                                                  sliderInput("select_cor_enrichment_pvaluecutoff", 
                                                              label = "The pvalue of enrichment:",
                                                              min = 0, max = 1, value = 1, step = 0.01),
                                                  sliderInput("select_cor_enrichment_qvaluecutoff", 
                                                              label = "The qvalue of enrichment:",
                                                              min = 0, max = 1, value = 1, step = 0.01)
                                                )),
                               
                               conditionalPanel("input.cPanels1 == 6",
                                                h4(strong("Download")),
                                                wellPanel(
                                                  textInput("fnamecorrelation", "filename", value = "correlation"),
                                                  downloadButton('Downloadcorrelationtable', 'Download correlation table'),
                                                  br(),
                                                  br(),
                                                  downloadButton('Downloadcorenrichment', 'Download enrichment plot'),
                                                  br(),
                                                  br(),
                                                  downloadButton('Downloadecorenrichmenttable', 'Download enrichment table')
                                                )),
                               
                               conditionalPanel("input.cPanels1 == 7",
                                                wellPanel(
                                                  h4(strong("Prediction")),
                                                  selectInput("select_prediction_filetype", 
                                                              label= "The file type:",
                                                              choices=c("array", "RNAseq"),
                                                              selected = "array",
                                                              multiple = FALSE),
                                                  selectInput("select_prediction_type", 
                                                              label= "Show the type of prediction:",
                                                              choices=c("miRNA-lnc","miRNA-mRNA","miRNA-lnc-mRNA"),
                                                              selected = "miRNA-mRNA",
                                                              multiple = FALSE),
                                                  sliderInput("select_calculate_pre_DEG_p", 
                                                              label = "The pvalue of differential genes:",
                                                              min = 0, max = 1, value = 0.05, step = 0.01),
                                                  sliderInput("select_calculate_pre_DEG_fc", 
                                                              label = "The fold change of differential genes:",
                                                              min = 0, max = 5, value = 1, step = 0.1),
                                                  # enrichment
                                                  selectInput("select_prediction_enrichment_methods", 
                                                              label= "The methods of enrichment:",
                                                              choices=c("KEGG","GO","Reactome"),
                                                              selected = "KEGG",
                                                              multiple = FALSE),
                                                  sliderInput("select_prediction_enrichment_pvaluecutoff", 
                                                              label = "The pvalue of enrichment:",
                                                              min = 0, max = 1, value = 1, step = 0.01),
                                                  sliderInput("select_prediction_enrichment_qvaluecutoff", 
                                                              label = "The qvalue of enrichment:",
                                                              min = 0, max = 1, value = 1, step = 0.01),
                                                )),
                               
                               conditionalPanel("input.cPanels1 == 7",
                                                h4(strong("Download")),
                                                wellPanel(
                                                  textInput("fnameprediction", "filename", value = "prediction"),
                                                  downloadButton('Downloadpredictiontable', 'Download prediction table'),
                                                  br(),
                                                  br(),
                                                  downloadButton('Downloadpreenrichment', 'Download enrichment plot'),
                                                  br(),
                                                  br(),
                                                  downloadButton('Downloadepreenrichmenttable', 'Download enrichment table')
                                                )),
                               
                               conditionalPanel("input.cPanels1 == 8",
                                                wellPanel(
                                                  h4(strong("Network")),
                                                  selectInput("select_ce_filetype", 
                                                              label= "The file type:",
                                                              choices=c("array", "RNAseq"),
                                                              selected = "array",
                                                              multiple = FALSE),
                                                  selectInput("select_layout_ce", 
                                                              label= "The algorithm of layout:",
                                                              choices=c("fruchterman.reingold",
                                                                        "circle",
                                                                        "reingold.tilford",
                                                                        "star",
                                                                        "with_drl",
                                                                        "random"),
                                                              selected = "fruchterman.reingold",
                                                              multiple = FALSE),
                                                  # selectInput("select_network_methods", 
                                                  #             label= "The method of correlation:",
                                                  #             choices=c("pearson","spearman"),
                                                  #             selected = "pearson",
                                                  #             multiple = FALSE),
                                                  selectInput("select_vertex.shape_ce", 
                                                              label= "The vertex shape:",
                                                              choices=c("circle","rectangle","none"),
                                                              selected = "circle",
                                                              multiple = FALSE),
                                                  sliderInput("select_network_DEG_p", 
                                                              label = "The pvalue of differential genes:",
                                                              min = 0, max = 1, value = 0.05, step = 0.01),
                                                  sliderInput("select_network_DEG_fc", 
                                                              label = "The fold change of differential genes:",
                                                              min = 0, max = 5, value = 1, step = 0.1),
                                                  sliderInput("select_network_cor_pvalue", 
                                                              label = "The pvalue of correlation:",
                                                              min = 0, max = 1, value = 0.05, step = 0.01),
                                                  sliderInput("select_network_cor_corvalue", 
                                                              label = "The value of correlation:",
                                                              min = -1, max = 1, value = -0.4, step = 0.01),
                                                  sliderInput("select_vertex.size_ce", 
                                                              label = "The value of vertex size:",
                                                              min = 1, max = 10, value = 7, step = 1),
                                                  sliderInput("select_vertex.label.cex_ce", 
                                                              label = "The value of vertex label size:",
                                                              min = 0.1, max = 1.5, value = 0.9, step = 0.1),
                                                  sliderInput("select_vertex.label.dist_ce", 
                                                              label = "The value of vertex label dist:",
                                                              min = 0, max = 2, value = 0, step = 0.1),
                                                  sliderInput("select_edge.width_ce", 
                                                              label = "The value of edge width:",
                                                              min = 0, max = 2, value = 0.2, step = 0.1),
                                                  # enrichment
                                                  selectInput("select_network_enrichment_methods", 
                                                              label= "The methods of enrichment:",
                                                              choices=c("KEGG","GO","Reactome"),
                                                              selected = "KEGG",
                                                              multiple = FALSE),
                                                  sliderInput("select_network_enrichment_pvaluecutoff", 
                                                              label = "The pvalue of enrichment:",
                                                              min = 0, max = 1, value = 1, step = 0.01),
                                                  sliderInput("select_network_enrichment_qvaluecutoff", 
                                                              label = "The qvalue of enrichment:",
                                                              min = 0, max = 1, value = 1, step = 0.01)
                                                )),
                               
                               conditionalPanel("input.cPanels1 == 8",
                                                h4(strong("Download")),
                                                wellPanel(
                                                  textInput("fnamenetwork", "filename", value = "Network"),
                                                  downloadButton('Downloadnetwork', 'Download network plot'),
                                                  br(),
                                                  br(),
                                                  downloadButton('Downloadceenrichment', 'Download enrichment'),
                                                  br(),
                                                  br(),
                                                  downloadButton('Downloadceenrichmentable', 'Download enrichment plot'),
                                                ))
                               
                        ),
                        
                        column(10,
                               tabsetPanel(
                                 tabPanel("Manual", htmlOutput("ReadMe1"), value =1),
                                 tabPanel("Heatmap plot", plotOutput("heatmap1", height= 800, width =800), DT::dataTableOutput("DEGheatmaptable",width = 800),value = 2),
                                 tabPanel("Volcano plot", plotOutput("volcano", height= 800, width = 1000), DT::dataTableOutput("DEGvolcanotable",width = 800),value = 3),
                                 tabPanel("Enrichment", plotOutput("enrichment", height= 600, width = 1000),DT::dataTableOutput("enrichmenttable",width = 800),value = 4),
                                 tabPanel("RNA classification", DT::dataTableOutput("classificationtable",width = 800),value = 5),
                                 tabPanel("ceRNA network1", DT::dataTableOutput("correlationtable",width = 1000),plotOutput("correlationEnrichment",height= 500,width = 1000),
                                          #DT::dataTableOutput("corenrichmenttable",height= 400,width =600),
                                          value = 6),                                 
                                 tabPanel("ceRNA network2", DT::dataTableOutput("predictiontable",width = 1000),plotOutput("predictionEnrichment", height= 500, width = 1000),value = 7),
                                 tabPanel("ceRNA network3",  textOutput("networktext"),plotOutput("network", height=800, width = 800),
                                          plotOutput("networkEnrichment", height= 500, width =1000),value = 8),
                                 id = "cPanels1"
                               )                
                               
                        ),
                        
                        column(12,
                               tags$head(tags$style(type="text/css", "
                                                    #loadmessage {
                                                    position: fixed;
                                                    bottom: 0px;
                                                    right: 0px;
                                                    width: 100%;
                                                    padding: 5px 0px 5px 0px;
                                                    text-align: center; 
                                                    font-weight: bold;
                                                    font-size: 100%;
                                                    color: #000000;
                                                    background-color: #b8b8b8;
                                                    z-index: 105;
                                                    }
                                                    ")),
                               conditionalPanel(condition="$('html').hasClass('shiny-busy')",
                                                tags$div("Loading...",id="loadmessage"))
                        )
                      )
             ),
             # another Panel
  )
)

# Define server
server <- function(input, output, session) {
  data_input1 <- reactive({
    if(input$file1 == 'Example1'){
      d1 <- read.csv("./example/GSE7621.csv",row.names = 1)
    }else if(input$file1 == 'Example1_RNAseq'){
      d1 <- read.csv("./example/GSE136666.csv",row.names = 1)
    }else if(input$file1 == 'load_my_own1'){
      inFile <- input$loadfile1
      if (is.null(inFile))
        return(NULL)
      else if(grepl(".xlsx", inFile[1])) { 
        d1 <-  read.xlsx(as.character(inFile$datapath), colNames = TRUE, rowNames = F)
      }
      else if(grepl(".csv", inFile[1])) { 
        d1 <-  read.csv(as.character(inFile$datapath), header = TRUE, sep = ",", stringsAsFactors = F, as.is = T, fill = T) 
      }
      else if(grepl(".txt", inFile[1])) { 
        d1 <-  read.table(as.character(inFile$datapath), header = TRUE, sep = "\t", stringsAsFactors = F, as.is = T, fill = T)
      }
      #openxlsx::write.xlsx(d2,"tmp.xlsx")
    }
    else 
      return(NULL)
    Dataset1 <- data.frame(d1)
    return(as.data.frame(Dataset1))
  })
  
  output$downloadEx1 <- downloadHandler( 
    filename <- function() {
      paste0('GSE7621','.csv')
    },
    content <- function(file) {
      ds1 <- data_input1()
      write.csv(ds1, file)
    }
  )
  
  output$downloadEx1_RNAseq <- downloadHandler( 
    filename <- function() {
      paste0('GSE136666','.csv')
    },
    content <- function(file) {
      ds1 <- data_input1()
      write.csv(ds1, file)
    }
  )
  
  data_input2 <- reactive({
    if(input$file2 == 'Example2'){
      d2 <- read.csv("./example/samplelist.csv")
    }else if(input$file2 == 'Example2_RNAseq'){
      d2 <- read.csv("./example/samplelist_RNAseq.csv")
    }else if(input$file2 == 'load_my_own2'){
      inFile <- input$loadfile2
      if (is.null(inFile))
        return(NULL)
      else if(grepl(".xlsx", inFile[1])) { 
        d2 <-  read.xlsx(as.character(inFile$datapath), colNames = TRUE, rowNames = F) 
      }
      else if(grepl(".csv", inFile[1])) { 
        d2 <-  read.csv(as.character(inFile$datapath), header = TRUE, sep = ",", stringsAsFactors = F, as.is = T, fill = T) 
      }
      else if(grepl(".txt", inFile[1])) { 
        d2 <-  read.table(as.character(inFile$datapath), header = TRUE, sep = "\t", stringsAsFactors = F, as.is = T, fill = T)
      }
    }
    else 
      return(NULL)
    Dataset2 <- data.frame(d2)
    return(as.data.frame(Dataset2))
  })
  
  output$downloadEx2 <- downloadHandler( 
    filename <- function() {
      paste0("samplelist",".csv")
    },
    content <- function(file) {
      ds2 <- data_input2()
      write.csv(ds2, file, row.names = F)
    }
  )
  
  output$downloadEx2_RNAseq <- downloadHandler( 
    filename <- function() {
      paste0("samplelist",".csv")
    },
    content <- function(file) {
      ds2 <- data_input2()
      write.csv(ds2, file, row.names = F)
    }
  )
  
  data_input3 <- reactive({
    if(input$file3 == 'Example3'){
      d3 <- read.csv("./data/GPL570.csv")
    }else if(input$file3 == 'Example3_RNAseq'){
      d3 <- read.csv("./data/Novaseq6000.csv")
    }else if(input$file3 == 'load_my_own3'){
      inFile <- input$loadfile3
      if (is.null(inFile))
        return(NULL)
      else if(grepl(".xlsx", inFile[1])) { 
        d3 <-  read.xlsx(as.character(inFile$datapath), colNames = TRUE, rowNames = F)
      }
      else if(grepl(".csv", inFile[1])) { 
        d3 <-  read.csv(as.character(inFile$datapath), header = TRUE, sep = ",", stringsAsFactors = F, as.is = T, fill = T) 
      }
      else if(grepl(".txt", inFile[1])) { 
        d3 = read.table(as.character(inFile$datapath), header = TRUE, sep = "\t", stringsAsFactors = F, as.is = T, fill = T) 
      }
    }
    else 
      return(NULL)
    Dataset3 <- data.frame(d3)
    return(as.data.frame(Dataset3))
  })
  
  output$downloadEx3 <- downloadHandler( 
    filename <- function() {
      paste0("GPL570",".csv")
    },
    content <- function(file) {
      ds3 <- data_input3()
      write.csv(ds3, file, row.names = F)
    }
  )
  
  output$downloadEx3_RNAseq <- downloadHandler( 
    filename <- function() {
      paste0("Novaseq6000",".csv")
    },
    content <- function(file) {
      ds3 <- data_input3()
      write.csv(ds3, file, row.names = F)
    }
  )
  
  RNAmap_input  <- reactive({read.csv("./data/RNAmap.csv")})
  miRNA_trans_input <- reactive({read.csv("./data/miRNA_database.csv",row.names = 1,header = T)})
  mi_lnc_input <- reactive({read.csv("./data/database-miRNA-LncRNA.csv",header = T)})
  mi_m_input <- reactive({read.csv("./data/database_mRNA_miRNAN.csv",header = T)})
  
  # one step panel
  data_input12 <- reactive({
    if(input$file12 == 'Example1'){
      d1 <- read.csv("./example/GSE7621.csv",row.names = 1)
    }
    else if(input$file12 == 'load_my_own1'){
      inFile <- input$loadfile12
      if (is.null(inFile))
        return(NULL)
      else if(grepl(".xlsx", inFile[1])) { 
        d1 <-  read.xlsx(as.character(inFile$datapath), colNames = TRUE, rowNames = F)
      }
      else if(grepl(".csv", inFile[1])) { 
        d1 <-  read.csv(as.character(inFile$datapath), header = TRUE, sep = ",", stringsAsFactors = F, as.is = T, fill = T) 
      }
      else if(grepl(".txt", inFile[1])) { 
        d1 <-  read.table(as.character(inFile$datapath), header = TRUE, sep = "\t", stringsAsFactors = F, as.is = T, fill = T)
      }
    }
    else 
      return(NULL)
    Dataset1 <- data.frame(d1)
    return(as.data.frame(Dataset1))
  })
  
  output$downloadEx12<- downloadHandler( 
    filename <- function() {
      paste0('GSE7621','.csv')
    },
    content <- function(file) {
      ds1 <- data_input12()
      write.csv(ds1, file)
    }
  )
  
  data_input22 <- reactive({
    if(input$file22 == 'Example2'){
      d2 <- read.csv("./example/samplelist.csv")
    }
    else if(input$file2 == 'load_my_own2'){
      inFile <- input$loadfile22
      if (is.null(inFile))
        return(NULL)
      else if(grepl(".xlsx", inFile[1])) { 
        d2 <-  read.xlsx(as.character(inFile$datapath), colNames = TRUE, rowNames = F) 
      }
      else if(grepl(".csv", inFile[1])) { 
        d2 <-  read.csv(as.character(inFile$datapath), header = TRUE, sep = ",", stringsAsFactors = F, as.is = T, fill = T) 
      }
      else if(grepl(".txt", inFile[1])) { 
        d2 <-  read.table(as.character(inFile$datapath), header = TRUE, sep = "\t", stringsAsFactors = F, as.is = T, fill = T)
      }
    }
    else 
      return(NULL)
    Dataset2 <- data.frame(d2)
    return(as.data.frame(Dataset2))
  })
  
  output$downloadEx22 <- downloadHandler( 
    filename <- function() {
      paste0("samplelist",".csv")
    },
    content <- function(file) {
      ds2 <- data_input22()
      write.csv(ds2, file)
    }
  )
  
  data_input32 <- reactive({
    if(input$file32 == 'Example3'){
      d3 <- read.csv("./data/GPL570.csv")
    }
    else if(input$file32 == 'load_my_own3'){
      inFile <- input$loadfile32
      if (is.null(inFile))
        return(NULL)
      else if(grepl(".xlsx", inFile[1])) { 
        d3 <-  read.xlsx(as.character(inFile$datapath), colNames = TRUE, rowNames = F)
      }
      else if(grepl(".csv", inFile[1])) { 
        d3 <-  read.csv(as.character(inFile$datapath), header = TRUE, sep = ",", stringsAsFactors = F, as.is = T, fill = T) 
      }
      else if(grepl(".txt", inFile[1])) { 
        d3 = read.table(as.character(inFile$datapath), header = TRUE, sep = "\t", stringsAsFactors = F, as.is = T, fill = T) 
      }
    }
    else 
      return(NULL)
    Dataset3 <- data.frame(d3)
    return(as.data.frame(Dataset3))
  })
  
  output$downloadEx32 <- downloadHandler( 
    filename <- function() {
      paste0("GPL570",".csv")
    },
    content <- function(file) {
      ds3 <- data_input32()
      write.csv(ds3, file)
    }
  )
  
  observe({
    dsnames1 <- colnames(data_input1())
  })
  
  #pheatmap
  heatmapforuse <- function(){
    data1 <- data_input1()
    data2 <- data_input2()
    data3 <- data_input3()
    DEGG <- DEGfunc(data1,data2,data3,
                    type=input$select_heatmap_filetype,
                    pvalue=input$select_heatmap_pvalue,
                    fc=input$select_heatmap_fc)
    DEGpheatmap(data1,
                data2,
                DEGG[[2]],
                Scale=input$select_heatmap_Scale,
                Clusterrows=input$select_heatmap_Clusterrows,
                Clustercols=input$select_heatmap_Clustercols,
                Annotation_legend=input$select_Annotation_legend,
                Showrownames=input$select_Showrownames
    )
   # plot(x=DEGG[[1]][1:100,1],y=DEGG[[1]][1:100,2])
  }
  
  output$heatmap1 <- renderPlot({
    print(heatmapforuse())
  })

  output$DownloadHeatmap <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fnameHeatmap)
      paste(pdf_file,".pdf", sep='')
    },
    content <- function(file) {
      pdf(file, height= 10, width=12,onefile = FALSE)
      print(heatmapforuse())
      dev.off()
      #file.copy(paste(pdf_file,'.pdf', sep='') ,file, overwrite=TRUE)
    }, contentType = 'image/pdf')
  
  output$DEGheatmaptable <- renderDataTable({
    data1 <- data_input1()
    data2 <- data_input2()
    data3 <- data_input3()
    DEGG <- DEGfunc(data1,data2,data3,
                    type=input$select_heatmap_filetype,
                    pvalue=input$select_heatmap_pvalue,
                    fc=input$select_heatmap_fc)
    DEGG[[2]]
  })

  output$DownloadDEGheatmaptable <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fnameHeatmap)
      paste0(pdf_file,'.csv')
    },
    content <- function(file) {
      data1 <- data_input1()
      data2 <- data_input2()
      data3 <- data_input3()
      DEGG <- DEGfunc(data1,data2,data3,
                      type=input$select_heatmap_filetype,
                      pvalue=input$select_heatmap_pvalue,
                      fc=input$select_heatmap_fc)
      write.csv(DEGG[[2]],file)
    })
    
  ##volcano
  volcanoforuse <- function(){
    data1 <- data_input1()
    data2 <- data_input2()
    data3 <- data_input3()
    DEGG <- DEGfunc(data=data1,group_list=data2,gpl=data3,
                    type=input$select_volcano_filetype,
                    pvalue=input$select_volcano_pvalue,
                    fc=input$select_volcano_fc)
    DEGvolcano(DEGG[[1]],
               pvalue=input$select_volcano_pvalue,
               fc=input$select_volcano_fc,
               Title=input$select_volcano_Title)
  }
  
  output$volcano <- renderPlot({
    volcanoforuse()
  })
  
  output$Downloadvolcano <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fnamevolcano)
      paste(pdf_file,".pdf", sep='')
    },
    content <- function(file) {
      pdf(file, height= 10, width=12,onefile = FALSE)
      print(volcanoforuse())
      dev.off()
      #file.copy(paste(pdf_file,'.pdf', sep='') ,file, overwrite=TRUE)
    }, contentType = 'image/pdf')
  
  output$DEGvolcanotable <- renderDataTable({
    data1 <- data_input1()
    data2 <- data_input2()
    data3 <- data_input3()
    DEGG <- DEGfunc(data1,data2,data3,
                    type=input$select_volcano_filetype,
                    pvalue=input$select_volcano_pvalue,
                    fc=input$select_volcano_fc)
    DEGG[[2]]
  })
  
  output$DownloadDEGvolcanotable <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fnamevolcano)
      paste0(pdf_file,'.csv')
    },
    content <- function(file) {
      data1 <- data_input1()
      data2 <- data_input2()
      data3 <- data_input3()
      DEGG <- DEGfunc(data1,data2,data3,
                      type=input$select_volcano_filetype,
                      pvalue=input$select_volcano_pvalue,
                      fc=input$select_volcano_fc)
      write.csv(DEGG[[2]],file)
    })
  # enrichment
  enrichmentforuse <- function(){
    data1 <- data_input1()
    data2 <- data_input2()
    data3 <- data_input3()
    DEGG <- DEGfunc(data1,data2,data3,
                    pvalue=input$select_enrichment_pvalue,
                    type=input$select_enrichment_filetype,
                    fc=input$select_enrichment_fc)
    DE <- DEGG[[2]]
    if (input$select_enrichment_genesymbol=="up"){
      DEu <- DE[DE$logFC >0,]
      enrichment_GO_KEGG_Reactome(DEu$GeneSymbol,
                                  pvaluecutoff=input$select_enrichment_pvaluecutoff,
                                  qvaluecutoff=input$select_enrichment_qvaluecutoff,
                                  methods=input$select_enrichment_methods)
    }else{
      DEd <- DE[DE$logFC <0,]
      enrichment_GO_KEGG_Reactome(DEd$GeneSymbol,
                                  pvaluecutoff=input$select_enrichment_pvaluecutoff,
                                  qvaluecutoff=input$select_enrichment_pvaluecutoff,
                                  methods=input$select_enrichment_methods)
    }
  }
  
  output$enrichment <- renderPlot({
    enrichmentforuse()
  })
  
  output$Downloadenrichment <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fnameenrichment)
      paste(pdf_file,".pdf", sep='')
    },
    content <- function(file) {
      pdf(file, height= 10, width=12,onefile = FALSE)
      enrichmentforuse()
      dev.off()
      #file.copy(paste(pdf_file,'.pdf', sep='') ,file, overwrite=TRUE)
    }, contentType = 'image/pdf')
  
  output$enrichmenttable <- renderDataTable({
    data.frame(enrichmentforuse())
  })
  
  output$Downloadenrichmenttable <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fnameenrichment)
      paste0(pdf_file,'.csv')
    },
    content <- function(file) {
      enrichment_table <- data.frame(enrichmentforuse())
      write.csv(enrichment_table,file,row.names = F)
    })
  
  # Classification
  classficationforuse <- function(){
    data1 <- data_input1()
    data2 <- data_input2()
    data3 <- data_input3()
    DEGG <- DEGfunc(data1,data2,data3,
                    type=input$select_classification_filetype,
                    pvalue=input$select_classification_pvalue,
                    fc=input$select_classification_fc)
    DE <- DEGG[[2]]
    RNAmap  <- RNAmap_input()
    data_all <- merge(DE,RNAmap,by.x="GeneSymbol",by.y="gene_name")
    if (input$select_classification_type=="All"){
      return(data_all)
    }else if(input$select_classification_type=="mRNA"){
      data_mRNA <- data_all[data_all$type=="mRNA",]
      return(data_mRNA)
    }else if(input$select_classification_type=="miRNA"){
      data_miRNA <- data_all[data_all$type=="miRNA",]
      return(data_miRNA)
    }else if(input$select_classification_type=="lncRNA"){
      data_lnc <-data_all[data_all$type=="lncRNA",]
      return(data_lnc)
    }
  }
  
  output$classificationtable <- renderDataTable({
    out <- classficationforuse()
    datatable(out,
              options = list(
                "pageLength" =20))
  })
  
  output$Downloadclassificationtable <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fnameclassification)
      paste0(pdf_file,'.csv')
    },
    content <- function(file) {
      fnameclassification_table <- classficationforuse()
      write.csv(fnameclassification_table,file,row.names = F)
    })
  
  # correlation
  correlationforusepre <- function(){
    data1 <- data_input1()
    data2 <- data_input2()
    data3 <- data_input3()
    RNAmap  <- RNAmap_input()
    DEGG <- DEGfunc(data1,data2,data3,
                    type=input$select_correlation_filetype,
                    pvalue=input$select_calculate_cor_DEG_p,
                    fc=input$select_calculate_cor_DEG_fc)
    DE <- DEGG[[2]]
    data_all <- merge(DE,RNAmap,by.x="GeneSymbol",by.y="gene_name")
    data_lnc <- data_all[data_all$type=="lncRNA",]
    data_mRNA <- data_all[data_all$type=="mRNA",]
    data_miRNA <- data_all[data_all$type=="miRNA",]
    DEG_cor <- DEG_cal_cor(data1,data3,data_lnc,data_mRNA,data_miRNA,
                           method=input$select_calculate_cor_methods,
                           corvalue=input$select_calculate_cor_corvalue,
                           pvalue=input$select_calculate_cor_pvalue)
  }
  data_input_cor  <- reactive(correlationforusepre())
  
  data_input_cor_enrichment_foruse <- function(){
    outcor <- data_input_cor()
    DEG_miRNA_all<- outcor[[3]]
    cha <- unique(c(as.character(DEG_miRNA_all[,1]),as.character(DEG_miRNA_all[,2])))
    enrichment_GO_KEGG_Reactome(cha,pvaluecutoff=input$select_cor_enrichment_pvaluecutoff,
                                qvaluecutoff=input$select_cor_enrichment_qvaluecutoff,
                                methods=input$select_cor_enrichment_methods)
  }
  data_input_cor_enrichment <- reactive(data_input_cor_enrichment_foruse())
  
  output$correlationtable <- renderDataTable({
    outcor <- data_input_cor()
    if (input$select_correlation_type=="all"){
      return(outcor[[3]])
    }else if(input$select_correlation_type=="miRNA-mRNA"){
      return(outcor[[2]])
    }else if(input$select_correlation_type=="miRNA-lncRNA"){
      return(outcor[[1]])
    }
    # datatable(outcor,
    #           options = list(
    #             "pageLength" =20))
  })
  
  output$correlationEnrichment <- renderPlot({
    data_input_cor_enrichment()
  })
  
  # output$corenrichmenttable <- renderDataTable({
  #   data.frame(data_input_cor_enrichment())
  # })
  output$Downloadcorrelationtable <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fnamecorrelation)
      paste0(pdf_file,'.csv')
    },
    content <- function(file) {
      fnamecorrelation_table <- data_input_cor()
      if (input$select_correlation_type=="all"){
        fnamecorrelation_table1 <-  fnamecorrelation_table[[3]]
      }else if(input$select_correlation_type=="miRNA-mRNA"){
        fnamecorrelation_table1 <-  fnamecorrelation_table[[2]]
      }else if(input$select_correlation_type=="miRNA-lncRNA"){
        fnamecorrelation_table1 <-  fnamecorrelation_table[[1]]
      }
      write.csv(fnamecorrelation_table1,file,row.names = F)
    })
  
  output$Downloadcorenrichment <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fnamecorrelation)
      paste(pdf_file,"_enrichment.pdf", sep='')
    },
    content <- function(file) {
      pdf(file, height= 10, width=12,onefile = FALSE)
      data_input_cor_enrichment_foruse()
      dev.off()
      #file.copy(paste(pdf_file,'.pdf', sep='') ,file, overwrite=TRUE)
    }, contentType = 'image/pdf')
  
  output$Downloadecorenrichmenttable <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fnamecorrelation)
      paste0(pdf_file,'_enrichment.csv')
    },
    content <- function(file) {
      enrichment_table <- data.frame(data_input_cor_enrichment())
      write.csv(enrichment_table,file,row.names = F)
    })

  # prediction
  predictionforuse <- function(){
    data1 <- data_input1()
    data2 <- data_input2()
    data3 <- data_input3()
    RNAmap  <- RNAmap_input()
    miRNA_trans <- miRNA_trans_input()
    mi_lnc <- mi_lnc_input()
    mi_m <- mi_m_input()
    mi_mRNA_lnc <- pre_database(data1,data2,data3,RNAmap,miRNA_trans,mi_lnc,mi_m,
                                calculate_DEG_p=input$select_calculate_pre_DEG_p,
                                calculate_DEG_fc=input$select_calculate_pre_DEG_fc,
                                pre_type=input$select_prediction_filetype
                                )
  }
  data_input_pre  <- reactive(predictionforuse())
  
  data_input_pre_enrichment_foruse <- function(){
    outcor <- data_input_pre()
    miRNA_lnc_mRNA<- outcor[["miRNA-lnc-mRNA"]]
    cha <- unique(c(as.character(miRNA_lnc_mRNA[,2]),as.character(miRNA_lnc_mRNA[,3])))
    enrichment_GO_KEGG_Reactome(cha,
                                pvaluecutoff=input$select_prediction_enrichment_pvaluecutoff,
                                qvaluecutoff=input$select_prediction_enrichment_qvaluecutoff,
                                methods=input$select_prediction_enrichment_methods)
  }
  data_input_pre_enrichment <- reactive(data_input_pre_enrichment_foruse())
  
  output$predictiontable <- renderDataTable({
    outcor <- data_input_pre()
    if (input$select_prediction_type=="miRNA-lnc"){
      return(outcor[["miRNA-lnc"]])
    }else if(input$select_prediction_type=="miRNA-mRNA"){
      return(outcor[["miRNA-mRNA"]])
    }else if(input$select_prediction_type=="miRNA-lnc-mRNA"){
      return(outcor[["miRNA-lnc-mRNA"]])
    }
  })
  
  output$predictionEnrichment <- renderPlot({
    data_input_pre_enrichment()
  })
  
  output$Downloadpredictiontable <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fnameprediction)
      paste0(pdf_file,'.csv')
    },
    content <- function(file) {
      fnameprediction_table1 <- data_input_pre()
      if (input$select_prediction_type=="miRNA-lnc"){
        fnameprediction_table <- fnameprediction_table1[["miRNA-lnc"]]
      }else if(input$select_prediction_type=="miRNA-mRNA"){
        fnameprediction_table <- fnameprediction_table1[["miRNA-mRNA"]]
      }else if(input$select_prediction_type=="miRNA-lnc-mRNA"){
        fnameprediction_table <- fnameprediction_table1[["miRNA-lnc-mRNA"]]
      }
      write.csv(fnameprediction_table,file,row.names = F)
    })
  
  output$Downloadpreenrichment<- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fnameprediction)
      paste(pdf_file,"_enrichment.pdf", sep='')
    },
    content <- function(file) {
      pdf(file, height= 10, width=12,onefile = FALSE)
      data_input_pre_enrichment_foruse()
      dev.off()
      #file.copy(paste(pdf_file,'.pdf', sep='') ,file, overwrite=TRUE)
    }, contentType = 'image/pdf')
  
  output$Downloadepreenrichmenttable <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fnameprediction)
      paste0(pdf_file,'_enrichment.csv')
    },
    content <- function(file) {
      enrichment_table <- data.frame(data_input_pre_enrichment())
      write.csv(enrichment_table,file,row.names = F)
    })

  # network
  networkforuse <- function(){
    data1 <- data_input1()
    data2 <- data_input2()
    data3 <- data_input3()
    RNAmap  <- RNAmap_input()
    miRNA_trans <- miRNA_trans_input()
    mi_lnc <-  mi_lnc_input()
    mi_m <- mi_m_input()
    cor_re <- calculate_cor_all(data=data1,
                                group_list=data2,
                                gpl=data3,
                                RNAmap,
                                miRNA_trans,
                                mi_lnc,
                                mi_m,
                                calculate_cor_DEG_p=input$select_network_DEG_p,
                                calculate_cor_DEG_fc=input$select_network_DEG_fc,
                                calculate_cor_pvalue=input$select_network_cor_pvalue,
                                calculate_cor_corvalue=input$select_network_cor_corvalue,
                                type=input$select_ce_filetype
                                )
    cor.mi_lnc <-cor_re[cor_re$genelist2 %in%mi_lnc$lnc_name,,drop=F]
    cor.mi_m <- cor_re[!cor_re$genelist2 %in%mi_lnc$lnc_name,,drop=F]
    #cor.mi_lnc <- cor.mi_lnc[cor.mi_lnc$genelist1 %in%DEG_miRNA_lnc$miRNA & cor.mi_lnc$genelist2 %in%DEG_miRNA_lnc$lncRNA, ]
    #cor.mi_m <- cor.mi_m[cor.mi_m$genelist1%in%DEG_miRNA_mRNA$miRNA & cor.mi_m$genelist2%in% DEG_miRNA_mRNA$mRNA,]
    if(nrow(cor.mi_lnc)==0|nrow(cor.mi_m)==0){
      mi_m_lnc <- data.frame()
    }else{
      mi_m_lnc <- merge(cor.mi_lnc[,1:2],cor.mi_m[,1:2],by="genelist1")
    }
    return(mi_m_lnc)
  }
  networkdata <- reactive(networkforuse())

  networkenrichmentforuse <- function(){
    mi_m_lnc <- networkdata()
    if(nrow(mi_m_lnc)==0){
      print("There is no result, please change your CUTOFF")
    }else{
      cha <- unique(c(as.character(mi_m_lnc[,1]),as.character(mi_m_lnc[,2]),as.character(mi_m_lnc[,3])))
      enrichment_GO_KEGG_Reactome(cha,pvaluecutoff=input$select_network_enrichment_pvaluecutoff,
                                  qvaluecutoff=input$select_network_enrichment_qvaluecutoff,
                                  methods=input$select_network_enrichment_methods)
    }
  }
  
  networkplot <- function(){
    mi_m_lnc <- networkdata()
    if(nrow(mi_m_lnc)==0){
      print("There is no result, please change your CUTOFF")
    }else{
      if (input$select_layout_ce=="fruchterman.reingold"){
        draw_ceRNA(mi_m_lnc, 
                   layout_ce=layout.fruchterman.reingold,
                   vertex.size_ce=input$select_vertex.size_ce,
                   vertex.shape_ce=input$select_vertex.shape_ce,
                   vertex.label.cex_ce=input$select_vertex.label.cex_ce,
                   vertex.label.dist_ce=input$select_vertex.label.dist_ce,
                   edge.width_ce=input$select_edge.width_ce,
                   edge.color_ce="gray",
                   vertex.label.color_ce="black")
      }else if(input$select_layout_ce=="circle"){
        draw_ceRNA(mi_m_lnc, 
                   layout_ce=layout.circle,
                   vertex.size_ce=input$select_vertex.size_ce,
                   vertex.shape_ce=input$select_vertex.shape_ce,
                   vertex.label.cex_ce=input$select_vertex.label.cex_ce,
                   vertex.label.dist_ce=input$select_vertex.label.dist_ce,
                   edge.width_ce=input$select_edge.width_ce,
                   edge.color_ce="gray",
                   vertex.label.color_ce="black")
      }else if(input$select_layout_ce=="reingold.tilford"){
        draw_ceRNA(mi_m_lnc, 
                   layout_ce=layout.reingold.tilford,
                   vertex.size_ce=input$select_vertex.size_ce,
                   vertex.shape_ce=input$select_vertex.shape_ce,
                   vertex.label.cex_ce=input$select_vertex.label.cex_ce,
                   vertex.label.dist_ce=input$select_vertex.label.dist_ce,
                   edge.width_ce=input$select_edge.width_ce,
                   edge.color_ce="gray",
                   vertex.label.color_ce="black")
      }else if(input$select_layout_ce=="star"){
        draw_ceRNA(mi_m_lnc, 
                   layout_ce=layout_as_star,
                   vertex.size_ce=input$select_vertex.size_ce,
                   vertex.shape_ce=input$select_vertex.shape_ce,
                   vertex.label.cex_ce=input$select_vertex.label.cex_ce,
                   vertex.label.dist_ce=input$select_vertex.label.dist_ce,
                   edge.width_ce=input$select_edge.width_ce,
                   edge.color_ce="gray",
                   vertex.label.color_ce="black")
      }else if(input$select_layout_ce=="with_drl"){
        draw_ceRNA(mi_m_lnc, 
                   layout_ce=layout_with_drl,
                   vertex.size_ce=input$select_vertex.size_ce,
                   vertex.shape_ce=input$select_vertex.shape_ce,
                   vertex.label.cex_ce=input$select_vertex.label.cex_ce,
                   vertex.label.dist_ce=input$select_vertex.label.dist_ce,
                   edge.width_ce=input$select_edge.width_ce,
                   edge.color_ce="gray",
                   vertex.label.color_ce="black")
      }else if(input$select_layout_ce=="random"){
        draw_ceRNA(mi_m_lnc, 
                   layout_ce=layout.random,
                   vertex.size_ce=input$select_vertex.size_ce,
                   vertex.shape_ce=input$select_vertex.shape_ce,
                   vertex.label.cex_ce=input$select_vertex.label.cex_ce,
                   vertex.label.dist_ce=input$select_vertex.label.dist_ce,
                   edge.width_ce=input$select_edge.width_ce,
                   edge.color_ce="gray",
                   vertex.label.color_ce="black")
      }
    }
  }
  output$networktext <- renderText({
    mi_m_lnc <- networkdata()
    if(nrow(mi_m_lnc)==0){
     info <- "There is no result, please change your CUTOFF"
     info
    }
  })
  
  output$network <- renderPlot({
    networkplot()
  })
  output$networkEnrichment <- renderPlot({
    networkenrichmentforuse()
  })
  
  output$Downloadnetwork <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fnamenetwork)
      paste(pdf_file,".pdf", sep='')
    },
    content <- function(file) {
      pdf(file, height= 10, width=12,onefile = FALSE)
      networkplot()
      dev.off()
    }, contentType = 'image/pdf')
  
  output$Downloadceenrichment <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fnamenetwork)
      paste(pdf_file,"_enrichment.pdf", sep='')
    },
    content <- function(file) {
      pdf(file, height= 10, width=12,onefile = FALSE)
      networkenrichmentforuse()
      dev.off()
      #file.copy(paste(pdf_file,'.pdf', sep='') ,file, overwrite=TRUE)
    }, contentType = 'image/pdf')
  
  output$Downloadceenrichmentable <- downloadHandler(
    filename <- function() {
      pdf_file <<- as.character(input$fnamenetwork)
      paste0(pdf_file,'_enrichment.csv')
    },
    content <- function(file) {
      enrichment_table <- data.frame(networkenrichmentforuse())
      write.csv(enrichment_table,file,row.names = F)
    })
  
  ## ReadMe
  output$ReadMe1 <- renderUI({
    str00 <- paste("&emsp;")
    str0 <- paste("Example")
    str1 <- paste("&emsp; 1. Input data: GSE7621 or GSE136666.csv")
    str2 <- paste("&emsp; 2. Easily click the buttons from table panel to get corresponding results")
    str3 <- paste("&emsp; 3. NOTICE: Please check that your file types are the same as the ones you clicked on, e.g. all are Array types")

    str21 <- paste("Format of your files")
    str211 <- paste("Array Data")
    str22 <- paste("&emsp; 1. Upload your file: expression matrix, group list and annotation platform information")
    str23 <- paste("&emsp; 2. Row names of the expression matrix are probe names")
    str24 <- paste("&emsp; 3. Column names are sample lists")
    str212 <- paste("High-throughput RNA sequencing data")
    str222 <- paste("&emsp; 1. Upload your file: expression matrix, group list and annotation platform information")
    str232 <- paste("&emsp; 2. Row names of the expression matrix are Ensembl numbers")
    str242 <- paste("&emsp; 3. Column names are sample lists")
    
    str31 <- paste("Module")
    str32 <- paste("&emsp; 1. Heatmap plot &emsp;")
    str33 <- paste("&emsp; 2. Volcanol plot &emsp;")
    str34 <- paste("&emsp; 3. Enrichment &emsp;")
    str35 <- paste("&emsp; 4. RNA classification &emsp;")
    str36 <- paste("&emsp; 5. ceRNA network1: ceRNA network based on PCC &emsp;")
    str37 <- paste("&emsp; 6. ceRNA network2: ceRNA network based on databases &emsp;")
    str38 <- paste("&emsp; 7. ceRNA network3: ceRNA network based on PCC and databases &emsp;")

    HTML(paste(str00,h4(strong(str0)), str1, str2, str3,
               str00,h4(strong(str21)),strong(str211),str22,str23,str24,str00,strong(str212),str222,str232,str242,
               str00,h4(strong(str31)),str32,str33,str34,str35,str36,str37,str38,
               sep = '<br/>'))
  })
}


# Run the application 
shinyApp(ui = ui, server = server)
