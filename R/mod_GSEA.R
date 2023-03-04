#' GSEA UI Function
#'
#' @description This module run Gen Set Enrichment Analysis on GO terms.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_GSEA_ui <- function(id){
  ns <- NS(id)
  tagList(
    uiOutput(ns('condPanel'))
  )
}

#' GSEA Server Functions
#'
#' @noRd
mod_GSEA_server <- function(id,inputParameter){
  library(DOSE)
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(GO.db)
  library(plotly)
  library(shinycssloaders)


  moduleServer( id, function(input, output, session){

    ns <- session$ns

    #Fonction renvoyant l'ui si la méthode choisi est GSEA
    output$condPanel<- renderUI({
      LL=c()
      if(inputParameter$Gsea_ORA =='GSEA'  ){


        LL[[1]]=fluidRow(box(width = 12,status = 'info',solidHeader = TRUE,collapsible=TRUE,title=h1('GSEA GO terms Settings',icon('gear')),
      
          fluidRow(
            column(width=1),
            column(width=3,
              selectInput(ns("onto"),"Choose GO annotation", choices=c("Biological Process"='BP', "Cellular Component"="CC","Molecular Function"="MF","ALL"="ALL" ))),
            column(width=3,
              sliderInput(ns("pvalCutOff"),label = "p-Value cutoff from output ",
                      min = 0, max =0.1, value = 0.05,step=0.01)
            ),
            column(width=3,
              sliderInput(ns("GoLevel"),label="Select Go Level",min=1,max=7,value=4,step=1)
            ),
            column(width=2,br(),br(),br()
              ,actionButton(ns("go"),"Start",class="buttS",icon("play")))
          )
        ))

        LL[[2]]=fluidRow(tabBox(width=12,height=NULL,
          tabPanel("Table",
            fluidRow(
              column(width=12,
                box(width = NULL,status = 'primary',solidHeader = TRUE,title=h1('Table of results',icon('table')),
                  DT::dataTableOutput(ns('tableau'))%>% withSpinner()
                )
              )
            )
          ),
          tabPanel("Plot",
            fluidRow(
              column(width=6,
                box(width = NULL,status = 'success',solidHeader = TRUE,title=h1('Dotplot',icon('chart-simple')),
                  uiOutput(ns('numberCategory')),
                  plotlyOutput(ns('dotplot'))%>% withSpinner()
                )
              ),
              column(width=6,
                box(width = NULL,status = 'success',solidHeader = TRUE,title=h1('GSEAplot',icon('chart-simple')),
                  uiOutput(ns('pathId')),
                  plotOutput(ns('gseaplot'))%>% withSpinner()
                )
              )
            )
          )
        ))
      }
      return(LL)  
    })
    # Fonction renvoyant à l'ui un slideerInput pour deteminer le nombre de categrories a afficher sur le dotPlot
    output$numberCategory<-renderUI({
        numericInput(ns("numberCat"), min=1,value=min(5,nrow(tabb())),max=nrow(tabb()),label="Write number of category to show",width="20%")
    })

    # Fonction renvoyant à l'ui les différent GO terms significativement enrichi sous la forme d'un select input
    output$pathId <- renderUI({

      L= selectInput(ns('pathId'),"Select Pathway Id to observe" ,choices = choix(),width="20%")
      return(L)

    })

# FUNCTION

    # Fonction renvoyant un objet GSEA enfonction des spécificité choisi
    geneList<-eventReactive(input$go,{
  
      req(inputParameter$fileOr,inputParameter$orDb,inputParameter$gsea_order,input$GoLevel)

      BiocManager::install(inputParameter$orDb,update=FALSE)
      library(inputParameter$orDb, character.only = TRUE)

      df = read.csv(inputParameter$fileOr[1,'datapath'], header=TRUE,sep=inputParameter$sep)
      organism = inputParameter$orDb
     



      if(inputParameter$gsea_order=="log2FoldChange") {original_gene_list <- df$log2FoldChange}
      if(inputParameter$gsea_order=="pval")   {original_gene_list <- df$padj}
      if(inputParameter$gsea_order=="Stat")   {original_gene_list <- df$stat}     

      row.names(df)=df[,inputParameter$row.names]
      names(original_gene_list) <- row.names(df)
      gene_list<-na.omit(original_gene_list)
      gene_list = sort(gene_list, decreasing = TRUE)

      gse <- gseGO(geneList=gene_list,
                   ont =input$onto,
                   keyType = 'ENSEMBL',
                   minGSSize = 3,
                   maxGSSize = 800,
                   pvalueCutoff =input$pvalCutOff,
                   verbose = FALSE,
                   OrgDb = organism,
                   pAdjustMethod = inputParameter$pAdjMethod)

      gse=goFilter(gse,input$GoLevel)

      return(gse)

    })

    # Fonction renvoyant les différent GO terms significativement enrichi
    choix <- reactive({
      name=geneList()$ID
      num  = c(1:length(name))
      choiceTable = data.frame(name, num)
      choix = setNames(as.numeric(choiceTable$num), choiceTable$name)
      return(choix)

    })

    
# OUTPUT
## TABLE
# Fonction renvoyant le tableau d'enrcihissement
    tabb <- reactive({
     
        req(inputParameter)
     

      table=geneList()

      table= cbind(table$ID,table$Description,table$enrichmentScore,table$NES,table$pvalue,table$p.adjust,table$qvalue)
      colnames(table)= c('ID','Description','enrichmentScore','NES','pvalue','p.adjust','qvalues')
      return(data.frame(table))


    })


    # Fonction renvoyant sous la forme d'un Datatable le tableau d'enrcihissement
    output$tableau <- DT::renderDataTable(DT::datatable(
      { tabb()},
      extensions = 'Buttons',
      caption="Table: GSEA table of GO terms",
      options = list(

        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel')
      ),

      class = "display"
    ))

## PLOT
    # Fonction renvoyant un dotplot
    output$dotplot <- renderPlotly({
      req(input$numberCat)
      withProgress(message = "Plotting dotPlot ...",{
        
        require(DOSE)
        ggplotly(dotplot(geneList(), showCategory=input$numberCat, split=".sign") + facet_grid(.~.sign))})

    })

    # Fonction renvoyant le GSEA plot
    output$gseaplot <- renderPlot({
      req(input$pathId)

      require(DOSE)
      enrichplot::gseaplot2(geneList(), title = geneList()$Description[as.numeric(input$pathId)], geneSetID = as.numeric(input$pathId))

    })

  })
}