#' ORA UI Function
#'
#' @description This module run Over Representation Analysis on GO terms.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_ORA_ui <- function(id){
  ns <- NS(id)
  tagList(

            uiOutput(ns('condPanel')),
  )
}

#' ORA Server Functions
#'
#' @noRd
mod_ORA_server <- function(id,inputParameter){

  moduleServer( id, function(input, output, session){
    library(shinycssloaders)
    library(plotly)

    ns <- session$ns

# UI OUTPUT

    #Fonction renvoyant l'ui si la méthode choisi est ORA
    output$condPanel<- renderUI({
      LL=c()

      if(inputParameter$Gsea_ORA=='ORA' ){
        LL[[1]]=fluidRow(box(width = 12,status = 'info',solidHeader = TRUE,collapsible=TRUE,title=h1('ORA GO terms Settings',icon('gear')),          
          fluidRow(
            column(width=1),
            column(width=3,
              sliderInput(ns("pvalI"),label = "p-Value cutoff from input ",
                          min = 0, max =0.1, value = 0.05,step=0.01),
              sliderInput(ns("log2I"),label = "Log2 Fold change cutoff from input ",
                          min = 0, max =5, value = 2,step=0.5)
            ),
            column(width=3,
              selectInput(ns("onto"),"Choose GO annotation", choices=c("Biological Process"='BP', "Cellular Component"="CC","Molecular Function"="MF","ALL"="ALL" )),
              br(),sliderInput(ns("GoLevel"),label="Select Go Level",min=1,max=7,value=4,step=1)
            ),
            column(width=3,
              sliderInput(ns("pvalCutOff"),label = "p-Value cutoff from output ",
                      min = 0, max =0.1, value = 0.05,step=0.01)
            ),
            column(width=2,
              br(),br(),br(),br(),br(),br(),br(),br(),br(),
              actionButton(ns("go"),"Start",class="buttS",icon("play"))
            )
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
          tabPanel("Plot",fluidRow(
            column(width=6,
              box(width = NULL,status = 'success',solidHeader = TRUE,title=h1('Dotplot',icon('chart-simple')),
                uiOutput(ns('numberCategrory')),
                plotlyOutput(ns('dotplot'))%>% withSpinner()
              )
            ),
            column(width=6,
              box(width = NULL,status = 'success',solidHeader = TRUE,title=h1('Barplot',icon('chart-simple')),
                uiOutput(ns('numberCategrory2')),
                plotlyOutput(ns('barplot'))%>% withSpinner()
              )
            ),
            column(width=6,
              box(width = NULL,status = 'success',solidHeader = TRUE,title=h1('Goplot',icon('chart-simple')),
                plotlyOutput(ns('goplot'))%>% withSpinner()
              )
            )
          ))
        ))
      }
      return(LL)
    })
    

    # Fonction renvoyant à l'ui un slideerInput pour deteminer le nombre de categrories a afficher sur le dotPlot
    output$numberCategrory <- renderUI({
        numericInput(ns("numberCat"), min=1,value=min(5,nrow(tabb())),max=nrow(tabb()),label="Write number of category to show",width="20%")
    })

    # Fonction renvoyant à l'ui un slideerInput pour deteminer le nombre de categrories a afficher sur le barPlot
    output$numberCategrory2 <- renderUI({
        numericInput(ns("numberCat2"), min=1,value=min(5,nrow(tabb())),max=nrow(tabb()),label="Write number of category to show",width="20%")
    })


# FUNCTION

    # Fonction renvoyant un objet ORA enfonction des spécificité choisi
    geneList<-eventReactive(input$go,{

      req(inputParameter$fileOr)
      req(inputParameter$orDb)
      BiocManager::install(inputParameter$orDb,update=FALSE)
      library(inputParameter$orDb, character.only = TRUE)

      df = read.csv(inputParameter$fileOr[1,'datapath'], header=TRUE,sep=inputParameter$sep)
      organism=inputParameter$orDb
      if(inputParameter$ora_order=="overexpressed")   {df = subset(df, log2FoldChange > 0)}
      if(inputParameter$ora_order=="underexpressed")  {df = subset(df, log2FoldChange < 0)}
   
     
      row.names(df) = df[,inputParameter$row.names]
      sig_genes_df = subset(df, padj < input$pvalI)
      genes = sig_genes_df$log2FoldChange
      names(genes) = sig_genes_df[,inputParameter$row.names]
      genes = na.omit(genes)
      genes = names(genes)[abs(genes) > input$log2I]

      go_enrich = enrichGO(gene = genes,
                            OrgDb = organism,
                            keyType = inputParameter$keytype,
                            readable = T,
                            ont = input$onto,
                            pvalueCutoff = input$pvalCutOff,
      )

      go_enrich=simplify(
        go_enrich,
        cutoff = 0.7,
        by = "p.adjust",
        select_fun = min,
        measure = "Wang",
        semData = NULL
      )

      go_enrich=gofilter(go_enrich,input$GoLevel)

      return(go_enrich)

    })

    # Fonction renvoyant le tableau de resultat de la methode d'enrichissement
    tabb <- reactive({

      table=cbind(geneList()$ID,geneList()$Description,geneList()$GeneRatio,geneList()$BgRatio,geneList()$pvalue,geneList()$p.adjust,geneList()$qvalue,geneList()$geneID)
      colnames(table)=cbind('ID','Description','GeneRatio','BgRatio','pvalue','p.adjust','qvalue','geneID')
      return( data.frame(table))
    })

# OUTPUT
## TABLE
    # Fonction renvoyant à l'ui le tableau de résultat d'enrichissement
    output$tableau <- DT::renderDataTable(
      DT::datatable({
        tabb()
      },
      extensions = 'Buttons',
      caption="Table: ORA table of GO terms",
      options = list(

        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel')
      ),  class = "display")
    )
## PLOT
    # Fonction renvoyant un dotplot
    output$dotplot <- renderPlotly({

      geneList = geneList()
      require(DOSE)
     ggplotly( dotplot(geneList, showCategory=input$numberCat))

    })

    # Fonction renvoyant un barplot
    output$barplot <- renderPlotly({

      geneList=geneList()
      ggplotly(barplot(geneList,
              drop = TRUE,
              showCategory = input$numberCat2,
              title = "GO Pathways",
              font.size = 8))

    })

    # Fonction renvoyant un goplot
    output$goplot <- renderPlotly({

      geneList=geneList()
      require(DOSE)
      ggplotly(goplot(geneList))

    })

  })
}