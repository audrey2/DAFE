#' GSEA UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_GSEA_ui <- function(id){
  ns <- NS(id)
  tagList(


  hr(),
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


        LL[[1]]=fluidRow(box(width = 12,status = 'info',solidHeader = TRUE,collapsible=TRUE,title=h1('GSEA GO terms settings',icon('gear')),
          hr(),
          fluidRow(
            column(width=3,
          selectInput(ns("onto"),"Choose GO annotation", choices=c("Biological Process"='BP', "Cellular Component"="CC","Molecular Function"="MF","ALL"="ALL" ))),
          column(width=3,
          sliderInput(ns("pvalCutOff"),label = "p-Value cutoff from input ",
                      min = 0, max =0.1, value = 0.05,step=0.01)),
          column(width=3,
          sliderInput(ns("GoLevel"),label="Select Go Level",min=1,max=7,value=4,step=1)),
        )
        ))


      LL[[3]]=fluidRow(tabBox(width=12,
        tabPanel("Table",
                 box(width = 12,status = 'primary',solidHeader = TRUE,title=h1('Table of results',icon('table')),
                 DT::dataTableOutput(ns('tableau'))%>% withSpinner())
        ),tabPanel("Plot", box(width = 6,status = 'success',solidHeader = TRUE,title=h1('Dotplot',icon('char-simple')),
                     plotlyOutput(ns('dotplot'))%>% withSpinner()),

                   box(width = 6,status = 'success',solidHeader = TRUE,title=h1('Gseaplot',icon('char-simple')),

                   uiOutput(ns('pathId')),plotOutput(ns('gseaplot'))%>% withSpinner()
        )),
        ))



}
      return(LL)}
    )

    # Fonction chargant la librairie correspondante à l'organisme choisi
    organism<- reactive({
      organism = inputParameter$orDb
     
      library(organism, character.only = TRUE)

      return(organism)
    })



    # Fonction renvoyant le tableau charger
    df <-reactive({

      req(inputParameter$fileOr)


      df = read.csv(inputParameter$fileOr[1,'datapath'], header=TRUE,sep=inputParameter$sep)

      return(df)

    })


    # Fonction renvoyant un objet GSEA enfonction des spécificité choisi
    geneList<-function(){
      df=df()

      organism=organism()


      if(inputParameter$gsea_order=="log2FC")
      {original_gene_list <- df$log2FC}
      if(inputParameter$gsea_order=="pval")
      {original_gene_list <- df$padj}
      if(inputParameter$gsea_order=="Stat")
      {original_gene_list <- df$stat}
      print(inputParameter$row.names)
      print(df[,inputParameter$row.names])
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
                   verbose = TRUE,
                   OrgDb = organism,
                   pAdjustMethod = inputParameter$pAdjMethod)

print('aaaaaaaaaaaaa')
      #gse=simplify(
       # gse,
        #cutoff = 0.7,
       # by = "p.adjust",
        #select_fun = min,
        #measure = "Wang",
        #semData = NULL
      #)
print('bbbbbbbbbbbbbbbbb')
      gse=goFilter(gse,input$GoLevel)
    print('ccccccccccccccccccccc')
      return(gse)

    }

    # Fonction renvoyant les différent GO terms significativement enrichi
    choix <- reactive({
      name=geneList()$ID
      num  = c(1:length(name))
      choiceTable = data.frame(name, num)
      choix = setNames(as.numeric(choiceTable$num), choiceTable$name)
      return(choix)

    })

    # Fonction renvoyant à l'ui les différent GO terms significativement enrichi sous la forme d'un select input
    output$pathId <- renderUI({

      L= selectInput(ns('pathId'),"Select Path Id to observe" ,choices = choix())
      return(L)

    })



    gg <- reactive({

      return(geneList())
    })
    # Fonction renvoyant un dotplot
    output$dotplot <- renderPlotly({

      withProgress(message = "Plotting dotPlot ...",{
        gbb=gg()
        require(DOSE)
        ggplotly(dotplot(gbb, showCategory=10, split=".sign") + facet_grid(.~.sign))})

    })

    # Fonction renvoyant le GSEA plot
    output$gseaplot <- renderPlot({

      gbb=geneList()
      require(DOSE)
      enrichplot::gseaplot2(gbb, title = gbb$Description[as.numeric(input$pathId)], geneSetID = as.numeric(input$pathId))

    })

    # Fonction renvoyant le tableau d'enrcihissement
    tabb <- reactive({
      table=gg()
      print(head(table))
      table= cbind(table$ID,table$Description,table$enrichmentScore,table$NES,table$pvalue,table$p.adjust,table$qvalue,table$rank)
      colnames(table)= c('ID','Description','enrichmentScore','NES','pvalue','p.adjust','qvalues','rank')
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





  return(input)

  })
}

## To be copied in the UI
# mod_GSEA_ui("GSEA_1")

## To be copied in the server
# mod_GSEA_server("GSEA_1")
