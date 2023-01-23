#' ORA UI Function
#'
#' @description A shiny Module.
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


    #Fonction renvoyant l'ui si la méthode choisi est ORA
    output$condPanel<- renderUI({
      LL=c()

      if(inputParameter$Gsea_ORA=='ORA' ){
        LL[[1]]=fluidRow(box(class='VIOLET',width = 12,

              h1('ORA GO terms Settings',icon('cogs')),
              hr(),
              fluidRow(
              column(width=3,
                selectInput(ns("onto"),"Choose GO annotation", choices=c("Biological Process"='BP', "Cellular Component"="CC","Molecular Function"="MF","ALL"="ALL" ))),
              column(width=3,
                sliderInput(ns("pvalCutOff"),label = "p-Value cutoff from input ",
                            min = 0, max =0.1, value = 0.05,step=0.01)),
              column(width=3,
                sliderInput(ns("GoLevel"),label="Select Go Level",min=1,max=7,value=4,step=1)),
              column(width=3,
                selectInput(ns("DEG")," Choose element",choices=c('both','over','under'))),)
        ))


      LL[[3]]=fluidRow(tabsetPanel(
        tabPanel("Table",
                 DT::dataTableOutput(ns('tableau'))%>% withSpinner()
        ),tabPanel("Plot",
                  uiOutput(ns('pathId')),h3('Dotplot'),  plotlyOutput(ns('dotplot'))%>% withSpinner(),h3('Goplot'),plotlyOutput(ns('goplot'))%>% withSpinner(),plotlyOutput(ns('barplot'))%>% withSpinner()
        ))
      )
}
      return(LL)}
    )
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

    # Fonction renvoyant un objet ORA enfonction des spécificité choisi
    geneList<-function(){

      req(inputParameter$fileOr)

      df=read.csv(inputParameter$fileOr[1,'datapath'], header=inputParameter$header,sep=inputParameter$sep)
      if(input$DEG=='over'){df=subset(df,log2FoldChange >0)}
      if(input$DEG=='under'){df=subset(df,log2FoldChange <0)}

      organism=inputParameter$orDb
      BiocManager::install(organism, character.only = TRUE)
      library(organism, character.only = TRUE)

      original_gene_list <- df$log2FoldChange
      names(original_gene_list) <- df$X
      gene_list<-na.omit(original_gene_list)
      gene_list = sort(gene_list, decreasing = TRUE)
      sig_genes_df = subset(df, padj < 0.05)
      genes <- sig_genes_df$log2FoldChange
      names(genes) <- sig_genes_df$X
      genes <- na.omit(genes)
      genes <- names(genes)[abs(genes) > 2]
      go_enrich <- enrichGO(gene = genes,
                            universe = names(gene_list),
                            OrgDb = organism,
                            keyType = 'ENSEMBL',
                            readable = T,
                            ont = "BP",
                            pvalueCutoff = 0.05,
                            qvalueCutoff = 0.10)



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

    }

    # Fonction renvoyant le tableau de resultat de la methode d'enrichissement
    tabb <- reactive({

      table=cbind(gg()$ID,gg()$Description,gg()$GeneRatio,gg()$BgRatio,gg()$pvalue,gg()$p.adjust,gg()$qvalue,gg()$geneID)
      colnames(table)=cbind('ID','Description','GeneRatio','BgRatio','pvalue','p.adjust','qvalue','geneID')
      return( data.frame(table))
    })

    # Fonction renvoyant à l'ui le tableau de résultat d'enrichissement
    output$tableau <- DT::renderDataTable(DT::datatable(
    {
        tabb()},
      extensions = 'Buttons',
      caption="Table: ORA table of GO terms",
      options = list(

        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel')
      ),

      class = "display"
    ))


    gg <- reactive({

      return(geneList())
    })

    # Fonction renvoyant un dotplot
    output$dotplot <- renderPlotly({

      gbb=gg()
      require(DOSE)
     ggplotly( dotplot(gbb, showCategory=10))

    })

    # Fonction renvoyant un barplot
    output$barplot <- renderPlotly({

      gbb=gg()
      ggplotly(barplot(gbb,
              drop = TRUE,
              showCategory = 10,
              title = "GO Biological Pathways",
              font.size = 8))

    })

    # Fonction renvoyant un goplot
    output$goplot <- renderPlotly({

      gbb=gg()
      require(DOSE)
      ggplotly(goplot(gbb))

    })

  })
}

## To be copied in the UI
# mod_ORA_ui("ORA_1")

## To be copied in the server
# mod_ORA_server("ORA_1")
