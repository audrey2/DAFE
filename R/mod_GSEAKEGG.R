#' GSEAKEGG UI Function
#'
#' @description This module run Gene Set Enrichment Analysis on Kegg pathway.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @author Audrey BEAUFILS
#'
#' @importFrom shiny NS tagList
mod_GSEAKEGG_ui <- function(id){
  library(shinyFiles)
  ns <- NS(id)
  tagList(
    uiOutput(ns('condPanel'))
  )
}

#' GSEAKEGG Server Functions
#'
#' @noRd
mod_GSEAKEGG_server <- function(id,inputParameter){
  moduleServer( id, function(input, output, session){
    library(clusterProfiler)
    library(enrichplot)
    library(pathview)
    library(shinyFiles)
    library(png)
    library(plotly)


    ns <- session$ns


# UI OUTPUT
    #Fonction renvoyant l'ui si la méthode choisi est GSEA
    output$condPanel<- renderUI({
      LL=c()

      if(inputParameter$Gsea_ORA == 'GSEA'){



        LL[[1]] = fluidRow(box(width = 12,status = 'info',solidHeader = TRUE,collapsible=TRUE,title=h1('GSEA KEGG Settings',icon('gear')),
          
          fluidRow(
            column(width=1),
            column(width=3,
              sliderInput(ns("pvalCutOff"),label = "p-Value cutoff from output ",
                     min = 0, max =0.1, value = 0.05,step=0.01),
            ),
            column(width=3,
              br(),br(),
              shinyDirButton(ns('folder'), 'Select a folder to save pathway plot', 'Please select a folder', FALSE)                  
            ),
            column(width=3,
              br(),br(),br(),
              actionButton(ns("go"),"Start",class="buttS",icon("play"))
            )
          )
        ))
        LL[[2]] = fluidRow(tabBox(width=12,height=NULL,
          tabPanel("Table",
            fluidRow(
              column(width=12,
                box(width = NULL,status = 'primary',solidHeader = TRUE,title=h1('Table of results',icon('table')),
                   DT::dataTableOutput(ns('tableau'))%>% withSpinner(color="#605CA8")
                )
              )
            )
          ),
          tabPanel("Plot",
            fluidRow(
              column(width=12,
                box(width = NULL,status = 'success',solidHeader = TRUE,title=h1('KEGG pathway',icon('chart-simple')),
                    column(width=3,uiOutput(ns('pathIdd'))),
                    column(width=6,plotlyOutput(ns('png'))%>% withSpinner(color="#CDCDE6"))
                )
              ),
              column(width=6,
                box(width = NULL,status = 'success',solidHeader = TRUE,title=h1('Dotplot',icon('chart-simple')),
                  uiOutput(ns('numberCategrorie')),
                  plotlyOutput(ns('dotplotKegg'))%>% withSpinner(color="#CDCDE6"),
                )
              ),
              column(width=6,
                box(width = NULL,status = 'success',solidHeader = TRUE,title=h1('KEGG GSEAplot',icon('chart-simple')),
                    uiOutput(ns('pathIdd2')),
                    plotOutput(ns('gseaplotKegg'))%>% withSpinner(color="#CDCDE6")
                )
              )
            ),
          )
        ))
      }
      return(LL)
    })

    # Fonction renvoyant à l'ui un slideerInput pour deteminer le nombre de categrories a afficher sur le dotPlot
    output$numberCategrorie <- renderUI({
       numericInput(ns("numberCat"), min=1,value=min(5,nrow(tabb())),max=nrow(tabb()),label="Write number of category to show",width="20%")

    })

    # Fonction renvoyant à l'ui les différent pathway significativement enrichi sous la forme d'un select input pour le pathway KEGG
    output$pathIdd <- renderUI({

      selectInput(ns('pathiId'),"Select Path Id to observe" ,choices = choix())
    })
    # Fonction renvoyant à l'ui les différent pathway significativement enrichi sous la forme d'un select input pour le gseaKEGG
    output$pathIdd2 <- renderUI({

      selectInput(ns('pathiId2'),"Select Path Id to observe" ,choices = choix(),width="20%")
    })

    shinyDirChoose(input, 'folder', root=c(root='~'), filetypes=c('png', 'txt'))


# FUNCTION

    # Fonction renvoyant un objet gseKEGG
    gene_list_kegg <- eventReactive(input$go,{
      
      req(inputParameter$fileOr)
      req(inputParameter$orDb)
      BiocManager::install(inputParameter$orDb,update=FALSE)
      library(inputParameter$orDb, character.only = TRUE)

      kegg_gene_list = kge()

      kegg_organism = switch(inputParameter$orDb,'org.Dm.eg.db'='dme','org.At.tair.db'='ath',
                                'org.Dr.eg.db'='dre','org.Ss.eg.db'='ssc','org.Cf.eg.db'='cfa',
                             'org.Ag.eg.db'='aga','org.EcSakai.eg.db'='ecs','org.Hs.eg.db'='hsa',
                             'org.Sc.sgd.db'='sce','org.Ce.eg.db'='sko','org.Gg.eg.db'='gga',
                             'org.EcK12.eg.db'='eco','org.Pt.eg.db'='ptr','org.Mxanthus.db'='mxa',
                             'org.Mm.eg.db'='mmu','org.Rn.eg.db'='rno','org.Bt.eg.db'='bta',
                             'org.Mmu.eg.db'='mvv','org.Xl.eg.db'='xla','org.Pf.plasmo.db'='pfa')


       gseaKegg <- gseKEGG(geneList     = kegg_gene_list,
                     organism     = kegg_organism,
                     nPerm        = 10000,
                     minGSSize    = 3,
                     maxGSSize    = 800,
                     pvalueCutoff = input$pvalCutOff,
                     pAdjustMethod = inputParameter$pAdjMethod,
                     keyType       = "ncbi-geneid")


      return(gseaKegg)
    })

   
    # Fonction renvoyant la liste des gènes triés par log2FoldChange decroissant
    kge <- eventReactive(input$go,{
      req(inputParameter$fileOr)
      df = read.csv(inputParameter$fileOr[1,'datapath'], header=TRUE,sep=inputParameter$sep)
      
      if(inputParameter$gsea_order=="log2FoldChange") {kegg_gene_list <- df$log2FoldChange}
      if(inputParameter$gsea_order=="pval")   {kegg_gene_list <- df$padj}
      if(inputParameter$gsea_order=="Stat")   {kegg_gene_list <- df$stat}
      
      df$X=df[,inputParameter$row.names]
      names(kegg_gene_list)=df$X
 

      ids<-bitr(names(kegg_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=inputParameter$orDb)
      dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]
      df2 = df[df$X %in% dedup_ids$ENSEMBL,]
      df2$Y = dedup_ids$ENTREZID

      if(inputParameter$gsea_order=="log2FoldChange")  {kegg_gene_list <- df2$log2FoldChange}
      if(inputParameter$gsea_order=="pval")    {kegg_gene_list <- df2$padj}
      if(inputParameter$gsea_order=="Stat")    {kegg_gene_list <- df2$stat}

      names(kegg_gene_list) <- df2$Y
      kegg_gene_list<-na.omit(kegg_gene_list)
      kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

      return(kegg_gene_list)
    })

    PATHID<- reactive({
      req(input$pathiId) 
      input$pathiId})
    PATHID2<- reactive({
      req(input$pathiId2) 
      input$pathiId2})

     # Fonction renvoyant les différent pathway significativement enrichi
    choix <- reactive({
      name = gene_list_kegg()$ID
      num  = c(1:length(name))
      choiceTable = data.frame(name, num)
      choix = setNames(as.numeric(choiceTable$num), choiceTable$name)

      return(choix)

    })

    # Fonction renvoyant le nom de l'image correspondant au patway KEGG
    pathView<- reactive({
      req(inputParameter$fileOr,input$folder)
      dir="~"
      for( i in 1: length(input$folder)){

        dir= paste0(dir,input$folder$path[i],"/")
      }
      setwd(dir)
      kegg_organism = switch(inputParameter$orDb,'org.Dm.eg.db'='dme','org.At.tair.db'='ath',
                           'org.Dr.eg.db'='dre','org.Ss.eg.db'='ssc','org.Cf.eg.db'='cfa',
                           'org.Ag.eg.db'='aga','org.EcSakai.eg.db'='ecs','org.Hs.eg.db'='hsa',
                           'org.Sc.sgd.db'='sce','org.Ce.eg.db'='sko','org.Gg.eg.db'='gga',
                           'org.EcK12.eg.db'='eco','org.Pt.eg.db'='ptr','org.Mxanthus.db'='mxa',
                           'org.Mm.eg.db'='mmu','org.Rn.eg.db'='rno','org.Bt.eg.db'='bta',
                           'org.Mmu.eg.db'='mvv','org.Xl.eg.db'='xla','org.Pf.plasmo.db'='pfa')

    
      pathview(gene.data=kge(), pathway.id=names(choix()[as.numeric(PATHID())]), species = kegg_organism, kegg.dir = dir)

      dir=paste0(dir,names(choix()[as.numeric(PATHID())]))

      dir=paste0(dir,".pathview.png")
      return(dir)

     })


# OUTPUT
## TABLE
    # Fonction renvoyant le tableau de resultat de l'enrichissement
    tabb <- reactive({
      
      if(nrow(gene_list_kegg())==0){
      table= cbind(gene_list_kegg()$ID,gene_list_kegg()$Description,gene_list_kegg()$enrichmentScore,gene_list_kegg()$NES,gene_list_kegg()$pvalue,gene_list_kegg()$p.adjust,gene_list_kegg()$qvalue)
      colnames(table)= c('ID','Description','enrichmentScore','NES','pvalue','p.adjust','qvalues')
      return(data.frame(table))}
      return(data.frame())
    })


    # Fonction renvoyant à l'ui le datatable des resultats
    output$tableau <- DT::renderDataTable(
      DT::datatable({
      tabb()},
        extensions = 'Buttons',
        caption="Table: GSEA table of KEGG Pathway",
        options = list(
            dom = 'Bfrtip',
            buttons = c('copy', 'csv', 'excel')
        ),  class = "display"
      )
    )

## PLOT
   
    output$dotplotKegg <- renderPlotly({
      req(input$numberCat)
      ggplotly( 
        dotplot(gene_list_kegg(), showCategory = input$numberCat, title = "Enriched Pathways" )
        )
    })
    output$gseaplotKegg <- renderPlot( 
      enrichplot::gseaplot2(gene_list_kegg(),  title = gene_list_kegg()$Description[as.numeric(PATHID2())],
          geneSetID =as.numeric( PATHID())
      )
    )
    output$png<- renderPlotly( {

        fig=plot_ly()
      fig = fig %>% add_trace(type="image", source= raster2uri(readPNG(pathView())), hoverinfo = "skip")
      return(fig)
    })

  })
}