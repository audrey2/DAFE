#' GSEAKEGG UI Function
#'
#' @description This module run Gene Set Enrichment Analysis on Kegg pathway.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
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
    library(cowplot)
    library(magick)
    library(ggimage)
    ns <- session$ns


    #Fonction renvoyant l'ui si la méthode choisi est GSEA
    output$condPanel<- renderUI({
      LL=c()

      if(inputParameter$Gsea_ORA == 'GSEA'){



        LL[[1]] = fluidRow(box(width = 12,status = 'info',solidHeader = TRUE,collapsible=TRUE,title=h1('GSEA KEGG settings',icon('gear')),
          
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
                   DT::dataTableOutput(ns('tableau'))%>% withSpinner()
                )
              )
            )
          ),
          tabPanel("Plot",
            fluidRow(
              column(width=12,
                box(width = NULL,status = 'success',solidHeader = TRUE,title=h1('KEGG pathway',icon('chart-simple')),
                    column(width=3,uiOutput(ns('pathIdd'))),
                    column(width=6,plotlyOutput(ns('png'))%>% withSpinner())
                )
              ),
              column(width=6,
                box(width = NULL,status = 'success',solidHeader = TRUE,title=h1('Dotplot',icon('chart-simple')),
                  uiOutput(ns('numberCategrorie')),
                  plotlyOutput(ns('dotplotKegg'))%>% withSpinner(),
                )
              ),
              column(width=6,
                box(width = NULL,status = 'success',solidHeader = TRUE,title=h1('KEGG GSEAplot',icon('chart-simple')),
                    uiOutput(ns('pathIdd2')),
                    plotOutput(ns('gseaplotKegg'))%>% withSpinner()
                )
              )
            ),
          )
        ))
      }
      return(LL)
    })

    # Fonction chargant la librairie correspondante à l'organisme choisi
    organism<- reactive({
      organism = inputParameter$orDb
     
      library(organism, character.only = TRUE)

      return(organism)
    })

    output$numberCategrorie <- renderUI({
       numericInput(ns("numberCat"), min=1,value=min(5,nrow(tabb())),max=nrow(tabb()),label="Write number of category to show",width="20%")

      })
    # Fonction renvoyant un objet gseKEGG
    gene_list_kegg <- eventReactive(input$go,{
      
      req(inputParameter$fileOr)
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

    # Fonction renvoyant les différent pathway significativement enrichi
    choix <- reactive({
      name = gene_list_kegg()$ID
      num  = c(1:length(name))
      choiceTable = data.frame(name, num)
      choix = setNames(as.numeric(choiceTable$num), choiceTable$name)

      return(choix)

    })

    # Fonction renvoyant à l'ui les différent pathway significativement enrichi sous la forme d'un select input
    output$pathIdd <- renderUI({

      selectInput(ns('pathiId'),"Select Path Id to observe" ,choices = choix())
    })
      output$pathIdd2 <- renderUI({

      selectInput(ns('pathiId2'),"Select Path Id to observe" ,choices = choix(),width="20%")
    })

    shinyDirChoose(input, 'folder', root=c(root='~'), filetypes=c('png', 'txt'))


    # Fonction renvoyant la liste des gènes triés par log2FC decroissant
    kge <- eventReactive(input$go,{
      req(inputParameter$fileOr)
      df = read.csv(inputParameter$fileOr[1,'datapath'], header=TRUE,sep=inputParameter$sep)
      
      if(inputParameter$gsea_order=="log2FC") {kegg_gene_list <- df$log2FC}
      if(inputParameter$gsea_order=="pval")   {kegg_gene_list <- df$padj}
      if(inputParameter$gsea_order=="Stat")   {kegg_gene_list <- df$stat}
      
      df$X=df[,inputParameter$row.names]
      names(kegg_gene_list)=df$X

      organism=organism()
      ids<-bitr(names(kegg_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=organism)
      dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]
      df2 = df[df$X %in% dedup_ids$ENSEMBL,]
      df2$Y = dedup_ids$ENTREZID

      if(inputParameter$gsea_order=="log2FC")  {kegg_gene_list <- df2$log2FC}
      if(inputParameter$gsea_order=="pval")    {kegg_gene_list <- df2$padj}
      if(inputParameter$gsea_order=="Stat")    {kegg_gene_list <- df2$stat}

      names(kegg_gene_list) <- df2$Y
      kegg_gene_list<-na.omit(kegg_gene_list)
      kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

      return(kegg_gene_list)
    })

    PATHID<- reactive(input$pathiId)
    PATHID2<- reactive(input$pathiId2)


### PLOT
   
    output$dotplotKegg <- renderPlotly(
      ggplotly( 
        dotplot(gene_list_kegg(), showCategory = input$numberCat, title = "Enriched Pathways" )
        )
    )
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


   # Fonction renvoyant le nom de l'image correspondant au patway KEGG
    pathView<- reactive({

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

    # Fonction renvoyant le tableau de resultat de l'enrichissement
    tabb <- reactive({
      table= cbind(gene_list_kegg()$ID,gene_list_kegg()$Description,gene_list_kegg()$enrichmentScore,gene_list_kegg()$NES,gene_list_kegg()$pvalue,gene_list_kegg()$p.adjust,gene_list_kegg()$qvalue,gene_list_kegg()$rank)
      colnames(table)= c('ID','Description','enrichmentScore','NES','pvalue','p.adjust','qvalues','rank')
      return(data.frame(table))
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
  })
}

## To be copied in the UI
# mod_GSEAKEGG_ui("GSEAKEGG_1")

## To be copied in the server
# mod_GSEAKEGG_server("GSEAKEGG_1")
