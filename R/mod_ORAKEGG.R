#' ORAKEGG UI Function
#'
#' @description This module run Over Representation Analysis on KEGG pathway.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_ORAKEGG_ui <- function(id){
  ns <- NS(id)
  tagList(

    uiOutput(ns('condPanel'))
  )
}

# ('RAKEGG Server Functions
#'
#' @noRd
mod_ORAKEGG_server <- function(id,inputParameter){
  library(shinycssloaders)
  moduleServer( id, function(input, output, session){
    library(clusterProfiler)
    library(enrichplot)
    library(pathview)
    library(shinyFiles)
    library(png)
    library(plotly)

    ns <- session$ns

# UI OUTPUT
    #Fonction renvoyant l'ui si la méthode choisi est ORA
    output$condPanel<- renderUI({
      LL=c()

      if(inputParameter$Gsea_ORA == 'ORA'){

        LL[[1]]=fluidRow(
          box(width = 12,status = 'info',solidHeader = TRUE,collapsible=TRUE,title=h1('ORA KEGG Settings',icon('gear')),
            fluidRow(
              column(width=1),
              column(width=3,
                sliderInput(ns("pvalI"),label = "p-Value cutoff from input ",
                          min = 0, max =0.1, value = 0.05,step=0.01),
                sliderInput(ns("log2I"),label="Log2 fold change cutoff from input ",min=0,max=5,value=1,step=0.1),
              ),
              column(width=3,
                sliderInput(ns("pvalCutOff"),label = "p-Value cutoff from output ",
                          min = 0, max =0.1, value = 0.05,step=0.01),
              ),
              column(width=3,
                br(),br(),
                 shinyDirButton(ns('folder'), 'Select a folder to save pathway plot', 'Please select a folder', FALSE),
              ),
              column(width=2,
                 br(),br(),br(),br(),br(),br(),br(),br(),br(),
                 actionButton(ns("go"),"Start",class="buttS",icon("play"))
              )
            )
          )
        )
        LL[[2]]= br()

        LL[[3]]=fluidRow(
          tabBox(width=12,height=NULL,
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
                    uiOutput(ns('pathIdd')),
                    plotlyOutput(ns('png'))%>% withSpinner()
                  )
                ),
                column(width=6,
                  box(width = NULL,status = 'success',solidHeader = TRUE,title=h1('Dotplot',icon('chart-simple')),
                    uiOutput(ns('numberCategrory')),
                    plotlyOutput(ns('dotplotK'))%>% withSpinner(),
                  )
                ),
                column(width=6,
                  box(width = NULL,status = 'success',solidHeader = TRUE,title=h1('Barplot',icon('chart-simple')),
                      uiOutput(ns('numberCategrory2')),
                      plotlyOutput(ns('barplotK'))%>% withSpinner()
                  )
                )
              )
            )
          )
        )
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

    # Fonction renvoyant à l'ui les différent pathway significativement enrichi sous la forme d'un select input

    output$pathIdd <- renderUI({

      L= selectInput(ns('pathiId'),"Select Path Id to observe" ,choices = choix(),width="20%")
      return(L)
    })

    PATHID<- reactive({
      req(input$pathiId)
      input$pathiId})


    shinyDirChoose(input, 'folder', root=c(root='~'), filetypes=c('png', 'txt'))

# FUNCTION 
    # Fonction renvoyant le dossier choisi par l'utilisateur
    pathway <- reactive({
      req(input$folder)
      a="~"
      for( i in 1: length(input$folder)){
        a= paste0(a,input$folder$path[i],"/")
      }
     return(a)
     })


    # Fonction renvoyant le nom de l'image correspondant au patway KEGG
    pathView<- reactive({
      req(inputParameter)
      req(inputParameter$orDb)
      dir=pathway()
      kegg_organism = switch(inputParameter$orDb,'org.Dm.eg.db'='dme','org.At.tair.db'='ath',
                           'org.Dr.eg.db'='dre','org.Ss.eg.db'='ssc','org.Cf.eg.db'='cfa',
                           'org.Ag.eg.db'='aga','org.EcSakai.eg.db'='ecs','org.Hs.eg.db'='hsa',
                           'org.Sc.sgd.db'='sce','org.Ce.eg.db'='sko','org.Gg.eg.db'='gga',
                           'org.EcK12.eg.db'='eco','org.Pt.eg.db'='ptr','org.Mxanthus.db'='mxa',
                           'org.Mm.eg.db'='mmu','org.Rn.eg.db'='rno','org.Bt.eg.db'='bta',
                           'org.Mmu.eg.db'='mvv','org.Xl.eg.db'='xla','org.Pf.plasmo.db'='pfa')

      
      setwd(dir)

      pathview(gene.data=kge(), pathway.id=names(choix()[as.numeric(PATHID())]), species = kegg_organism, kegg.dir = dir)

      a=paste0(dir,names(choix()[as.numeric(PATHID())]))
      a=paste0(a,".pathview.png")
     
      return(a)
    })

    choix <- reactive({
      name=gene_list_kegg()$ID
      num  = c(1:length(name))
      choiceTable = data.frame(name, num)
      choix = setNames(as.numeric(choiceTable$num), choiceTable$name)

      return(choix)

    })


    # Fonction renvoyant un objet enrichKEGG
    gene_list_kegg<-eventReactive(input$go,{
      kegg_gene_list=kge()
      kegg_genes <- names(kegg_gene_list)[abs(kegg_gene_list) > input$log2I]

      kegg_organism = switch(inputParameter$orDb,'org.Dm.eg.db'='dme','org.At.tair.db'='ath',
                             'org.Dr.eg.db'='dre','org.Ss.eg.db'='ssc','org.Cf.eg.db'='cfa',
                             'org.Ag.eg.db'='aga','org.EcSakai.eg.db'='ecs','org.Hs.eg.db'='hsa',
                             'org.Sc.sgd.db'='sce','org.Ce.eg.db'='sko','org.Gg.eg.db'='gga',
                             'org.EcK12.eg.db'='eco','org.Pt.eg.db'='ptr','org.Mxanthus.db'='mxa',
                             'org.Mm.eg.db'='mmu','org.Rn.eg.db'='rno','org.Bt.eg.db'='bta',
                             'org.Mmu.eg.db'='mvv','org.Xl.eg.db'='xla','org.Pf.plasmo.db'='pfa')

    
      enrichKEGG(
        gene=kegg_genes,
        universe=names(kegg_gene_list),organism=kegg_organism, pvalueCutoff =input$pvalCutOff, keyType = "ncbi-geneid")

    })


   kge<- eventReactive(input$go,{
       req(inputParameter$fileOr)
      req(inputParameter$orDb)
      BiocManager::install(inputParameter$orDb,update=FALSE)
      library(inputParameter$orDb, character.only = TRUE)

      df=read.csv(inputParameter$fileOr[1,'datapath'], header=TRUE,sep=inputParameter$sep)
      organism=inputParameter$orDb
      if(inputParameter$ora_order=="overexpressed"){
        df=subset(df,log2FoldChange >0)
        #df=subset(df,padj < input$pvalI)
      }
      if(inputParameter$ora_order=="underexpressed"){
        df=subset(df,log2FoldChange <0)
        #df=subset(df, padj < input$pvalI)
      }
      if(inputParameter$ora_order=="both"){
       #df=subset(df, padj < input$pvalI)
     }
  
      df$X=df[,inputParameter$row.names]
      
   
     
      original_gene_list <- df$log2FoldChange
      names(original_gene_list) <- df[,inputParameter$row.names]
      ids<-bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=organism)
      dedup_idsE = ids[!duplicated(ids[c("ENSEMBL")]),]
      dedup_ids = dedup_idsE[!duplicated(dedup_idsE[c("ENTREZID")]),]
      df2 = df[df$X %in% dedup_ids$ENSEMBL,]
      df2$Y = dedup_ids$ENTREZID
      kegg_gene_list <- df2$log2FoldChange
      names(kegg_gene_list) <- df2$Y
      kegg_gene_list<-na.omit(kegg_gene_list)
   
      return(kegg_gene_list)


    })


# OUTPUT
## TABLE
     
    # Fonction renvoyant le tableau de resultat de l'enrichissement
    tabb <- reactive({
      table=cbind(gene_list_kegg()$ID,gene_list_kegg()$Description,gene_list_kegg()$GeneRatio,gene_list_kegg()$BgRatio,gene_list_kegg()$pvalue,gene_list_kegg()$p.adjust,gene_list_kegg()$qvalue,gene_list_kegg()$geneID)
      colnames(table)=cbind('ID','Description','GeneRatio','BgRatio','pvalue','p.adjust','qvalue','geneID')
      return( data.frame(table))
    })

    # Fonction renvoyant à l'ui le datatable es resultats
    output$tableau <- DT::renderDataTable(
      DT::datatable({ tabb()},

        extensions = 'Buttons',
        caption="Table: ORA table of KEGG pathway ",
        options = list(
          dom = 'Bfrtip',
          buttons = c('copy', 'csv', 'excel')
        ),class = "display"
      )
    )
## PLOT

    # Fonction renvoyant l'image KEGG à l'ui
    output$png<- renderPlotly( {
   
      fig=plot_ly()
      fig = fig %>% add_trace(type="image", source= raster2uri(readPNG(pathView())), hoverinfo = "skip")
      return(fig)
    })

 
    # Fonctions renvoyant les plot
    output$dotplotK <- renderPlotly({

      gbb=gene_list_kegg()
      require(DOSE)
      ggplotly(dotplot(gbb, showCategory=input$numberCat))

    })

    output$barplotK <- renderPlotly({

      gbb=gene_list_kegg()


      ggplotly(barplot(gbb,

              showCategory = input$numberCat2,
              title = "GO Biological Pathways",
              font.size = 8))

    })


  }
)}
