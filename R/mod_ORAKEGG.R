#' ORAKEGG UI Function
#'
#' @description A shiny Module.
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


    #Fonction renvoyant l'ui si la méthode choisi est ORA
    output$condPanel<- renderUI({
      LL=c()

      if(inputParameter$Gsea_ORA == 'ORA'){



        LL[[1]]=fluidRow(box(width=13,class='VIOLET',

          h1('ORA KEGG settings',icon('cogs')),
          hr(),
          fluidRow(
            column(width=2),
            column(width=3,

              sliderInput(ns("pvalCutOff"),label = "p-Value cutoff from input ",
                          min = 0, max =0.1, value = 0.05,step=0.01),

            ),
            column(width=3,
                   selectInput(ns("DEG")," Choose element",choices=c('both','over','under')),

            ),
            column(width=3,
                   shinyDirButton(ns('folder2'), 'Select a folder to save pathway plot', 'Please select a folder', FALSE),

            )

          )

        ))
        LL[[2]]= br()

        LL[[3]]=fluidRow(tabsetPanel(
          tabPanel("Table",
                   DT::dataTableOutput(ns('tableau'))%>% withSpinner()
          ),tabPanel("Plot",
                     uiOutput(ns('path_Id')),
                     imageOutput(ns('png'))%>% withSpinner(),
                     plotlyOutput(ns('dotplotK'))%>% withSpinner(),
                     plotlyOutput(ns('barplotK'))%>% withSpinner(),
          ),
        ))



      }
      return(LL)}
    )
    shinyDirChoose(input, 'folder2', root=c(root='~'), filetypes=c('png', 'txt'))

    # Fonction renvoyant un objet enrichKEGG
    geneListK<-reactive({

      req(inputParameter$fileOr)

      df= read.csv(inputParameter$fileOr[1,'datapath'], header=inputParameter$header,sep=inputParameter$sep)
      if(input$DEG=='over'){df=subset(df,log2FoldChange >0)}
      if(input$DEG=='under'){df=subset(df,log2FoldChange <0)}

      organism = inputParameter$orDb
      original_gene_list <- df$log2FoldChange
      names(original_gene_list) <- df$X
      ids<-bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=organism)
      dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]
      df2 = df[df$X %in% dedup_ids$ENSEMBL,]
      df2$Y = dedup_ids$ENTREZID
      kegg_gene_list <- df2$log2FoldChange
      names(kegg_gene_list) <- df2$Y
      kegg_gene_list<-na.omit(kegg_gene_list)
      kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
      kegg_sig_genes_df = subset(df2, padj < 0.05)
      kegg_genes <- kegg_sig_genes_df$log2FoldChange
      names(kegg_genes) <- kegg_sig_genes_df$Y
      kegg_genes <- na.omit(kegg_genes)
      kegg_genes <- names(kegg_genes)[abs(kegg_genes) > 2]

      kegg_organism = switch(inputParameter$orDb,'org.Dm.eg.db'='dme','org.At.tair.db'='ath',
                             'org.Dr.eg.db'='dre','org.Ss.eg.db'='ssc','org.Cf.eg.db'='cfa',
                             'org.Ag.eg.db'='aga','org.EcSakai.eg.db'='ecs','org.Hs.eg.db'='hsa',
                             'org.Sc.sgd.db'='sce','org.Ce.eg.db'='sko','org.Gg.eg.db'='gga',
                             'org.EcK12.eg.db'='eco','org.Pt.eg.db'='ptr','org.Mxanthus.db'='mxa',
                             'org.Mm.eg.db'='mmu','org.Rn.eg.db'='rno','org.Bt.eg.db'='bta',
                             'org.Mmu.eg.db'='mvv','org.Xl.eg.db'='xla','org.Pf.plasmo.db'='pfa')



      kk= enrichKEGG(
        gene=kegg_genes,
        universe=names(kegg_gene_list),organism=kegg_organism, pvalueCutoff =input$pvalCutOff, keyType = "ncbi-geneid")

      return(kk)

    })


    kge<- reactive({
      df= read.csv(inputParameter$fileOr[1,'datapath'], header=inputParameter$header,sep=inputParameter$sep)
      if(input$DEG=='over'){df=subset(df,log2FoldChange >0)}
      if(input$DEG=='under'){df=subset(df,log2FoldChange <0)}


      original_gene_list <- df$log2FoldChange
      names(original_gene_list) <- df$X
      gene_list<-na.omit(original_gene_list)
      gene_list = sort(gene_list, decreasing = TRUE)
    })
    # Fonction renvoyant l'image KEGG à l'ui
    output$png<- renderImage( {list(src = pathView(),
                                    contentType = 'image/png',
                                    width = 400,
                                    height = 300,
                                    alt = "This is alternate text")
    }, deleteFile = FALSE)

    # Fonction renvoyant le dossier choisi par l'utilisateur
    pathway <- reactive({

      a="~"
      for( i in 1: length(input$folder2)){

        a= paste0(a,input$folder2$path[i],"/")}
      return(a)

    })
    # Fonction renvoyant le nom de l'image correspondant au patway KEGG
    pathView<- reactive({


      dir=pathway()
      kegg_organism = switch(inputParameter$orDb,'org.Dm.eg.db'='dme','org.At.tair.db'='ath',
                             'org.Dr.eg.db'='dre','org.Ss.eg.db'='ssc','org.Cf.eg.db'='cfa',
                             'org.Ag.eg.db'='aga','org.EcSakai.eg.db'='ecs','org.Hs.eg.db'='hsa',
                             'org.Sc.sgd.db'='sce','org.Ce.eg.db'='sko','org.Gg.eg.db'='gga',
                             'org.EcK12.eg.db'='eco','org.Pt.eg.db'='ptr','org.Mxanthus.db'='mxa',
                             'org.Mm.eg.db'='mmu','org.Rn.eg.db'='rno','org.Bt.eg.db'='bta',
                             'org.Mmu.eg.db'='mvv','org.Xl.eg.db'='xla','org.Pf.plasmo.db'='pfa')


      d=pathview(gene.data=kge(), pathway.id=names(choix()[as.numeric(PATHID())]), species = kegg_organism, kegg.dir = dir )

      a=paste0(dir,names(choix()[as.numeric(PATHID())]))
      a=paste0(a,".pathview.png")

      return(a)

    })


    # Fonctions renvoyant les plot
    output$dotplotK <- renderPlotly({

      gbb=geneListK()
      require(DOSE)
      ggplotly(dotplot(gbb, showCategory=10))

    })

    output$barplotK <- renderPlotly({

      gbb=geneListK()


      ggplotly(barplot(gbb,

              showCategory = 10,
              title = "GO Biological Pathways",
              font.size = 8))

    })







    # Fonction renvoyant à l'ui les différent pathway significativement enrichi sous la forme d'un select input
    output$path_Id <- renderUI({

      L= selectInput(ns('ID_path'),"Select Path Id to observe" ,choices = choix())
      return(L)

    })
    PATHID<- reactive(input$ID_path)

    # Fonction renvoyant les différent pathway significativement enrichi
  choix <- reactive({

    name=geneListK()$ID

    num  = c(1:length(name))


    choiceTable = data.frame(name, num)
    choix = setNames(as.numeric(choiceTable$num), choiceTable$name)
    return(choix)

  })
  # Fonction renvoyant l'image KEGG à l'ui
  output$png<- renderImage( {list(src = pathView(),
                                  contentType = 'image/png',
                                  width = 400,
                                  height = 300,
                                  alt = "This is alternate text")
  }, deleteFile = FALSE)

  # Fonction renvoyant le tableau de resultat de l'enrichissement
  tabb <- reactive({


    table=cbind(geneListK()$ID,geneListK()$Description,geneListK()$GeneRatio,geneListK()$BgRatio,geneListK()$pvalue,geneListK()$p.adjust,geneListK()$qvalue,geneListK()$geneID)
    colnames(table)=cbind('ID','Description','GeneRatio','BgRatio','pvalue','p.adjust','qvalue','geneID')
    return( data.frame(table))
  })

  # Fonction renvoyant à l'ui le datatable es resultats
  output$tableau <- DT::renderDataTable(DT::datatable(
    { tabb()},

    extensions = 'Buttons',
  caption="Table: ORA table of KEGG pathway ",

    options = list(

      dom = 'Bfrtip',
      buttons = c('copy', 'csv', 'excel')
    ),

    class = "display"
  ))



  }
)}

## To be copied in the UI
# mod_ORAKEGG_ui("ORAKEGG_1")

## To be copied in the server
# mod_ORAKEGG_server("ORAKEGG_1")
