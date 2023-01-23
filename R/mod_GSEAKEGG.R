#' GSEAKEGG UI Function
#'
#' @description A shiny Module.
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

    ns <- session$ns


    #Fonction renvoyant l'ui si la méthode choisi est GSEA
    output$condPanel<- renderUI({
      LL=c()

      if(inputParameter$Gsea_ORA == 'GSEA'){



        LL[[1]] = fluidRow(box(width=13,class='VIOLET',

          h1('GSEA KEGG settings',icon('cogs')),
          hr(),
          fluidRow(
            column(width=3.5),
            column(width=3,
                   sliderInput(ns("pvalCutOff"),label = "p-Value cutoff from output ",
                               min = 0, max =0.1, value = 0.05,step=0.01),
            ),
            column(width=3,
                  shinyDirButton(ns('folder'), 'Select a folder to save pathway plot', 'Please select a folder', FALSE),
            )
          )
        ))

        LL[[2]] =  br()
        LL[[3]] = fluidRow(tabsetPanel(
          tabPanel("Table",
                   DT::dataTableOutput(ns('tableau'))%>% withSpinner()
          ),tabPanel("Plot",
                     uiOutput(ns('pathIdd')),
                     imageOutput(ns('png'))%>% withSpinner(),
                     plotlyOutput(ns('dotplotKegg'))%>% withSpinner(),
                     plotlyOutput(ns('gseaplotKegg'))%>% withSpinner(),
          ),
        ))



      }
      return(LL)}
    )

    # Fonction chargant la librairie correspondante à l'organisme choisi
    organism<- reactive({
      organism = inputParameter$orDb
      BiocManager::install(organism, character.only = TRUE)
      library(organism, character.only = TRUE)

      return(organism)
    })


    # Fonction renvoyant un objet gseKEGG
    gene_list_kegg <- reactive({
      req(inputParameter$fileOr)
      organism = organism()
      df= read.csv(inputParameter$fileOr[1,'datapath'], header=inputParameter$header,sep=inputParameter$sep)

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

      kegg_organism = switch(inputParameter$orDb,'org.Dm.eg.db'='dme','org.At.tair.db'='ath',
                                'org.Dr.eg.db'='dre','org.Ss.eg.db'='ssc','org.Cf.eg.db'='cfa',
                             'org.Ag.eg.db'='aga','org.EcSakai.eg.db'='ecs','org.Hs.eg.db'='hsa',
                             'org.Sc.sgd.db'='sce','org.Ce.eg.db'='sko','org.Gg.eg.db'='gga',
                             'org.EcK12.eg.db'='eco','org.Pt.eg.db'='ptr','org.Mxanthus.db'='mxa',
                             'org.Mm.eg.db'='mmu','org.Rn.eg.db'='rno','org.Bt.eg.db'='bta',
                             'org.Mmu.eg.db'='mvv','org.Xl.eg.db'='xla','org.Pf.plasmo.db'='pfa')


       kk2 <- gseKEGG(geneList     = kegg_gene_list,
                     organism     = kegg_organism,
                     nPerm        = 10000,
                     minGSSize    = 3,
                     maxGSSize    = 800,
                     pvalueCutoff = input$pvalCutOff,
                     pAdjustMethod = inputParameter$pAdjMethod,
                     keyType       = "ncbi-geneid")


      return(kk2)




    })
    # Fonction renvoyant les différent pathway significativement enrichi
    choix <- reactive({
      name=gene_list_kegg()$ID
      num  = c(1:length(name))
      choiceTable = data.frame(name, num)
      choix = setNames(as.numeric(choiceTable$num), choiceTable$name)

      return(choix)

    })

    # Fonction renvoyant à l'ui les différent pathway significativement enrichi sous la forme d'un select input

  output$pathIdd <- renderUI({

    L= selectInput(ns('pathiId'),"Select Path Id to observe" ,choices = choix())
    return(L)

  })


  shinyDirChoose(input, 'folder', root=c(root='~'), filetypes=c('png', 'txt'))


  # Fonction renvoyant la liste des gènes triés par log2FC decroissant
  kge <- reactive({
    organism = organism()
    df= read.csv(inputParameter$fileOr[1,'datapath'], header=inputParameter$header,sep=inputParameter$sep)

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

    return(kegg_gene_list)})

    PATHID<- reactive(input$pathiId)



    # Fonction renvoyant les plot à l'ui
    output$dotplotKegg <- renderPlotly(ggplotly( dotplot(gene_list_kegg(), showCategory = 10, title = "Enriched Pathways" ) ))
    output$gseaplotKegg <- renderPlotly(ggplotly(gseaplot(gene_list_kegg(), by = "all", title = gene_list_kegg()$Description[as.numeric(PATHID())], geneSetID =as.numeric( PATHID()))))


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
     for( i in 1: length(input$folder)){

       a= paste0(a,input$folder$path[i],"/")}
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

    setwd(dir)
    pathview(gene.data=kge(), pathway.id=names(choix()[as.numeric(PATHID())]), species = kegg_organism, kegg.dir = dir,limit = list(gene=max(abs(kge())), cpd=1 ))

    a=paste0(dir,names(choix()[as.numeric(PATHID())]))
    print(a)
    a=paste0(a,".pathview.png")
    print(a)


    return(a)

   })

  # Fonction renvoyant le tableau de resultat de l'enrichissement
  tabb <- reactive({
    table= cbind(gene_list_kegg()$ID,gene_list_kegg()$Description,gene_list_kegg()$enrichmentScore,gene_list_kegg()$NES,gene_list_kegg()$pvalue,gene_list_kegg()$p.adjust,gene_list_kegg()$qvalue,gene_list_kegg()$rank,gene_list_kegg()$leading_edge,gene_list_kegg()$core_enrichment)
    colnames(table)= c('ID','Description','enrichmentScore','NES','pvalue','p.adjust','qvalues','rank','leadingEdge','coreEnrichment')
    return(data.frame(table))

  })


  # Fonction renvoyant à l'ui le datatable es resultats
  output$tableau <- DT::renderDataTable(DT::datatable(
  { tabb()},

    extensions = 'Buttons',
    caption="Table: GSEA table of KEGG Pathway",
    options = list(

      dom = 'Bfrtip',
      buttons = c('copy', 'csv', 'excel')
    ),

    class = "display"
  ))








  })
}

## To be copied in the UI
# mod_GSEAKEGG_ui("GSEAKEGG_1")

## To be copied in the server
# mod_GSEAKEGG_server("GSEAKEGG_1")
