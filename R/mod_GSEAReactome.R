#' GSEAReactome UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_GSEAReactome_ui <- function(id){
  ns <- NS(id)
  tagList(
    uiOutput(ns('condPanel'))
  )
}

#' GSEAReactome Server Functions
#'
#' @noRd
mod_GSEAReactome_server <- function(id,inputParameter){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
    library(ReactomePA)
    ns <- session$ns

    #Fonction renvoyant l'ui si la méthode choisi est GSEA
    output$condPanel<- renderUI({
      LL = c()

      if(inputParameter$Gsea_ORA == 'GSEA'){


        LL[[1]] = fluidRow(box(width=13, class='VIOLET',

          h1('GSEA Reactome settings',icon('cogs')),
          hr(),
          fluidRow(

            column(width=4),
            column(width=5,
              sliderInput(ns("pvalCutOff"),label = "p-Value cutoff from input ",
                      min = 0, max =0.1, value = 0.05,step=0.01)
            )
          )
        ))
        LL[[2]] = br()

        LL[[3]] = fluidRow( tabsetPanel(

          tabPanel("Table",
            DT::dataTableOutput(ns('tableau'))%>% withSpinner()
          ),
          tabPanel("Plot",
            plotlyOutput(ns('dotplotR'))%>% withSpinner(),
            plotlyOutput(ns('gseaplotR'))%>% withSpinner(),
            plotlyOutput(ns('pathwayR'))%>% withSpinner()
          ),
        ))

      }
      return(LL)
    })


    shinyDirChoose(input, 'folder2', root=c(root='~'), filetypes=c('png', 'txt'))

    # Fonction renvoyant un objet gsepathway
    geneListR<-reactive({
      req(inputParameter$fileOr)
      df = read.csv(inputParameter$fileOr[1,'datapath'], header=inputParameter$header, sep=inputParameter$sep)

      organism = inputParameter$orDb
      original_gene_list <- df$log2FoldChange
      names(original_gene_list) <- df$X
      ids<-bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=organism)
      dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]
      df2 = df[df$X %in% dedup_ids$ENSEMBL,]
      df2$Y = dedup_ids$ENTREZID
      react_gene_list <- df2$log2FoldChange

      names(react_gene_list) <- df2$Y
      react_gene_list<-na.omit(react_gene_list)
      react_gene_list = sort(react_gene_list, decreasing = TRUE)
      react_sig_genes_df = subset(df2, padj < 0.05)
      react_genes <- react_sig_genes_df$log2FoldChange
      names(react_genes) <- react_sig_genes_df$Y
      react_genes <- na.omit(react_genes)
      react_genes <- names(react_genes)[abs(react_genes) > 2]



      kk=gsePathway(react_gene_list,
                 pvalueCutoff = 0.2,
                 pAdjustMethod = "BH",
                 verbose = FALSE)

      return(kk)

    })


    kge<- reactive({
      df= read.csv(inputParameter$fileOr[1,'datapath'], header=inputParameter$header,sep=inputParameter$sep)

      organism = inputParameter$orDb
      original_gene_list <- df$log2FoldChange

      names(original_gene_list) <- df$X

      ids<-bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=organism)

      dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]


      df2 = df[df$X %in% dedup_ids$ENSEMBL,]

      df2$Y = dedup_ids$ENTREZID

      react_gene_list <- df2$log2FoldChange

      names(react_gene_list) <- df2$Y

      react_gene_list<-na.omit(react_gene_list)

      react_gene_list = sort(react_gene_list, decreasing = TRUE)
      duplicated_names <- duplicated(names(react_gene_list))

      react_gene_list=react_gene_list[!duplicated_names]

      return(react_gene_list)


    })



    # Fonctions renvoyant les graphiques à l'ui

    output$dotplotR <- renderPlotly({

      gbb=geneListR()
      require(DOSE)
      ggplotly(dotplot(gbb, showCategory=10))

    })

   output$gseaplotKegg <- renderPlotly(ggplotly(gseaplot(kge(), by = "all", title = gene_list_kegg()$Description[as.numeric(PATHID())], geneSetID =as.numeric( PATHID()))))


   # Fonction renvoyant le nom de l'image correspondant au pathway Reactome
    output$pathwayR <- renderPlotly({

      gbb=geneListR()


      kg=switch(inputParameter$orDb,'org.Dm.eg.db'='dme','org.At.tair.db'='ath',
                'org.Dr.eg.db'='dre','org.Ss.eg.db'='ssc','org.Cf.eg.db'='cfa',
                'org.Ag.eg.db'='aga','org.EcSakai.eg.db'='ecs','org.Hs.eg.db'='hsa',
                'org.Sc.sgd.db'='sce','org.Ce.eg.db'='sko','org.Gg.eg.db'='gga',
                'org.EcK12.eg.db'='eco','org.Pt.eg.db'='ptr','org.Mxanthus.db'='mxa',
                'org.Mm.eg.db'='mmu','org.Rn.eg.db'='rno','org.Bt.eg.db'='bta',
                'org.Mmu.eg.db'='mvv','org.Xl.eg.db'='xla','org.Pf.plasmo.db'='pfa')
      react=switch(kg,'dme'='fly' ,'ath'='arabidopsis','dre'='zebrafish','ssc'='pig','cfa'='canine','ecs'='ecolik12','hsa'='human','sce'='yeast','sko'='celegans','gga'='chicken','eco'='ecolik12','mmu'='mouse','rno'='rat','bta'='bovine','xla'='xenopus')


      ggplotly(viewPathway(gbb$Description[1],
                           readable = TRUE,
                           organism=react,
                           foldChange = kge()))

    })

    # Fonction renvoyant le tableau de resultat de l'enrichissement
    tabb <- reactive({
      table= cbind(geneListR()$ID,geneListR()$Description,geneListR()$enrichmentScore,geneListR()$NES,geneListR()$pvalue,geneListR()$p.adjust,geneListR()$qvalue,geneListR()$rank,geneListR()$leading_edge,geneListR()$core_enrichment)
      colnames(table)= c('ID','Description','enrichmentScore','NES','pvalue','p.adjust','qvalues','rank','leadingEdge','coreEnrichment')
      return(data.frame(table))

    })

    # Fonction renvoyant à l'ui la datatable des resultats
    output$tableau <- DT::renderDataTable(DT::datatable(
      { tabb()},

      extensions = 'Buttons',
      caption="Table: GSEA table of Reactome Pathway",


      options = list(

        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel')
      ),

      class = "display"
    ))













  })
}

## To be copied in the UI
# mod_GSEAReactome_ui("GSEAReactome_1")

## To be copied in the server
# mod_GSEAReactome_server("GSEAReactome_1")
