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
organism<- reactive({
      organism = inputParameter$orDb
     
      library(organism, character.only = TRUE)

      return(organism)
    })


    shinyDirChoose(input, 'folder2', root=c(root='~'), filetypes=c('png', 'txt'))



    kge<- reactive({
      req(inputParameter$fileOr)
      df=read.csv(inputParameter$fileOr[1,'datapath'], header=TRUE,sep=inputParameter$sep)
      
      if(inputParameter$gsea_order=="log2FC")
      {kegg_gene_list <- df$log2FC}
      if(inputParameter$gsea_order=="pval")
      {kegg_gene_list <- df$padj}
      if(inputParameter$gsea_order=="Stat")
      {kegg_gene_list <- df$stat}
      df$X=df[,inputParameter$row.names]
      names(kegg_gene_list)=df$X

      organism=organism()
      print(organism)
      ids<-bitr(names(kegg_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=organism)
        # remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]

# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 = df[df$X %in% dedup_ids$ENSEMBL,]

# Create a new column in df2 with the corresponding ENTREZ IDs
df2$Y = dedup_ids$ENTREZID

# Create a vector of the gene unuiverse
 if(inputParameter$gsea_order=="log2FC")
      {kegg_gene_list <- df2$log2FC}
      if(inputParameter$gsea_order=="pval")
      {kegg_gene_list <- df2$padj}
      if(inputParameter$gsea_order=="Stat")
      {kegg_gene_list <- df2$stat}


# Name vector with ENTREZ ids
names(kegg_gene_list) <- df2$Y

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

    return(kegg_gene_list)


    })
    # Fonction renvoyant un objet gsepathway
    geneListR<-reactive({
      req(inputParameter$fileOr)
     
Reactome_gene_list = kge()



react_organism = switch(inputParameter$orDb,'org.Dm.eg.db'='fly',
                                'org.Dr.eg.db'='zebrafish','org.Hs.eg.db'='human',
                             'org.Sc.sgd.db'='yeast','org.Ce.eg.db'='celegans',
                             'org.Mm.eg.db'='mouse','org.Rn.eg.db'='rat')


      kk=gsePathway(Reactome_gene_list,
                 pvalueCutoff = 0.2,
                 pAdjustMethod = "BH",
                 verbose = FALSE,organism=react_organism)

      return(kk)

    })


    



    # Fonctions renvoyant les graphiques à l'ui

    output$dotplotR <- renderPlotly({

      gbb=geneListR()
      require(DOSE)
      ggplotly(dotplot(gbb, showCategory=10))

    })

   output$gseaplotKegg <- renderPlotly(ggplotly(gseaplot(kge(), by = "all", title = gene_list_Reactome()$Description[as.numeric(PATHID())], geneSetID =as.numeric( PATHID()))))


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
      table= cbind(geneListR()$ID,geneListR()$Description,geneListR()$enrichmentScore,geneListR()$NES,geneListR()$pvalue,geneListR()$p.adjust,geneListR()$qvalue,geneListR()$rank)
      colnames(table)= c('ID','Description','enrichmentScore','NES','pvalue','p.adjust','qvalues','rank')
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
