#' ORAReactome UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_ORAReactome_ui <- function(id){
  ns <- NS(id)
  tagList(
    uiOutput(ns('condPanel'))
  )
}

#' ORAReactome Server Functions
#'
#' @noRd
mod_ORAReactome_server <- function(id,inputParameter){
  moduleServer( id, function(input, output, session){
    library(ReactomePA)
    ns <- session$ns


    #Fonction renvoyant l'ui si la méthode choisi est ORA
    output$condPanel<- renderUI({
      LL=c()

      if(inputParameter$Gsea_ORA == 'ORA'){

        LL[[1]] = fluidRow(box(width=13,class='VIOLET',

          h1('ORA Reactome settings',icon('cogs')),
          hr(),
          fluidRow(
            column(width=2),
            column(width=3,

              sliderInput(ns("pvalCutOff"),label = "p-Value cutoff from input ",
                          min = 0, max =0.1, value = 0.05,step=0.01),

            ),
            column(width=3,sliderInput(ns("log2FC"),label="LOg2FC fold change ",min=0,max=5,value=1,step=0.1)
                   

            ),
            column(width=3
                   
            )
          )

        ))

        LL[[2]] = br()
        LL[[3]] = fluidRow(tabsetPanel(
          tabPanel("Table",
                   DT::dataTableOutput(ns('tableau'))%>% withSpinner()
          ),tabPanel("Plot",


                     plotlyOutput(ns('dotplotR'))%>% withSpinner(),
                     plotlyOutput(ns('barplotR'))%>% withSpinner(),
                     plotlyOutput(ns('pathwayR'))%>% withSpinner()
          ),
        ))

      }
      return(LL)}
    )


    shinyDirChoose(input, 'folder2', root=c(root='~'), filetypes=c('png', 'txt'))


    # Fonction renvoyant un objet enrichpathway
    geneListR<-reactive({
  req(inputParameter$fileOr)

      df=read.csv(inputParameter$fileOr[1,'datapath'], header=TRUE,sep=inputParameter$sep)
      organism=inputParameter$orDb
      if(inputParameter$ora_order=="overexpressed"){print('over')
      df=subset(df,log2FC >0)
      df=subset(df,padj< input$pvalCutOff)}
      if(inputParameter$ora_order=="underexpressed"){print('under')
      print(head(df))
      df=subset(df,log2FC <0)
      print(head(df))
      df=subset(df, padj < input$pvalCutOff)
       print(head(df))}
      if(inputParameter$ora_order=="both"){ df=subset(df, padj < input$pvalCutOff)}
  
      df$X=df[,inputParameter$row.names]
 
      library(organism, character.only = TRUE)
      original_gene_list <- df$log2FC



      names(original_gene_list) <- df[,inputParameter$row.names]

      ids<-bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=organism)
        # remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
        dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]

      # Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
      df2 = df[df$X %in% dedup_ids$ENSEMBL,]

      # Create a new column in df2 with the corresponding ENTREZ IDs
      df2$Y = dedup_ids$ENTREZID

      kegg_gene_list <- df2$log2FC

      names(kegg_gene_list) <- df2$Y

      kegg_gene_list<-na.omit(kegg_gene_list)
      kegg_genes <- names(kegg_gene_list)[abs(kegg_gene_list) > input$log2FC]

react_organism = switch(inputParameter$orDb,'org.Dm.eg.db'='fly',
                                'org.Dr.eg.db'='zebrafish','org.Hs.eg.db'='human',
                             'org.Sc.sgd.db'='yeast','org.Ce.eg.db'='celegans',
                             'org.Mm.eg.db'='mouse','org.Rn.eg.db'='rat')

    
      kk= enrichPathway(
        gene=kegg_genes,pvalueCutoff = input$pvalCutOff, readable=TRUE,organism=react_organism

       )

      return(kk)

    })


   kge<- reactive({
       req(inputParameter$fileOr)

      df=read.csv(inputParameter$fileOr[1,'datapath'], header=TRUE,sep=inputParameter$sep)
      organism=inputParameter$orDb
      if(inputParameter$ora_order=="overexpressed"){print('over')
      df=subset(df,log2FC >0)
      df=subset(df,padj< input$pvalCutOff)}
      if(inputParameter$ora_order=="underexpressed"){print('under')
      print(head(df))
      df=subset(df,log2FC <0)
      print(head(df))
      df=subset(df, padj < input$pvalCutOff)
       print(head(df))}
      if(inputParameter$ora_order=="both"){ df=subset(df, padj < input$pvalCutOff)}
  
      df$X=df[,inputParameter$row.names]
      
   
      library(organism, character.only = TRUE)
      original_gene_list <- df$log2FC



      names(original_gene_list) <- df[,inputParameter$row.names]

      ids<-bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=organism)
        # remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
      dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]

      # Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
      df2 = df[df$X %in% dedup_ids$ENSEMBL,]

      # Create a new column in df2 with the corresponding ENTREZ IDs
      df2$Y = dedup_ids$ENTREZID

      kegg_gene_list <- df2$log2FC

      names(kegg_gene_list) <- df2$Y

      kegg_gene_list<-na.omit(kegg_gene_list)

      

      
  


      kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
      duplicated_names <- duplicated(names(kegg_gene_list))


      kegg_gene_list=kegg_gene_list[!duplicated_names]

      return(kegg_gene_list)


    })

# Fonctions renvoyant les graphiques à l'ui
    output$dotplotR <- renderPlotly({

      gbb=geneListR()
      require(DOSE)
      ggplotly(dotplot(gbb, showCategory=10))

    })

    output$barplotR <- renderPlotly({

      gbb=geneListR()


      ggplotly(barplot(gbb,

                       showCategory = 10,
                       title = "Barplot Recatome Pathways",
                       font.size = 8))

    })


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


      table=cbind(geneListR()$ID,geneListR()$Description,geneListR()$GeneRatio,geneListR()$BgRatio,geneListR()$pvalue,geneListR()$p.adjust,geneListR()$qvalue,geneListR()$geneID)
      colnames(table)=cbind('ID','Description','GeneRatio','BgRatio','pvalue','p.adjust','qvalue','geneID')
      return( data.frame(table))
    })

    # Fonction renvoyant à l'ui la datatable des resultats
    output$tableau <- DT::renderDataTable(DT::datatable(
      { tabb()},

      extensions = 'Buttons',
      caption="Table: ORA table of Reactome Pathway",
      options = list(

        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel')
      ),

      class = "display"
    ))



  }
  )
}

## To be copied in the UI
# mod_ORAReactome_ui("ORAReactome_1")

## To be copied in the server
# mod_ORAReactome_server("ORAReactome_1")
