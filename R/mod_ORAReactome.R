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


# UI OUTPUT
    #Fonction renvoyant l'ui si la méthode choisi est ORA
    output$condPanel<- renderUI({
      LL=c()

      if(inputParameter$Gsea_ORA == 'ORA'){
        LL[[1]] = fluidRow(box(width = 12,status = 'info',solidHeader = TRUE,collapsible=TRUE,title=h1('ORA Reactome Settings',icon('gear')),
          fluidRow(
            column(width=1),
            column(width=3,
              sliderInput(ns("pvalI"),label = "p-Value cutoff from input ",
                          min = 0, max =0.1, value = 0.05,step=0.01),
            ),
            column(width=3,sliderInput(ns("log2I"),label = "Log2 Fold change cutoff from input ",
                          min = 0, max =5, value = 2,step=0.5)
            ),
            column(width=3, sliderInput(ns("pvalCutOff"),label = "p-Value cutoff from output ",
                          min = 0, max =0.1, value = 0.05,step=0.01)
            ),
            column(width=2,br(),br(),br(),actionButton(ns("go"),"Start",class="buttS",icon("play")))
          )
        ))

        LL[[2]] = fluidRow( tabBox(width=12,height=NULL,
          tabPanel("Table",fluidRow(
            column(width=12,
              box(width = NULL,status = 'primary',solidHeader = TRUE,title=h1('Table of results',icon('table')),
                DT::dataTableOutput(ns('tableau'))%>% withSpinner(color="#605CA8")
              )
            )
          )),
          tabPanel("Plot",
            fluidRow(
            column(width=12,
            box(width = NULL,status = 'success',solidHeader = TRUE,title=h1('Reactome pathway',icon('chart-simple')),
                  uiOutput(ns('pathIdd')),
                  uiOutput(ns('titlePath')),
                 plotOutput(ns('pathwayR'))%>% withSpinner(color="#CDCDE6")
            )),
            column(width=6,
              box(width = NULL,status = 'success',solidHeader = TRUE,title=h1('Dotplot',icon('chart-simple')),
                uiOutput(ns('numberCategrory')),
                plotlyOutput(ns('dotplotR'))%>% withSpinner(color="#CDCDE6")
              )
            ),
            column(width=6,
              box(width = NULL,status = 'success',solidHeader = TRUE,title=h1('Barplot',icon('chart-simple')),
                uiOutput(ns('numberCategrory2')),
                plotlyOutput(ns('barplotR'))%>% withSpinner(color="#CDCDE6")
              )
            
            
            )
          ))
        ))
      }
      return(LL)
    })


    shinyDirChoose(input, 'folder2', root=c(root='~'), filetypes=c('png', 'txt'))

    # Fonction renvoyant à l'ui les différent pathway significativement enrichi sous la forme d'un select input pour le pathway KEGG
    output$pathIdd <- renderUI({

      selectInput(ns('pathiId'),"Select Path Id to observe" ,choices = choix(),width="20%")
    })
    PATHID<- reactive({
      req(input$pathiId)
      input$pathiId})

  # Fonction renvoyant le nom du pathway reactome observe   
  output$titlePath <- renderUI({
      num_path=as.numeric( PATHID())
      h4(geneListR()$Description[num_path])
      })
       # Fonction renvoyant à l'ui un slideerInput pour deteminer le nombre de categrories a afficher sur le dotPlot
    output$numberCategrory <- renderUI({
       numericInput(ns("numberCat"), min=1,value=min(5,nrow(tabb())),max=nrow(tabb()),label="Write number of category to show",width="20%")

    })
        output$numberCategrory2 <- renderUI({
       numericInput(ns("numberCat2"), min=1,value=min(5,nrow(tabb())),max=nrow(tabb()),label="Write number of category to show",width="20%")

    })

# FUNCTION 
    # Fonction renvoyant un objet enrichpathway
    geneListR<-eventReactive(input$go,{
      
      kegg_gene_list=kge()
      kegg_genes = names(kegg_gene_list)[abs(kegg_gene_list) > input$log2I]

      react_organism = switch(inputParameter$orDb,'org.Dm.eg.db'='fly',
                             'org.Dr.eg.db'='zebrafish','org.Hs.eg.db'='human',
                             'org.Sc.sgd.db'='yeast','org.Ce.eg.db'='celegans',
                             'org.Mm.eg.db'='mouse','org.Rn.eg.db'='rat')

      oraPath= enrichPathway(
        gene=kegg_genes,pvalueCutoff = input$pvalCutOff, readable=TRUE,organism=react_organism
      )
      return(oraPath)

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
        df=subset(df,padj< input$pvalCutOff)
      }
      if(inputParameter$ora_order=="underexpressed"){
        df=subset(df,log2FoldChange <0)
        df=subset(df, padj < inputParameter$pvalI)
      }
      if(inputParameter$ora_order=="both"){ df=subset(df, padj < input$pvalCutOff)}
  
      df$X=df[,inputParameter$row.names]
 
    
      original_gene_list = df$log2FoldChange

      names(original_gene_list) = df[,inputParameter$row.names]

      ids=bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=organism)
      dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]

      df2 = df[df$X %in% dedup_ids$ENSEMBL,]
      df2$Y = dedup_ids$ENTREZID

      kegg_gene_list = df2$log2FoldChange

      names(kegg_gene_list) = df2$Y

      kegg_gene_list=na.omit(kegg_gene_list)

      return(kegg_gene_list)


    })
   # Fonction renvoyant les différent GO terms significativement enrichi
    choix <- reactive({

      name=geneListR()$ID
      num  = c(1:length(name))
      choiceTable = data.frame(name, num)
      choix = setNames(as.numeric(choiceTable$num), choiceTable$name)

      return(choix)

    })

# OUTPUT
## TABLE
  # Fonction renvoyant le tableau de resultat de l'enrichissement
    tabb <- reactive({
      table=cbind(geneListR()$ID,geneListR()$Description,geneListR()$GeneRatio,geneListR()$BgRatio,geneListR()$pvalue,geneListR()$p.adjust,geneListR()$qvalue,geneListR()$geneID)
      colnames(table)=cbind('ID','Description','GeneRatio','BgRatio','pvalue','p.adjust','qvalue','geneID')
      return( data.frame(table))
    })

    # Fonction renvoyant à l'ui la datatable des resultats
    output$tableau <- DT::renderDataTable(
      DT::datatable({
       tabb()
      },
      extensions = 'Buttons',
      caption="Table: ORA table of Reactome Pathway",
      options = list(
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel')
      ), class = "display"
    ))
## PLOT  
  # Fonctions renvoyant le dot plot
    output$dotplotR <- renderPlotly({
      req(input$numberCat)
      require(DOSE)
      ggplotly(dotplot(geneListR(), showCategory=input$numberCat))

    })
    # Fonctions renvoyant le barplot
    output$barplotR <- renderPlotly({
      ggplotly(barplot(geneListR(),

                       showCategory = input$numberCat2,
                       title = "Barplot Recatome Pathways",
                       font.size = 8))

    })


    # Fonction renvoyant le nom de l'image correspondant au pathway Reactome
    output$pathwayR <- renderPlot({
      req(inputParameter)
      kg=switch(inputParameter$orDb,'org.Dm.eg.db'='dme','org.At.tair.db'='ath',
                'org.Dr.eg.db'='dre','org.Ss.eg.db'='ssc','org.Cf.eg.db'='cfa',
                'org.Ag.eg.db'='aga','org.EcSakai.eg.db'='ecs','org.Hs.eg.db'='hsa',
                'org.Sc.sgd.db'='sce','org.Ce.eg.db'='sko','org.Gg.eg.db'='gga',
                'org.EcK12.eg.db'='eco','org.Pt.eg.db'='ptr','org.Mxanthus.db'='mxa',
                'org.Mm.eg.db'='mmu','org.Rn.eg.db'='rno','org.Bt.eg.db'='bta',
                'org.Mmu.eg.db'='mvv','org.Xl.eg.db'='xla','org.Pf.plasmo.db'='pfa')
      react=switch(kg,'dme'='fly' ,'ath'='arabidopsis','dre'='zebrafish','ssc'='pig','cfa'='canine','ecs'='ecolik12','hsa'='human','sce'='yeast','sko'='celegans','gga'='chicken','eco'='ecolik12','mmu'='mouse','rno'='rat','bta'='bovine','xla'='xenopus')
      num_path=as.numeric( PATHID())
      
     viewPathway(geneListR()$Description[num_path+1],
                           readable = TRUE,
                           organism=react,
                           foldChange = kge())
      
      

    })    
  })
}
