#' GSEAReactome UI Function
#'
#' @description This module runs Over Representation Analysis on Reactome pathway.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @author Audrey BEAUFILS
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

    library(ReactomePA)

    ns <- session$ns

# UI OUTPUT
    #Fonction renvoyant l'ui si la méthode choisi est GSEA
    output$condPanel<- renderUI({
      LL = c()

      if(inputParameter$Gsea_ORA == 'GSEA'){


        LL[[1]] = fluidRow(box(width = 12,status = 'info',solidHeader = TRUE,collapsible=TRUE,title=h1('GSEA Reactome Settings',icon('gear')),
          fluidRow(
            column(width=1),
            column(width=3,
              sliderInput(ns("pvalCutOff"),label = "p-Value cutoff from output ",
                      min = 0, max =0.1, value = 0.05,step=0.01)
            ),
            column(width=2,
              br(),br(),br(),
              actionButton(ns("go"),"Start",icon("play"))
            )
          )
        ))
        LL[[2]] = fluidRow( tabBox(width=12,height=NULL,
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
                box(width = NULL,status = 'success',solidHeader = TRUE,title=h1('Reactome pathway',icon('chart-simple')),
                  uiOutput(ns('pathIdd2')),
                  uiOutput(ns('titlePath')),
                  plotOutput(ns('pathwayR'))%>% withSpinner(color="#CDCDE6")
                )
              ),
              column(width=6,
                box(width = NULL,status = 'success',solidHeader = TRUE,title=h1('Dotplot',icon('chart-simple')),
                  uiOutput(ns('numberCategory')),
                  plotlyOutput(ns('dotplotR'))%>% withSpinner(color="#CDCDE6")
                )
              ),
              column(width=6,
                box(width = NULL,status = 'success',solidHeader = TRUE,title=h1('GSEAplot',icon('chart-simple')),
                  uiOutput(ns('pathIdd')),
                  plotOutput(ns('gseaplotR'))%>% withSpinner(color="#CDCDE6")
                )
              )
            )
          )
        ))
      }
      return(LL)
    })

    shinyDirChoose(input, 'folder2', root=c(root='~'), filetypes=c('png', 'txt'))

    output$numberCategory <- renderUI({
        numericInput(ns("numberCat"), min=1,value=min(5,nrow(tabb())),max=nrow(tabb()),label="Write number of category to show",width="20%")

    })
    # Fonction renvoyant à l'ui les différent pathway significativement enrichi sous la forme d'un select input pour le pathway KEGG
    output$pathIdd <- renderUI({

      selectInput(ns('pathiId'),"Select Path Id to observe" ,
        choices = choix(),width="20%")
    })
    PATHID<- reactive({
      req(input$pathiId)
      input$pathiId})

    # Fonction renvoyant le nom du pathway reactome observe   
    output$titlePath <- renderUI({
      num_path=as.numeric( PATHID2())
      h4(geneListR()$Description[num_path])
    })
    
    # Fonction renvoyant à l'ui les différent pathway significativement enrichi sous la forme d'un select input pour le pathway KEGG
    output$pathIdd2 <- renderUI({

      selectInput(ns('pathiId2'),"Select Path Id to observe" ,choices = choix(),width="20%")
    })
    PATHID2<- reactive(input$pathiId2)



#FUNCTION
    # Fonction renvoyant la liste des gènes triés par log2FoldChange decroissant
    kge <- eventReactive(input$go,{
      req(inputParameter$fileOr)
      req(inputParameter$orDb)
      BiocManager::install(inputParameter$orDb,update=FALSE)
      library(inputParameter$orDb, character.only = TRUE)

      df = read.csv(inputParameter$fileOr[1,'datapath'], header=TRUE,sep=inputParameter$sep)
      
      if(inputParameter$gsea_order=="log2FoldChange") {kegg_gene_list = df$log2FoldChange}
      if(inputParameter$gsea_order=="pval")   {kegg_gene_list = df$padj}
      if(inputParameter$gsea_order=="Stat")   {kegg_gene_list = df$stat}
      
      df$X=df[,inputParameter$row.names]
      names(kegg_gene_list)=df$X


      ids=bitr(names(kegg_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=inputParameter$orDb)
      dedup_idsE = ids[!duplicated(ids[c("ENSEMBL")]),]
      dedup_ids = dedup_idsE[!duplicated(dedup_idsE[c("ENTREZID")]),]
      df2 = df[df$X %in% dedup_ids$ENSEMBL,]

      df2$Y = dedup_ids$ENTREZID
      

      if(inputParameter$gsea_order=="log2FoldChange")  {kegg_gene_list = df2$log2FoldChange}
      if(inputParameter$gsea_order=="pval")    {kegg_gene_list = df2$padj}
      if(inputParameter$gsea_order=="Stat")    {kegg_gene_list = df2$stat}

      names(kegg_gene_list) = df2$Y
      kegg_gene_list=na.omit(kegg_gene_list)
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


      gseaPath=gsePathway(Reactome_gene_list,
                 pvalueCutoff = input$pvalCutOff,
                 pAdjustMethod = inputParameter$pAdjMethod,
                 verbose = FALSE,organism=react_organism)

      return(gseaPath)

    })
       # Fonction renvoyant les différent GO terms significativement enrichi
    choix <- reactive({

      name=geneListR()$ID
      num  = c(1:length(name))
      choiceTable = data.frame(name, num)
      choix = setNames(as.numeric(choiceTable$num), choiceTable$name)

      return(choix)

    })

#OUTPUT
## TABLE

    # Fonction renvoyant le tableau de resultat de l'enrichissement
    tabb <- reactive({

      table = cbind(geneListR()$ID,geneListR()$Description,geneListR()$enrichmentScore,geneListR()$NES,geneListR()$pvalue,geneListR()$p.adjust,geneListR()$qvalue)
      colnames(table)= c('ID','Description','enrichmentScore','NES','pvalue','p.adjust','qvalues')
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




##PLOT 
 
    output$dotplotR <- renderPlotly({
      req(input$numberCat)
      gseaPath=geneListR()
      require(DOSE)
      ggplotly(dotplot(gseaPath, showCategory=input$numberCat))

    })

    output$gseaplotR<- renderPlot(
      enrichplot::gseaplot2(geneListR(),
          title = geneListR()$Description[as.numeric(PATHID())],
          geneSetID =as.numeric( PATHID())
      )
    )
    


   # Fonction renvoyant le nom de l'image correspondant au pathway Reactome
    output$pathwayR <- renderPlot({
      req(inputParameter)
      gseaPath=geneListR()
      print('aaaaa')
      print(gseaPath)

      organismR=switch(inputParameter$orDb,'org.Dm.eg.db'='fly','org.At.tair.db'='arabidopsis',
                'org.Dr.eg.db'='zebrafish','org.Ss.eg.db'='pig','org.Cf.eg.db'='canine',
                'org.Ag.eg.db'=NA,'org.EcSakai.eg.db'='ecolik12','org.Hs.eg.db'='human',
                'org.Sc.sgd.db'='yeast','org.Ce.eg.db'='celegans','org.Gg.eg.db'='chicken',
                'org.EcK12.eg.db'='ecolik12','org.Pt.eg.db'=NA,'org.Mxanthus.db'=NA,
                'org.Mm.eg.db'='mouse','org.Rn.eg.db'='rat','org.Bt.eg.db'='bovine',
                'org.Mmu.eg.db'=NA,'org.Xl.eg.db'='xenopus','org.Pf.plasmo.db'=NA)
      
      print(organismR)
     viewPathway(gseaPath$Description[as.numeric(PATHID2())],
             readable = TRUE,
             organism=organismR,
             foldChange = kge())

    })

 
  })
}