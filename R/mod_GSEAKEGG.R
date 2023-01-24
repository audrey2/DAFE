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
    library(cowplot)
    library(magick)
    library(ggimage)
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
                     plotlyOutput(ns('png'))%>% withSpinner(),
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
     
      library(organism, character.only = TRUE)

      return(organism)
    })


    # Fonction renvoyant un objet gseKEGG
    gene_list_kegg <- reactive({
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

    return(kegg_gene_list)})

    PATHID<- reactive(input$pathiId)



    # Fonction renvoyant les plot à l'ui
    output$dotplotKegg <- renderPlotly(ggplotly( dotplot(gene_list_kegg(), showCategory = 10, title = "Enriched Pathways" ) ))
    output$gseaplotKegg <- renderPlotly(ggplotly(gseaplot(gene_list_kegg(), by = "all", title = gene_list_kegg()$Description[as.numeric(PATHID())], geneSetID =as.numeric( PATHID()))))


# Fonction renvoyant l'image KEGG à l'ui
  img=reactive({img=ggdraw() +draw_image(pathView()) 
    return(img)})

   output$png<- renderPlotly( {
   

#Constants 

img_width = 1600 

img_height = 900 

scale_factor = 0.5 



# Add invisible scatter trace. 

# This trace is added to help the autoresize logic work. 

fig <- plot_ly(width=img_width * scale_factor, 

               height=img_height * scale_factor 

) %>% 

  add_trace( x= c(0, img_width * scale_factor), 

             y= c(0, img_height * scale_factor), 

             type = 'scatter',  mode = 'markers', alpha = 0) 


# Configure axes 

xconfig <- list( 

  title = "", 

  zeroline = FALSE, 

  showline = FALSE, 

  showticklabels = FALSE, 

  showgrid = FALSE, 

  range = c(0, img_width * scale_factor) 

) 


yconfig <- list( 

  title = "", 

  zeroline = FALSE, 

  showline = FALSE, 

  showticklabels = FALSE, 

  showgrid = FALSE, 

  range = c(0, img_height * scale_factor), 

  scaleanchor="x" 

) 


fig <- fig %>% layout(xaxis = xconfig, yaxis = yconfig) 


# Add image 


fig <- fig %>% layout( 

  images = list(  

    list(  

      source = raster2uri(readPNG(pathView())),  

      x=0, 

      sizex=img_width * scale_factor, 

      y=img_height * scale_factor, 

      sizey=img_height * scale_factor, 

      xref="x", 

      yref="y", 

      opacity=1.0, 

      layer="below", 

      sizing="stretch" 

    )  

  )) 


# Configure other layout 


m = list(r=0, l=0, b=0, t=0) 

fig <- fig %>% layout(margin = m) %>%

  layout(plot_bgcolor='#e5ecf6',  

          xaxis = list(  

            zerolinecolor = '#ffff',  

            zerolinewidth = 2,  

            gridcolor = 'ffff'),  

          yaxis = list(  

            zerolinecolor = '#ffff',  

            zerolinewidth = 2,  

            gridcolor = 'ffff')  

          )

fig
    })

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
    pathview(gene.data=kge(), pathway.id=names(choix()[as.numeric(PATHID())]), species = kegg_organism, kegg.dir = dir)

    a=paste0(dir,names(choix()[as.numeric(PATHID())]))
    print(a)
    a=paste0(a,".pathview.png")
    print(a)


    return(a)

   })

  # Fonction renvoyant le tableau de resultat de l'enrichissement
  tabb <- reactive({
    table= cbind(gene_list_kegg()$ID,gene_list_kegg()$Description,gene_list_kegg()$enrichmentScore,gene_list_kegg()$NES,gene_list_kegg()$pvalue,gene_list_kegg()$p.adjust,gene_list_kegg()$qvalue,gene_list_kegg()$rank)
    colnames(table)= c('ID','Description','enrichmentScore','NES','pvalue','p.adjust','qvalues','rank')
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
