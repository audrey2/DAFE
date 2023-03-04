#' parameter_EA UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_parameter_EA_ui <- function(id){
  ns <- NS(id)
  library(shinyFiles)
  tagList(
    box( width=12, status = 'info',solidHeader = TRUE,collapsible = TRUE,title=h1('Settings',icon('gear')),

fluidRow(

  column(width=3,fileInput(ns("fileOr"),"BrowseFile"),numericInput(ns("row.names"), "Column number of gene ID", min = 1, max = 10, step = 1, value = 1)),
  column(width=3, selectInput(ns("sep"),"Select the separator of files",choices=list("Comma"=',',"Tabulation"='\t',"Semicolon"=';'),selected = ","),br(), selectInput(ns("keytype"),"Select the key type",choices=c('ENSEMBL','ENTREZID'))),
  column(width=3,
         selectInput(ns("orDb"),"choose the Organism DB", choices=c('Fly'='org.Dm.eg.db','Human'='org.Hs.eg.db','Mouse'='org.Mm.eg.db','Rat'='org.Rn.eg.db',
                                                                    'Arabidopsis'='org.At.tair.db','Yeast'='org.Sc.sgd.db',
                                                                    'Zebrafish'='org.Dr.eg.db','Worm'='org.Ce.eg.db','Bovine'='org.Bt.eg.db',
                                                                    'Pig'='org.Ss.eg.db','Chicken'='org.Gg.eg.db','Rhesus'='org.Mmu.eg.db',
                                                                    'Canine'='org.Cf.eg.db','E coli strai'='org.EcK12.eg.db','Xenopus'='org.Xl.eg.db',
                                                                    'Anopheles'='org.Ag.eg.db','Chimp'='org.Pt.eg.db','Malaria'='org.Pf.plasmo.db',
                                                                    'E coli strai'='org.EcSakai.eg.db','Myxococcus x'='org.Mxanthus.db')),
         br(),
         selectInput(ns("Gsea_ORA")," Choose Method",choices=c("GSEA","ORA"))),
 #
  column(width=3,

         selectInput(ns("pAdjMethod"),"Choose Padj Methode", choices=c( "BH","holm", "hochberg", "hommel", "bonferroni", "BY", "fdr","none" )),
         br(),
         uiOutput(ns("subsetData"))
         
         ),


)),
box(width=12,status='primary',solidHeader = TRUE,title=(h1("Table",icon('table'))),
  DT::dataTableOutput(ns("file")))

)

}

#' parameter_EA Server Functions
#'
#' @noRd
mod_parameter_EA_server <- function(id){
  moduleServer( id, function(input, output, session){
    library(shinyFiles)
    library(DT)
    ns <- session$ns
    df <-reactive({
      # reading in data from deseq2

      req(input$fileOr)
      df = read.csv(input$fileOr[1,'datapath'], header=TRUE,sep=input$sep)


      return(df)

    })

    output$subsetData<- renderUI({
      L=c()
        if(input$Gsea_ORA=='ORA'){
          L=selectInput(ns("ora_order"),"Choose data to keep", choices=c("both","underexpressed","overexpressed" ))
          }
      else{L=selectInput(ns("gsea_order"),"Choose how sort data", choices=c( "pval","log2FoldChange", "Stat" ))}
      return(L)
    })
          
    # Fonction renvoyant à l'ui le tableau charger
    output$file <- DT::renderDataTable(DT::datatable(
      { df()},

      extensions = 'Buttons',
      caption="Table: View differential analysis table load",
      options = list(

        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel')
      ),

      class = "display"
    ))

  return(input)

  })
}

## To be copied in the UI
# mod_parameter_EA_ui("parameter_EA_1")

## To be copied in the server
# mod_parameter_EA_server("parameter_EA_1")
