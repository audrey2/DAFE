#' home UI Function
#'
#' @description This module displays the home page
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @author Audrey BEAUFILS
#'
#' @importFrom shiny NS tagList 
mod_home_ui <- function(id){
  ns <- NS(id)
  tagList(
    fluidRow(
    box( width=12, title=h1(strong('Presentation'),icon('display')),status = 'success',solidHeader = TRUE,collapsible = TRUE,
         
         fluidRow(
           column(width=12,
            strong(" Author :     ",style="font-size:18px;"),"Audrey Beaufils",br(),

            strong("Date :     ",style="font-size:18px;")," March 2022",
            br(),
            strong("Contact :     ",style="font-size:18px;")," audrey.beaufils@univ-rouen.fr",
              br(),strong("GitHub :     ",style="font-size:18px;"),a("Audrey2", href="https://github.com/audrey2/"),
              br(),
             strong("Description Application :      ",style="font-size:18px;"),"DAFE is web RShiny application which is able to run differential analysis and fonctional enrichment on",a("GO terms",href="http://geneontology.org/"),
              " Pathways",a("KEGG",href="https://www.genome.jp/kegg/")," and Pathways",a("Reactome",href="https://reactome.org/")


            ),
           
         ),
    )),
    br(),
    fluidRow(
    box(width=6,title=h1(strong("Differential Analysis Module"),icon('file-code')),status = 'info', solidHeader=TRUE,

      fluidRow(column(width=12,h3('Input specificity'),

              "The input must have this form, without header",
              tableOutput(ns('input_DA')),
              br(),
              h3('Packages necessary'),
              "Normalisation an Differential Analysis use packages : ", a("DESeq2",href="https://bioconductor.org/packages/release/bioc/html/DESeq2.html"),
              br(),
              "The display of file use packages :",a("ggplot",href="https://cran.r-project.org/web/packages/ggplot2/index.html"),a("ggpubr",href="https://cran.r-project.org/web/packages/ggpubr/index.html"),a("ggplotly",href="https://cran.r-project.org/web/packages/plotly/index.html"),
              br(),
              h3("Output specificity"),
              "All plot and all table are dowloadable",br(),
              "The normalized counts table has the following forms:",br(),
              tableOutput(ns('output_DA_norm')),
              "The differential analysis table has the following forms:",br(),
              tableOutput(ns('output_DA_DA'))




        )
      )


    ),
    box(width=6,title=h1(strong("Fonctionnal Enrichment Module"),icon('file-code')),status = 'info', solidHeader=TRUE,

      fluidRow(column(width=12,h3('Input specificity'),
      "The input must have this form, it must be have more columns but at least these names",
          tableOutput(ns('input_EF')),
          h3('Packages necessary'),
              "Gene Set Enrichment and Over Representation Analysis use packages : ", a("clusterProfiler",href="https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html"),
              br(),
              "The Reactome pathways enrichment use packages : ", a("ReactomePA",href="https://bioconductor.org/packages/release/bioc/html/reactomePA.html"),
              br(),
              "The display of file use packages :",a("ggplot",href="https://cran.r-project.org/web/packages/ggplot2/index.html"),a("enrichplot",href="https://bioconductor.org/packages/release/bioc/html/enrichplot.html"),a("ggplotly",href="https://cran.r-project.org/web/packages/plotly/index.html"),
              br(),
              h3('Output specificity '),
              "All plot and all table are dowloadable",br(),
              "GSEA tables have the following forms",br(),
              tableOutput(ns('output_GSEA')),
              br(),
              "ORA tables have the following forms",br(),
              tableOutput(ns('output_ORA')),


  )))))
}
    
#' home Server Functions
#'
#' @noRd 
mod_home_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns



    output$input_DA <- renderTable({
      
      tab= matrix(ncol=2)
      colnames(tab)=c('Gene_name','Counts')
      tab
      }
    )

     output$output_DA_norm <-renderTable({
        
      tab= matrix(nrow=1,ncol=5)
      colnames(tab)=c('name','Counts1_R1_norm','Counts1_R2_norm','Counts2_R1_norm','Counts2_R2_norm')
      tab
    })

      output$output_DA_DA <- renderTable({
      
       tab= matrix(nrow=1,ncol=7)
      colnames(tab)=c('name','baseMean','log2FoldChange','lfcSE','stat','pvalue','padj')
      tab
      }
    )
      output$input_EF <- renderTable({
      
       tab= matrix(nrow=1,ncol=6)
      colnames(tab)=c('name','baseMean','log2FoldChange','stat','pvalue','padj')
      tab
      }
    )

    output$output_GSEA <- renderTable({
      
       tab= matrix(nrow=1,ncol=7)
      colnames(tab)=c('ID','Description','enrichmentScore','NES','pvalue','p.adjust','qvalues')
      tab
      }
    )
      
      output$output_ORA <- renderTable({
      
       tab= matrix(nrow=1,ncol=8)
      colnames(tab)=c('ID','Description','GeneRatio','BgRatio','pvalue','p.adjust','qvalue','geneID')
      tab
      }
    )
  })
}
    
## To be copied in the UI
# mod_home_ui("home_1")
    
## To be copied in the server
# mod_home_server("home_1")
