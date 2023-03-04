#' Input_Norm UI Function
#'
#' @description This module displays parameter for differential Analysis
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @author Audrey BEAUFILS
#'
#' @importFrom shiny NS tagList
mod_Input_Norm_ui <- function(id){
  ns <- NS(id)
  library(shinydashboard)
  library(shinyWidgets)
  tagList(


   box(title=h1("Settings",icon('gear')),width=12,solidHeader=TRUE, collapsible = TRUE,



     actionButton(ns("exemple"),label="Use data test",icon("file"),class="butt"),


     uiOutput(ns("title")),


     uiOutput(ns("browseFiles")),
     hr(),
     fluidRow(column(width=10),column(width=2,
     actionButton(ns("check"),label="View/Check File",icon("check"),class="butt"))),

      status = "info",
     icon = icon("gear"))
   ,



    uiOutput(ns("table_file"))

    )


}

#' Input_Norm Server Functions
#'
#' @noRd
mod_Input_Norm_server <- function(id,inputNorm){
  moduleServer( id, function(input, output, session){

    ns <- session$ns
    library(ggplot2)
    library(stringr)
    library(aod)
    library(sgof)
    library(stats)
    library(plotly)
    library(reshape2)
    library(shinyBS)




    output$title <- renderUI({
      L=c()
      if(input$exemple== FALSE) {
        L[[1]]=   h3(HTML('&ensp;or'))
        L[[2]]=hr()
        L[[3]]= column(width=4, textInput(ns("title"), label="Write the title of Experience", value="My Experience"),
    numericInput(ns("nb_facteur"), "Number of condition", min = 1, max = 10, step = 1, value = 2))
      }
      return(L)
    })

    # Fonction affichant l'interface utilisateur en fonction des paramètre spécifié
    output$table_file <- renderUI({
    L=c()
        L[[1]]= box(width=12,style='margin-top: 10px;padding: 10px;',status='primary',title = h1('Table : View of all counts concatenate',icon('table')),solidHeader = TRUE,

                    DT::dataTableOutput(ns("viewAllComptage"))
        )

      return(L)
    })


    AllComptage <-reactive( {

      if(input$check != 0){

        if(tabCheck(input) == 1 ) {
          data=data.frame(tabData(input))
          return(data)
        }
      }
      return(NULL)
    } )

    #Fonction renvoyant un datatable del'intégralité des comptages concatener dans une seule table'
    output$viewAllComptage <- DT::renderDataTable(
      DT::datatable(
        { AllComptage()},
        extensions = 'Buttons',
        options = list(
          dom = 'Bfrtip',
          buttons = c('copy', 'csv', 'excel')
        ),class = "display"
      )
    )

    #Fonction renvoyant le nombre de bowfiles en fonction du nombre de condition spécifié
    output$browseFiles <- renderUI({

      boxFiles = c()
      if(input$exemple== FALSE) {
      for(i in 1:input$nb_facteur) {

        boxFiles[[i]] = column(width=4,
          if(i==1) { textInput(ns(paste0("condition", i)), label=" Write Condition Name",value = paste0('Control')) }

          else{ textInput(ns(paste0("condition",i)),label=" Write Condition Name",value= paste0('Condition',i-1)) },

          bsTooltip(id = "box1", title = "The files in the first box must be the controle condition"),
          fileInput(ns(paste0("file", i)), label="Choose count Files", accept = c(".csv",".txt",".count",".tsv",".counts"), multiple=TRUE)
        )
      }}
      return(boxFiles)
    })






    return(input)

  })
}

## To be copied in the UI
# mod_Input_Norm_ui("Input_Norm_1")

## To be copied in the server
# mod_Input_Norm_server("Input_Norm_1")
