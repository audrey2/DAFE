#' home UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_home_ui <- function(id){
  ns <- NS(id)
  tagList(
    box( width=4, title=h1('Presentation',icon('cogs')),status = 'info',solidHeader = TRUE,collapsible = TRUE,
         
         fluidRow(
           column(width=12," Author:  Audrey Beaufils",br(),"Date : mars2022",br(),"Contact : audrey_beaufils@hotmail.fr"),
           
         ),
    ),
  )
}
    
#' home Server Functions
#'
#' @noRd 
mod_home_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns
 
  })
}
    
## To be copied in the UI
# mod_home_ui("home_1")
    
## To be copied in the server
# mod_home_server("home_1")
