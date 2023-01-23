#' input_parameter_norm UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_input_parameter_norm_ui <- function(id){
  ns <- NS(id)

}

#' input_parameter_norm Server Functions
#'
#' @noRd
mod_input_parameter_norm_server <- function(id){
  moduleServer( id, function(input, output, session){
    ns <- session$ns


  print(input)
  return(input)
  })
}

## To be copied in the UI
# mod_input_parameter_norm_ui("input_parameter_norm_1")

## To be copied in the server
# mod_input_parameter_norm_server("input_parameter_norm_1")
