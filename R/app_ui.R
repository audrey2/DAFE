#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd

app_ui <- function(request) {
  library(shinydashboard)
  library(shinyWidgets)
  tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),
    # List the first level UI elements here
    dashboardPage(
      dashboardHeader(tags$li(class = "dropdown",
                              tags$style(".main-header {max-height: 100px}"),
                              tags$style(".main-header .logo {height: 100px}")),


        title = h3(HTML("Analyse différentielle et <br/>Enrichissement fonctionnel")),titleWidth = 275

      ),
      dashboardSidebar(
        width=275,
        sidebarMenu(
          br(),br(),
          h2(HTML("&ensp; BULK RNA-SEQ")),
          menuItem("Home",tabName = "home"),
          h3(HTML("&ensp;Differential Analysis")),
          menuItem("Input",tabName="input-da"),
          menuItem("Replica quality control",tabName = "ReplicaQuality"),
          menuItem("Normalization and Differential Analysis",tabName = "NormalisationAD"),
          h3(HTML("&ensp;  Enrichment Analysis ")),
          menuItem("Input",tabName = "input-ea"),
          menuItem("Enrichment Analysis",tabName = "eAnalysis",
            menuSubItem("Go terms",tabName = "ea-go"),
            menuSubItem("Pathway Kegg",tabName = "ea-kegg"),
            menuSubItem("Pathway Reactome",tabName = "ea-reactome")

                   )
        )
      ),
      dashboardBody(
        tabItems(
          tabItem(tabName = "input-da",
                  mod_Input_Norm_ui("Input_Norm_1")
                  ),
          tabItem(tabName = "ReplicaQuality",
                  mod_Replica_Quality_ui("Replica_Quality_1")
                  ),
          tabItem(tabName="NormalisationAD",
                  mod_Normalisation_AD_ui("Normalisation_AD_1")
                  ),
          tabItem(tabName="input-ea",
                  mod_parameter_EA_ui("parameter_EA_1")
                  ),
          tabItem(tabName="ea-go",
                  mod_GSEA_ui("GSEA_1"),
                  mod_ORA_ui("ORA_1")
                  ),
          tabItem(tabName="ea-kegg",
                  mod_GSEAKEGG_ui("GSEAKEGG_1"),
                  mod_ORAKEGG_ui("ORAKEGG_1")
                  ),
          tabItem(tabName="ea-reactome",
                  mod_ORAReactome_ui("ORAReactome_1"),
                  mod_GSEAReactome_ui("GSEAReactome_1")
                  )
        )
      )
      )
    )

}


#' @import shiny
golem_add_external_resources <- function(){


  tags$head(
    golem::activate_js(),
    golem::favicon(),
    tags$title("golemdashboard")
    # Add here all the external resources
    # If you have a custom.css in the inst/app/www
    # Or for example, you can add shinyalert::useShinyalert() here
    #tags$link(rel="stylesheet", type="text/css", href="www/custom.css")

  )
}
