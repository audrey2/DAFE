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
    dashboardPage(skin="purple",
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
          menuItem("Home",tabName = "home",icon=shiny::icon('home')),
          h3(HTML("&ensp;Differential Analysis"),icon=icon('magnifying-glass')),
          menuItem("Input",tabName="input-da",icon=icon('upload')),
          menuItem("Replica quality control",tabName = "ReplicaQuality",icon=icon('vial-circle-check')),
          menuItem("Normalization and Differential Analysis",tabName = "NormalisationAD",icon=icon('arrows-rotate')),
          h3(HTML("&ensp;  Enrichment Analysis "),icon=icon('magnifying-glass')),
          menuItem("Input",tabName = "input-ea",icon=icon('upload')),
          menuItem("Enrichment Analysis",tabName = "eAnalysis",icon=icon('arrows-rotate'),
            menuSubItem("Go terms",tabName = "ea-go",icon=icon('network-wired')),
            menuSubItem("Pathway Kegg",tabName = "ea-kegg",icon=icon('network-wired')),
            menuSubItem("Pathway Reactome",tabName = "ea-reactome",icon=icon('network-wired'))

                   )
        )
      ),
      dashboardBody(height='2000px',
            tags$style(HTML("
              
              .box.box-solid.box-info>.box-header {
                color:#ffffff;
                background:#CDCDE6;
              }
              .box.box-solid.box-info{
                border-bottom-color:#666666;
                border-left-color:#666666;
                border-right-color:#666666;
                border-top-color:#666666;
              }
              .box.box-solid.box-success>.box-header {
                color:#ffffff;
                background:#605CA8;
              }
              .box.box-solid.box-success{
                border-bottom-color:#666666;
                border-left-color:#666666;
                border-right-color:#666666;
                border-top-color:#666666;
              }
            ")),
        tabItems(
          tabItem(tabName="home",
                  mod_home_ui("home_1")
        ),
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

  )
}
