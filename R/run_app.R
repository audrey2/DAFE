#' Run the Shiny Application
#'
#' @param ... arguments to pass to golem_opts.
#' See `?golem::get_golem_options` for more details.
#' @inheritParams shiny::shinyApp
#'
#' @export
#' @importFrom shiny shinyApp
#' @importFrom golem with_golem_options
run_app <- function(
  onStart = function(){
   print('Install necessary package') 
    if(!"BiocManager" %in% installed.packages()){
      install.packages("BiocManager")
    }
    #BiocManager::install("DESeq2",update=FALSE)
    #BiocManager::install("clusterProfiler",update=FALSE)
    #BiocManager::install("DOSE",update=FALSE)
    #BiocManager::install("ReactomePA",update=FALSE)
    #BiocManager::install("enrichplot",update=FALSE)
    #BiocManager::install("pathview",update=FALSE)
    #BiocManager::install("GO.db",update=FALSE)
    #BiocManager::install("ggpubr",update=FALSE)
    #BiocManager::install("heatmaply",update=FALSE)
  },
  options = list(),
  enableBookmarking = NULL,
  uiPattern = "/",
  ...
) {
  with_golem_options(
    app = shinyApp(
      ui = app_ui,
      server = app_server,
      onStart = onStart,
      options = options,
      enableBookmarking = enableBookmarking,
      uiPattern = uiPattern
    ),
    golem_opts = list(...)
  )
}
