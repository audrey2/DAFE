#' Normalisation_AD UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_Normalisation_AD_ui <- function(id){
  ns <- NS(id)
  library(shinyBS)
  tagList(
    box(width=12, status='info',title = h1('Settings',icon('cogs')),solidHeader = TRUE,collapsible=TRUE,

      fluidRow(
        column(width=3,
          selectInput(ns("choixTest"),"Choose the test",choices=c("Wald","LRT")),
          bsTooltip(ns("choixTest"), title = "Hypothese Testing <br> Wald Test ( default) pairwise comparisons <br> LRT which is used to identify any genes that show change in expression across the different levels "),
        ),
        column(width=3,
          selectInput(ns("choixFitType"),label=h4("Choose the Fit type",
                                                  bsButton(ns("q1"), label = "", icon = icon("question"), style = "info", size = "extra-small"))


                                                  ,choices=c('local','parametric','mean'),


                      ),

          bsTooltip(ns("q1"),title="The type of fitting of dispersions to the mean intensity <br>parametric( default )fit a dispersion-mean relation a robust gamma-family GLM <br>local  fit a local regression of log dispersions over log base mean <br>mean - use the mean of gene-wise dispersion estimates"),
        ),
        column(width=3,
          uiOutput(ns("cond1")),
        ),
        column(width=3,
          selectInput(ns("choixPadj"),"choose the methode of padj", choices=c('holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none')),
          bsTooltip(ns("choixPadj"),title="the method to use for adjusting pvalue")
        ),
      column(width=9),
      column(width=3,
        actionButton(ns("goPlot"),"Start",class="buttS",icon("play"))
      ),
    ),
  ),

  tabBox(width=12,height="600px",
  id="tabBox",
    tabPanel("Plot ",

      box(width=4, class="VIOLET",style='margin-bottom:10px;',status='success',title = h4('Boxplot of raw counts',icon('chart-simple')),solidHeader = TRUE,plotlyOutput(ns("boxplot"))%>% withSpinner()),
      box(width=4, class="VIOLET",style='margin-bottom:10px;',status='success',title = h4('Boxplot of normalized counts',icon('chart-simple')),solidHeader = TRUE,plotlyOutput(ns("boxplotn"))%>% withSpinner()),
      box(width=4, class="VIOLET",style='margin-bottom:10px;',status='success',title = h4('Graph of PCA',icon('chart-simple')),solidHeader = TRUE,plotlyOutput(ns("pcaVsd"))%>% withSpinner()),

    ),
    tabPanel("Table",
      box(width=12,class="GREEN",status='primary',title = h1('Table of normalized count',icon('table')),solidHeader = TRUE, DT::dataTableOutput(ns("tableNorm"))%>% withSpinner()),
      box(width=12,class="GREEN",status='primary',title=h1("Table of differential analysis results",icon('table')),solidHeader = TRUE, DT::dataTableOutput(ns("tableResult"))%>% withSpinner()),

    ),
    tabPanel("Dispersion",
       mod_Dispersion_analysis_ui("Dispersion_analysis_1")
    )
  )
)
}

#' Normalisation_AD Server Functions
#'
#' @noRd
mod_Normalisation_AD_server <- function(id,inputInfo,inputReplicat,DDS,TAB_RES){
  library(DESeq2)
  global<- reactiveValues()

  mod_Dispersion_analysis_server("Dispersion_analysis_1",TAB_RES,input)

  moduleServer( id, function(input, output, session){

    ns <- session$ns
    mod_Dispersion_analysis_server("Dispersion_analysis_1",TAB_RES,input)


  # Fonction renvoyant à l'interface utilisateurla liste des conditions
    output$cond1 <- renderUI({

      choix=heatCondition()
      L=c()
      L[[1]]= selectInput(ns("selectCond1"), label = "Choose condition to observe ",
                          choices = choix[2:length(choix)],
                          selected=1)

      return(L)

    })
    # Fonction renvoyant la liste des conditions
    heatCondition <- function(){
      num  = c()
      name = c()

      for (i in 1:inputInfo$nb_facteur) {
        num[i] = i
        name[i] = inputInfo[[paste0('condition', i)]]
      }

      choiceTable = data.frame(name, num)
      choix = setNames(as.numeric(choiceTable$num), choiceTable$name)

      return(choix)

    }

    # Fonction renvoyant le boxplot des comptages brutes
    boxPlot <- eventReactive(input$goPlot, {

      withProgress(message = "Plotting brut boxplot ...",{
        data = tabData(inputInfo)
        data = data[, as.numeric(inputReplicat$selectKeep) + 1]
        data = melt(data)

        if(inputInfo$exemple != 0) {data = data[2:3]}

        colnames(data) = c('condition','count')
        data['count'] = as.numeric(unlist(data['count']))
        data['count'] = log10(data['count'])
        yp = subset(data, log10(count)>0)
        sts = boxplot.stats(yp$count)$stats

        plot = ggplot(data, aes(x=condition, y=count)) +
          geom_boxplot(outlier.shape = NA) +
          coord_cartesian(ylim = c(sts[2]/2, max(sts)*1.05))

      })
      return(plot)
    })



# Fonction renvoyant le boxplot des comptages normalises
    boxPlotNorm <- eventReactive(input$goPlot, {

      dds = DDS

      withProgress(message = "Normalization ...",{
        dds2 = estimateSizeFactors(dds)
        normalized_counts = counts(dds2, normalized = TRUE)
      })

      dataNorm = melt(normalized_counts)
      dataNorm = dataNorm[c(2,3)]

      withProgress(message = "Plotting Boxplot Normalized ...", {
        colnames(dataNorm) = c('condition', 'count')
        dataNorm['count'] = log10(dataNorm['count'])
        plot = ggplot(dataNorm, aes(x = condition, y = count)) +
          geom_boxplot()

        return(plot)
      })
    })

    # Fonction renvoyant le plot de l'ACP par une transforme VST des comptages normalises
    pcaVSD <- eventReactive(input$goPlot,{

      dds = DDS
      withProgress(message = "Plotting PCA plot ...", {

        vsd = varianceStabilizingTransformation(dds, blind = TRUE)
        data = plotPCA(vsd, intgroup = c("condition"), returnData = TRUE)
        percentVar = round(100 * attr(a, "percentVar"))

        plotPCA = ggplot(data, aes(PC1, PC2, color = condition)) +
          geom_point(size = 3) +
          xlab(paste0("PC1: ", percentVar[1], "% variance")) +
          ylab(paste0("PC2: ", percentVar[2], "% variance"))

        return(plotPCA)
      })
    })


    # Rendu sur l'interface utilisateur des plot
    output$boxplotn <- renderPlotly({ boxPlotNorm() })
   output$boxplot <- renderPlotly({ boxPlot() })

   output$pcaVsd <- renderPlotly({ pcaVSD() })

    global$dds <- reactive(DDS())

    #Fonction renvoyant la table des comptages normalisés
    tabNorm <-eventReactive(input$goPlot,{
      withProgress(message = "Normalization ...", {

        factor = estimateSizeFactors(DDS)
        normTable = counts(factor, normalized = TRUE)
        name = rownames(normTable)
        normTable = cbind(name, normTable)
        rownames(normTable) = NULL

        return(normTable)
      })
    })


    #### ------render Table -------
    output$tableNorm<- DT::renderDataTable(

      DT::datatable({data.frame(tabNorm())},
                    extensions = 'Buttons',
                    caption="Table: All counts normalized",
                    options = list(

                      dom = 'Bfrtip',
                      buttons = c('copy', 'csv', 'excel')
                    ),

                    class = "display"
      ))
  # Fonction renvoyant à l'interface utilisateurs la table de resultat de l'analyse diffrentielle
    output$tableResult<- DT::renderDataTable({
      table = results(TAB_RES,
                      pAdjustMethod = input$choixPadj,
                      name = resultsNames(TAB_RES)[as.numeric(isolate(input$selectCond1))])
      name = rownames(table)
      table = cbind(name, table)
      rownames(table) = NULL

      return(DT::datatable({data.frame(table)},
                           extensions = 'Buttons',
                           caption="Table: Analysis differential table",
                           options = list(

                             dom = 'Bfrtip',
                             buttons = c('copy', 'csv', 'excel')
                           ),

                           class = "display"
      ))

    })

      # Fonction renvoyant les condition utilise pour l'analyse differentielle
    output$tabInformation <- renderText({
      num = isolate(input$selectCond1)
      if(is.null(TAB_RES)) { return("") }
      else { return(resultsNames(TAB_RES)[as.numeric(num)]) }

    })






  return(input)

  })


}

## To be copied in the UI
# mod_Normalisation_AD_ui("Normalisation_AD_1")

## To be copied in the server
# mod_Normalisation_AD_server("Normalisation_AD_1")