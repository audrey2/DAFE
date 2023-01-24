#' Dispersion_analysis UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_Dispersion_analysis_ui <- function(id){
  ns <- NS(id)
  tagList(

    #column(width=3,



     #      box(width=12,status='info',solidHeader=TRUE,title=h3("Experience Information",icon('info')),

      #         verbatimTextOutput(ns("infoDiff")))


  #  ),
    column(width=12,
           box(width=6,status='success',solidHeader=TRUE,title=h1("Volcano Plot",icon('chart-simple')),
  fluidRow(column(width=1,
               dropdownButton(inline=TRUE,icon = icon('gear'),status = 'warning',width="300px",
                              tags$h3("Option of treshold"),
                              sliderInput(ns("ts_FC"), label = "log2 FoldChange  cutoff from input",
                                          min = 0, max = 5, value = 1,step=0.1),
                              sliderInput(ns("ts_padj"),label = "p-Value cutoff from output ",
                                          min = 0, max =0.1, value = 0.01,step=0.01))),column(width=11,

               plotlyOutput(ns("volcano"))%>% withSpinner()))),

        box(width=6,status='success',solidHeader=TRUE,title=h1("MA plot",icon('chart-simple')),plotlyOutput(ns("maplot"))%>% withSpinner()),
   ),
   box(width=12,status='primary',title = h1('Table of selected results',icon('table')),solidHeader = TRUE, DT::dataTableOutput(ns("degTableZoomed"))%>% withSpinner()),


  )
}

#' Dispersion_analysis Server Functions
#'
#' @noRd
mod_Dispersion_analysis_server <- function(id,TAB_RES,inputNorm){
  loadNamespace("DT")

  moduleServer( id, function(input, output, session){
    ns <- session$ns

    RES <- reactive({

      resD=  results(TAB_RES)
    name=rownames(resD)

    res2=cbind(name,resD)
    rownames(res2)=NULL
    df=data.frame(res2)

    return(df)
    })
    ranges <- reactiveValues(x = NULL, y = NULL)

    # Fonction generant le volcano plot
    Volcano <- reactive({


      withProgress(message = "Plotting Volcano ...",{
      tsPadj = as.numeric(input$ts_padj)
      tsFC = as.numeric(input$ts_FC)
      table = RES()

      table$diffexpressed = "NO"

      table$diffexpressed[table$log2FoldChange > tsFC & table$padj < tsPadj] = "UP"

      table$diffexpressed[table$log2FoldChange < -tsFC & table$padj < tsPadj] = "DOWN"

      cols=c("gold", "black","blue")

      if(nrow(subset(table,diffexpressed = 'NO'))   == 0) { cols = c('gold' , 'blue')}
      if(nrow(subset(table,diffexpressed = 'DOWN')) == 0) { cols = c('black', 'blue')}
      if(nrow(subset(table,diffexpressed = 'UP'))   ==0)  { cols = c('black', 'gold')}

      plot = ggplot(data=table, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, tooltip=name)) +
        geom_point(size=1) + theme_minimal()+
        geom_vline(xintercept=c(-tsFC, tsFC), col="red") +
        geom_hline(yintercept=-log10(tsPadj), col="red")+
        coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE)+
        scale_color_manual(values=cols)
})
      return(plot)
    })


    # Fonction generant le MA plot
    MA<- reactive({

      tsPadj = as.numeric(input$ts_padj)
      tsFC = as.numeric(input$ts_FC)
      table = RES()

      table$diffexpressed = "NO"
      table$diffexpressed[table$log2FoldChange > tsFC & table$padj < tsPadj] = "UP"
      table$diffexpressed[table$log2FoldChange < -tsFC & table$padj < tsPadj] = "DOWN"

      cols=c("gold", "black","blue")

      if(nrow(subset(table,diffexpressed = 'NO'))   == 0) { cols = c('gold' , 'blue')}
      if(nrow(subset(table,diffexpressed = 'DOWN')) == 0) { cols = c('black', 'blue')}
      if(nrow(subset(table,diffexpressed = 'UP'))   ==0)  { cols = c('black', 'gold')}


      plot = ggplot(data=table, aes(x=log2(baseMean), y=log2FoldChange, col=diffexpressed, tooltip=name)) +
        geom_point(size=1) + theme_minimal() +
        geom_hline(yintercept=c(-tsFC,tsFC), col="red") +
        scale_color_manual(values=cols)

      return(plot)

    })


    # Affichage du nombre de transcrits differentiellemnt exprimes

    output$infoDiff<- renderText({

      tsPadj = as.numeric(input$ts_padj)
      tsFC = as.numeric(input$ts_FC)
      table = RES()

      table$diffexpressed = "NO"
      table$diffexpressed[table$log2FoldChange > tsFC & table$padj < tsPadj] = "UP"
      table$diffexpressed[table$log2FoldChange < -tsFC & table$padj < tsPadj] = "DOWN"

      info = paste0("You have ",nrow(table)," transcrits :\n    - ",
                    nrow(subset(table,diffexpressed == 'NO')), " No diffexpressed \n   - ",
                    nrow(subset(table,diffexpressed == 'DOWN')), " Down regulated \n   - ",
                    nrow(subset(table,diffexpressed == 'UP')), " Up regulated \n",sep="")
      x = xlims()
      y = ylims()

      if(!is.null(x)) {

        table = subset(table,log2FoldChange > x[1] & log2FoldChange < x[2])
        table = subset(table,-log10(padj) > y[1] & -log10(padj) < y[2])

        table$diffexpressed = "NO"
        table$diffexpressed[table$log2FoldChange > tsFC & table$padj < tsPadj] = "UP"
        table$diffexpressed[table$log2FoldChange < -tsFC & table$padj < tsPadj] = "DOWN"


        info = paste0(info, "\nIn your selection you have :\n   -  ",
                      nrow(subset(table,diffexpressed == 'NO')), " No diffexpressed \n   - ",
                      nrow(subset(table,diffexpressed == 'DOWN')), " Down regulated \n   - ",
                      nrow(subset(table,diffexpressed == 'UP')), " Up regulated \n",sep = "")

      }

      return(info)
    })

    # Affichage Plot
    output$maplot  <- renderPlotly({ ggplotly(MA(), source = "MA")})
    output$volcano <- renderPlotly({ ggplotly(Volcano(), source = "Volcano")})
    xlims <- function(){
      zoom = event_data("plotly_relayout", "Volcano")

      if(is.null(zoom) || names(zoom[1]) %in% c("xaxis.autorange", "width")) {
        xlim = NULL
      }
      else {
        xmin = zoom$`xaxis.range[0]`
        xmax = zoom$`xaxis.range[1]`
        xlim = c(xmin, xmax)
      }
      return(xlim)
    }
    # Fonction renvoyant les bornes x du volcanoplot apres zoom
    ylims<- function(){
      zoom = event_data("plotly_relayout", "Volcano")

      if(is.null(zoom) || names(zoom[1]) %in% c("yaxis.autorange", "height")) {
        ylim = NULL
      }
      else {
        ymin = zoom$`yaxis.range[0]`
        ymax = zoom$`yaxis.range[1]`
        ylim = c(ymin, ymax)
      }
      return(ylim)
    }
    output$degTableZoomed <- DT::renderDataTable({

      x = xlims()
      y = ylims()
      table = data.frame(RES())

      if(!is.null(x)){
        table = subset(table, log2FoldChange > x[1] & log2FoldChange < x[2])
        table = subset(table, -log10(padj) > y[1] & -log10(padj) < y[2])
      }

      return(DT::datatable({data.frame(table)},extensions = 'Buttons',
                           caption="Table: Analysis differential table of selected transcripts",
                           options = list(

                             dom = 'Bfrtip',
                             buttons = c('copy', 'csv', 'excel')
                           ),

                           class = "display"
      ))

    })



  })
}

## To be copied in the UI
# mod_Dispersion_analysis_ui("Dispersion_analysis_1")

## To be copied in the server
# mod_Dispersion_analysis_server("Dispersion_analysis_1")
