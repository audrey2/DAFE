#' Replica_Quality UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_Replica_Quality_ui <- function(id){
  ns <- NS(id)
  library(shinycssloaders)
  library(plotly)
  library(colourpicker)
  tagList(
    box( width=12, title=h1('Setting',icon('cogs')),status = 'info',solidHeader = TRUE,collapsible = TRUE,

      fluidRow(
        column(width=3, uiOutput(ns("choicesCondition"))),
        column(width=3, uiOutput(ns("keepingReplica"))),
        column(width=3,br(),br(),
        actionButton(ns("goHeat"),"Start",class="buttS",icon("play")))
      ),
    ),
    br(),br(),br(),br(),
    box(width=12,title=h1('Heatmap',icon('chart-simple')),status = 'success',solidHeader = TRUE,height='800px',
      br(),
      column(width=3,
      dropdownButton(inline=TRUE,icon = icon('gear'),status = 'warning',
                     tags$h3("Personnalize"),
                     column(width=12, colourInput(ns("col1"), "Choose first color","#ade6f4")),
                     column(width=12, colourInput(ns("col2"), "Choose second color","#45d7b7")),
                     actionButton(ns("goHeat2"),"Start",class="buttS",icon("play")))),
      br(),
      column(width=6,      plotlyOutput(ns("heatMap"))),
      column(width=3)
    )
  )
}


#' Replica_Quality Server Functions
#'
#' @noRd
mod_Replica_Quality_server <- function(id,inputNorm){

  moduleServer( id, function(input, output, session){
    ns <- session$ns
    library(plotly)
    library(heatmaply)

    heatCondition <- reactive({
      num  = c()
      name = c()

      for (i in 1:inputNorm$nb_facteur) {
        num[i] = i
        name[i] = inputNorm[[paste0('condition', i)]]
      }

      choiceTable = data.frame(name, num)
      choix = setNames(as.numeric(choiceTable$num), choiceTable$name)

      return(choix)
    })


    # Fonction qui renvoie tout les replicats des conditions garder par l'utilisateur
    heatReplica <- reactive({

      conditionKeep = heatCondition()
      replicaNames = getReplicaNames(inputNorm)
      numReplica = c()

      for( i in 1: length(replicaNames)) {
        for(j in 1: inputNorm$nb_facteur) {
          sub = substr(replicaNames[i], 1, nchar(names(conditionKeep[j])))
          if(names(conditionKeep[j]) == sub) { numReplica[i] = conditionKeep[j] }
        }
      }

      num  = c()
      name = c()
      cpt = 1

      for (i in 1:length(replicaNames)) {
        if(numReplica[i] %in% input$select1) {
          num[cpt] = toString(i)
          name[cpt] = replicaNames[i]
          cpt = cpt + 1
        }
      }

      choiceTable = data.frame(name, num)
      choix = setNames(as.numeric(choiceTable$num), choiceTable$name)

      return( choix )

    })

    # Fonction qui renvoie l'integralite des conditions de l'experience sous forme d'un selectInput
    output$choicesCondition <- renderUI({

      choix=heatCondition()
      L=c()

      L[[1]]= selectInput(ns("select1"), label = h4("Select condition to observe"),
                          choices = choix, selected = choix,
                          multiple = TRUE)
      return(L)

    })

    # Fonction qui renvoie l'integralite des replicats de l'experience sous forme d'un selectInput
    output$keepingReplica <- renderUI({

      if(inputNorm$exemple==0) { choix = heatReplica() }
      else {
        choix = c("Control R1", "Control R2", "Control R3", "Condition1 R1", "Condition1 R2", "Condition1 R3")
        num = c(1, 2, 3, 4, 5, 6)
        choix = setNames(as.numeric(num), choix)
      }
      L = c()
      L[[1]] =  selectInput(ns("selectKeep"), label = h4("Select Replica to keep"),
                            choices = choix, selected = choix, multiple = TRUE)

      return(L)
    })


    # Fonction qui creer une heatmap de distance après go de l'utilisateur
    heat<-eventReactive((input$goHeat | input$goHeat2),{
      validate(
        need(input$selectKeep != "", "Please click on start")
      )
      withProgress(message = "Plotting heatMap ...", {

        chosenReplica = as.numeric(input$selectKeep) +1
        data = tabData(inputNorm)[, chosenReplica]

        dataT = t(data.frame(data))
        print(input$col1)
        print(input$col2)
        cols = colorRampPalette(c(input$col1, input$col2))(7)
        distance = dist(dataT, method = "euclidian")
        matrice = as.matrix(distance)

        figure = heatmaply(matrice,colors=cols)%>% layout(height = 700, width = 800)

      })

      return(figure)
    })
    

    
    # Affichage UI Heatmap
    output$heatMap <- renderPlotly(heat())



  return(input)
  })
}


## To be copied in the UI
# mod_Replica_Quality_ui("Replica_Quality_1")

## To be copied in the server
# mod_Replica_Quality_server("Replica_Quality_1")
