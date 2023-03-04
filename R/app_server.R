#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#'@author Audrey BEAUFILS
app_server <- function(input, output, session) {




  # Fonction renvoyant la liste des conditions
  heatCondition <- function(){
    req(inputInfo$nb_facteur)
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



  # Fonction qui renvoie un objet DESeq DataSet
   DDS <- reactive({
    req(inputInfo,inputReplicat$selectKeep,inputInfo$exemple)
    data = tabData(inputInfo)
    row.names(data) = data[, 1]
    data = data[, as.numeric(inputReplicat$selectKeep) + 1]
    data = as.matrix(data,)
    replica = getReplicaNames(inputInfo)[as.numeric(inputReplicat$selectKeep)]
    condition = heatCondition()
    replicatName = c()
    conditionName = c()

    if(inputInfo$exemple == 0) {

      for( i in 1: length(replica)) {
        for(j in 1: inputInfo$nb_facteur) {

          sub = substr(replica[i], 1, nchar(names(condition[j])))

          if(names(condition[j]) == sub) {

            conditionName[i] = names(condition[j])
            replicatName[i] = replica[i]
          }
        }
      }

      coldata = cbind(conditionName,replicatName)
      coldata = data.frame(coldata,row.names = 2)
      colnames(coldata) = c('condition')

    }
    else { # Pour les donnees test

      coldata = c("Control", "Control", "Control", "Condition1", "Condition1", "Condition1")
      coldata = cbind(coldata,c("Control R1", "Control R2", "Control R3", "Condition1 R1", "Condition1 R2", "Condition1 R3"))
      coldata = data.frame(coldata, row.names = 2)
      colnames(coldata) = c('condition')
      class(data) = "numeric"

    }

    coldata$condition = factor(coldata$condition, levels = unique(coldata$condition))

    mat=data


    ddds = DESeqDataSetFromMatrix(countData = mat,
                                  colData = coldata,
                                  design = ~condition)



    return(ddds)
  })

  # Fonction renvoyant le tabelau de rsultat de l'analyse differentielle
  tabRes <-reactive( {
    req(inputNorm$choixFitType,inputNorm$choixTest)
    
    withProgress(message = "Differential Analysis ..." ,{
      dds = DDS()

      if(inputNorm$choixTest == "LRT") {
        table = DESeq(dds,fitType = inputNorm$choixFitType, test = inputNorm$choixTest, reduced=~1)
      }
      else {
        table = DESeq(dds,fitType = inputNorm$choixFitType,test = inputNorm$choixTest)
      }
      return(table)
    })
  })

  mod_home_server("home_1")
  inputInfo=mod_Input_Norm_server("Input_Norm_1")
  inputReplicat=mod_Replica_Quality_server("Replica_Quality_1",inputInfo)
  inputNorm=mod_Normalisation_AD_server("Normalisation_AD_1",inputInfo,inputReplicat,DDS(),tabRes())
  inputParameter=mod_parameter_EA_server("parameter_EA_1")




  mod_GSEA_server("GSEA_1",inputParameter)
  mod_ORA_server("ORA_1",inputParameter)
  mod_ORAKEGG_server("ORAKEGG_1",inputParameter)
  mod_GSEAKEGG_server("GSEAKEGG_1",inputParameter)
  mod_GSEAReactome_server("GSEAReactome_1",inputParameter)
  mod_ORAReactome_server("ORAReactome_1",inputParameter)




}
