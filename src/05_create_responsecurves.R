#--------------------------------------------
#-----------To do: specify project-----------
#--------------------------------------------
#specify project name
projectname<-"Buffered_occurrences_FINAL"

#For which model do you want to create response curves?
model<-"Global"
#model<-"European"

#--------------------------------------------
#-----------  Load packages  ----------------
#--------------------------------------------
packages <- c("viridis", "dplyr", "grid", "here", "qs","terra", "sf", "ggplot2","RColorBrewer","magick","patchwork",
              "ape", "geoR", "raster", "pdp", "purrr", "gbm"
)

for(package in packages) {
  print(package)
  if( ! package %in% rownames(installed.packages()) ) { install.packages( package ) }
  library(package, character.only = TRUE)
}


#--------------------------------------------
#---  Load right version of caretEnsemble  --
#--------------------------------------------
desired_version <- "2.0.3"

# Check if caretEnsemble is installed
if ("caretEnsemble" %in% rownames(installed.packages())) {
  # Get the current version of caretEnsemble
  current_version <- packageVersion("caretEnsemble")
  # Compare current version with the desired version
  if (as.character(current_version) != desired_version) {
    # Uninstall the current version if it's not the desired version
    remove.packages("caretEnsemble")
    # Install the specific version
    devtools::install_github("zachmayer/caretEnsemble@2.0.3")
    # Load 
    library(caretEnsemble)
  } else {
    library(caretEnsemble)
  }
  
} else {
  # If caretEnsemble is not installed, install the specific version
  devtools::install_github("zachmayer/caretEnsemble@2.0.3")
  # Load 
  library(caretEnsemble)
  rm(current_version, desired_version)
}

#--------------------------------------------
#------- Source helper functions     --------
#--------------------------------------------
source("./src/helper_functions.R")



#--------------------------------------------
#--------- Create output folder -------------
#--------------------------------------------
project_path <- file.path("./data/projects",projectname)
if (model=="Global"){
  ModelOutput<-paste0(project_path,"/GlobalModelOutput")
  model<-qread(paste0(project_path,"/GlobalModelOutput", "/GlobalModelOutput.qs"))
} else {
  ModelOutput<-paste0(project_path,"/EuropeanModelOutput")
  model<-qread(paste0(ModelOutput, "/EuropeanModelOutput.qs"))
}
ensemble_model<-model$ensemble_model
ResponseCurvesOutput<-paste0(ModelOutput,"/ResponseCurvesOutput" )

if (!dir.exists(ResponseCurvesOutput)) {
  dir.create(ResponseCurvesOutput, recursive = TRUE)
}



### Get variable importance of best european model

variableImportance<-varImp(ensemble_model)
kable(variableImportance,digits=2,caption="Variable Importance") %>%
  kable_styling(bootstrap_options = c("striped"))
write.csv(variableImportance,file = paste0(ResponseCurvesOutput, "/varImp.csv"))


### Generate and export response curves in order of variable importance
topPreds <- variableImportance[with(variableImportance,order(-overall)),]
varNames<-rownames(topPreds)
## combine predictions from each model for each variable
## train data needs to be the training data used in the individual models used to build the ensemble model. This info can be extracted from the best ensemble model (ie. bestModel)
model.train<-ensemble_model$models[[1]]$trainingData





gbm.partial.list<-lapply(varNames,partial_gbm)
#glm.partial.list<-lapply(varNames,partial_glm)
rf.partial.list<-lapply(varNames,partial_rf)
mars.partial.list<-lapply(varNames,partial_mars)


#names(glm.partial.list)<-varNames
names(gbm.partial.list)<-varNames
names(rf.partial.list)<-varNames
names(mars.partial.list)<-varNames

#glm.partial.df<-as.data.frame(glm.partial.list)
gbm.partial.df<-as.data.frame(gbm.partial.list)
rf.partial.df<-as.data.frame(rf.partial.list)
mars.partial.df<-as.data.frame(mars.partial.list)

predx<-data.frame()
predy<-data.frame()

for (i in varNames){
  predx <- rbind(predx, as.data.frame(paste(i,i,sep=".")))
  predy<- rbind(predy,as.data.frame(paste(i,"yhat",sep=".")))
}
names(predx)<-""
names(predy)<-""

predx1<-t(predx)
predy1<-t(predy)


#glm.partial.df$data<-'GLM'
gbm.partial.df$data<-'GBM'
rf.partial.df$data<-'RF'
mars.partial.df$data<-'MARS'

all_dfs<-rbind.data.frame(gbm.partial.df,rf.partial.df,mars.partial.df)
allplots<-purrr::map2(predx1,predy1, ~responseCurves(.x,.y))

#export plots as PNGs
for(i in seq_along(allplots)){
  png(paste0(ResponseCurvesOutput,"/",i,".png"),width = 5, height = 5, units = "in",res=300)
  print(allplots[[i]])
  dev.off()
}



### Plot response curves

par(mfrow=c(3,4))
for(i in seq_along(allplots)){
  print(allplots[[i]])
}


