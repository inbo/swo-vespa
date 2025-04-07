#--------------------------------------------
#------To do: specify project & dataset------
#--------------------------------------------
#specify project name
projectname<-"Buffered_occurrences_FINAL"


#--------------------------------------------
#-----------    Load packages      ----------
#--------------------------------------------
options("rgdal_show_exportToProj4_warnings"="none")

packages <- c( "dplyr", "here", "qs","terra", "sf", "ggplot2","RColorBrewer","magick","patchwork","grid", "tidyterra", "viridisLite",
               "sp", "raster", "dismo", "caret", "kableExtra", "earth", "Formula", "plotmo", "plotrix", "mapview"
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
#--------- Create output folder -------------
#--------------------------------------------
project_path <- file.path("./data/projects",projectname)
GlobalModelOutput<-paste0(project_path,"/GlobalModelOutput")
EuropeanModelOutput<-paste0(project_path,"/EuropeanModelOutput")

if (!dir.exists(EuropeanModelOutput)) {
  dir.create(EuropeanModelOutput, recursive = TRUE)
}

#--------------------------------------------
#----------  To do: specify country  --------
#--------------------------------------------
#If you'd like to predict for another country, change the shapefile
country_name<-"Belgium"
country<-sf::st_read(here("./data/external/GIS/Belgium/belgium_boundary.shp"))
country_ext<-terra::ext(country) 
country_vector <- terra::vect(country) #Convert to a SpatVector, used for masking


#--------------------------------------------
#------- Source helper functions     --------
#--------------------------------------------
source("./src/helper_functions.R")


#--------------------------------------------
#---------   Load shape of Europe   ---------
#--------------------------------------------
euboundary<-st_read(here("./data/external/GIS/Europe/EUROPE.shp")) 

#--------------------------------------------
#--------Load global climate rasters --------
#--------------------------------------------
globalclimrasters <- list.files((here("./data/external/climate/trias_CHELSA")),pattern='tif',full.names = T) #import CHELSA data
globalclimpreds_terra <- terra::rast(globalclimrasters)


#--------------------------------------------
#-------- Load European habitat rasters -----
#--------------------------------------------
habitat_files<-list.files((here("./data/external/habitat")),pattern='tif',full.names = T)
habitat_stack<-rast(habitat_files[c(1:5,7)]) #Distance to water (layer 6) has another extent


#--------------------------------------------
#-------- Load European climate rasters -----
#--------------------------------------------
rmiclimrasters <- list.files((here("./data/external/climate/rmi_corrected")),pattern='tif',full.names = T) 
rmiclimrasters 
rmiclimpreds<-rast(rmiclimrasters) 


#---------------------------------------------
#----- Remove NA pixels from predictors ------
#---------------------------------------------
#This is to avoid that some layers have NA while others have values in certain pixels
#First mask pixels in the rasterstack where at least one layer has NA
na_mask_rmiclimpreds <- app(rmiclimpreds, function(x) any(is.na(x)))
na_mask_habitat_stack<- app(habitat_stack, function(x) any(is.na(x)))
rmiclimpreds<- mask(rmiclimpreds, na_mask_rmiclimpreds, maskvalue=1)
habitat_stack<- mask(habitat_stack, na_mask_habitat_stack, maskvalue=1)

#Second mask rmiclimpreds with habitat_stack and vice versa
rmiclimpreds<-mask(rmiclimpreds, habitat_stack[[1]])
habitat_stack<-mask(habitat_stack, rmiclimpreds[[1]])



#--------------------------------------------
#------------- Load species data -----------
#--------------------------------------------
taxa_info<-read.csv2(paste0(project_path,"/taxa_info.csv"))

accepted_taxonkeys<-taxa_info%>%
  pull(speciesKey)%>%
  unique()


#--------------------------------------------
#----------- Start SDM modelling   ----------
#--------------------------------------------
key<-1311477

#--------------------------------------------
#--------Extract species-specific data  -----
#--------------------------------------------
#Extract species name
species<-taxa_info%>%
  filter(speciesKey==key)%>%
  pull(acceptedScientificName)%>%
  unique()

#Define taxonkey
taxonkey<- key

#Read in globalmodels object that was stored as part of  script 02_fit_global_model
globalmodels<-qread( paste0(GlobalModelOutput,"/GlobalModelOutput.qs"))

#Extract different data objects stored in globalmodels

global.occ<-qread(paste0(project_path,"/thinned_global_occurrences.qs"))
biasgrid_sub<-terra::unwrap(globalmodels$biasgrid)
global_model<-terra::unwrap(globalmodels$model_predictions)
model_accuracy<-globalmodels$model_accuracy


#--------------------------------------------
#-------Prepare occurrence dataset-----------
#--------------------------------------------
global.occ.LL.cleaned<-global.occ %>%
  dplyr::select(c(decimalLongitude,decimalLatitude))


#--------------------------------------------
#------ Remove duplicates per grid cell -----
#--------------------------------------------
global.occ.LL.cleaned$cell<-terra::cellFromXY(globalclimpreds_terra , global.occ.LL.cleaned) #Indicate for each occurrence point in which cell of the raster it falls
global.occ.LL.cleaned<-global.occ.LL.cleaned[!is.na(global.occ.LL.cleaned$cell),]
unique_occurrences <- !duplicated(global.occ.LL.cleaned$cell)# Identify unique occurrences
global.occ.LL.cleaned <- global.occ.LL.cleaned[unique_occurrences, 1:2] # Subset the occurrence points to keep only one occurrence per raster cell 

global.occ.LL.cleaned<- terra::extract(globalclimpreds_terra, global.occ.LL.cleaned, xy = T, ID=F)%>%
  dplyr::filter(rowSums(is.na(.[, 1:(ncol(.) - 2)])) == 0)%>% #Keep rows that do not have any NA values in column 1- 3rd last 
  dplyr::select(c(x,y))%>%
  dplyr::rename(decimalLongitude=x,
                decimalLatitude=y) #Extract climatic values of occurrence points from each raster layer and remove occurrence points that fall in cells with NA values in at least one rasterlayer

#Convert to sf dataframe
global.occ.LL.cleaned$species<- rep(1,length(global.occ.LL.cleaned$decimalLongitude)) #adds columns indicating species presence (1) needed for modeling
global.occ.sf<-st_as_sf(global.occ.LL.cleaned, coords=c("decimalLongitude", "decimalLatitude"), crs=4326, remove= FALSE)



#-----------------------------------------------
#----- Create subset of European records -------
#-----------------------------------------------
#Check for occurrences that fall within Europe
eu_occ <- global.occ.sf[st_intersects(global.occ.sf, euboundary, sparse = FALSE), ] %>%
  dplyr::select(decimalLatitude, decimalLongitude, species) %>%
  dplyr::filter(!is.na(decimalLatitude))%>%
  sf::st_transform(crs=st_crs(rmiclimpreds))

# Convert to crs of rmiclimpreds
eu_occ<-eu_occ%>%
  st_coordinates()%>%
  cbind(., eu_occ)%>%
  select(-c(decimalLatitude, decimalLongitude))

#Only keep occurrences in pixels that have predictor data (not NA's)
extracted_value <- terra::extract(rmiclimpreds[[1]], vect(eu_occ))
eu_occ$extracted_value<-extracted_value[,2]
eu_occ <- eu_occ[!is.na(eu_occ$extracted_value), ]

# Keep XY coordinates
euocc<-eu_occ%>%
  st_coordinates()


#--------------------------------------------
#----- Clip biasgrid to European extent -----
#--------------------------------------------
ecoregions_eu<-terra::crop(biasgrid_sub, euboundary)
biasgrid_eu <- project(ecoregions_eu, rmiclimpreds) #reproject the ecoregions raster to match the spatial properties of rmi climpreds


#--------------------------------------------
#----------- Visualize biasgrid  ------------
#--------------------------------------------
biasgrid_map<-ggplot()+ 
  geom_sf(data = euboundary,  colour = "black", fill = NA)+
  geom_spatraster(data=biasgrid_eu)+
  scale_fill_continuous(na.value = "transparent",low = "blue", high = "orange")+
  labs(x="Longitude", y="Latitude")+
  theme_bw()
biasgrid_map
ggsave(filename = "biasgrid_map.png", plot = biasgrid_map, 
       device = "png", width =16 , height = 20, path= EuropeanModelOutput)

#--------------------------------------------
#Mask areas of high suitability in global model
#--------------------------------------------
#Create a mask of the global_model rasterlayer, containing only areas that are predicted to contain occurrences
#CUTOFFm<- global_model >= model_accuracy$threshold

#Mask the global_model layer with this occurrence layer (i.e., convert all areas where occurrences are predicted to NA)
#CUTOFFglobal_mask<-mask(global_model,m,maskvalue=TRUE)
#CUTOFFglobal_masked_proj<-project(global_mask,biasgrid_eu)
global_model_proj<-project(global_model, biasgrid_eu)

#New: mask with one of the environmental layers to make sure no pseudoabsences are generated outside the environmental layers
global_model_proj<-terra::mask(global_model_proj, rmiclimpreds[[1]]) 


#--------------------------------------------
#--Create sampling area for pseudoabsences --
#--------------------------------------------
# Combine areas of low predicted habitat suitability with bias grid to exclude low sampled areas and areas of high suitability
#CUTOFFpseudoSamplingArea<-mask(global_masked_proj,biasgrid_eu)
pseudoSamplingArea<-mask(global_model_proj,biasgrid_eu)

# make sure PA's are not sampled within cells with occurrences
alldata_Europe<-qread( paste0(project_path,"/Europe_occurrences.qs"))
alldata_Europe_sf <- st_as_sf(alldata_Europe, coords = c("decimalLongitude", "decimalLatitude"), crs=4326)%>%
  st_transform(crs=st_crs(rmiclimpreds))
alldata_Europe_raster <- rast(ext(pseudoSamplingArea), resolution = res(pseudoSamplingArea))
alldata_Europe_raster <- rasterize(alldata_Europe_sf, alldata_Europe_raster, field = 1, background = NA)
crs(alldata_Europe_raster) <- crs(rmiclimpreds)
pseudoSamplingArea <- mask(pseudoSamplingArea, alldata_Europe_raster, maskvalue=1)

pseudoSamplingArea_reverse<- 1-pseudoSamplingArea

#--------------------------------------------
#-------- Plot pseudosampling area ----------
#--------------------------------------------
pseudoSamplingArea_proj<-project(pseudoSamplingArea, global_model)
pseudosamplingarea_plot <- global_model
values(pseudosamplingarea_plot) <- NA
values(pseudosamplingarea_plot)[!is.na(values(global_model))] <- TRUE
values(pseudosamplingarea_plot)[!is.na(values(global_model)) & !is.na(values(pseudoSamplingArea_proj))] <- FALSE

png(paste0(EuropeanModelOutput,"/pseudoSamplingArea.png"), width = 800, height = 600)  
plot(pseudosamplingarea_plot)
dev.off() 
plot(pseudosamplingarea_plot)

#--------------------------------------------
#Randomly generate pseudoabsences in pseudoSamplingArea
#--------------------------------------------
# set number of pseudoabsences equal to the number of presences
numb.eu.pseudoabs<-nrow(euocc)

# Generate pseudoabsences 10 times, store in a list with 10 datasets and names them X1-X10
setlist<-seq(1,10,1)
set.seed(120)
#CUTOFFpseudoabs_pts <- lapply(setlist, generate_pseudoabs, mask = raster(pseudoSamplingArea), alternative_mask = raster(global_masked_proj), n = numb.eu.pseudoabs, p = euocc)
pseudoabs_pts <- lapply(setlist, generate_pseudoabs, mask = raster(pseudoSamplingArea_reverse), alternative_mask = raster(global_model_proj), n = numb.eu.pseudoabs, p = euocc, weighted=TRUE)

names(pseudoabs_pts) <- paste0("X", setlist)


#--------------------------------------------
#Prepare presence-absence dataset for modelling
#--------------------------------------------
# extract data from environmental predictors for absences
pseudoabs_pts1<-lapply(pseudoabs_pts, function(x) terra::extract(rmiclimpreds,x, ID=FALSE))

# add occ column with value 0 (indicating absences)
pseudoabs_pts2<-lapply(pseudoabs_pts1, function(x) add.occ(x,0))

# extract environmental data for eu presences and add presence indicator (1)
presence<-as.data.frame(euocc)
names(presence)<- c("x","y")
presence1<-terra::extract(rmiclimpreds,presence, ID=FALSE)
occ<-rep(1,nrow(presence1))
presence1<-cbind(presence1,occ)

# join each pseudoabsence set with presences 
eu_presabs.pts<-lapply(pseudoabs_pts2,  function(x) rbind(x, presence1))
eu_presabs.coord<-lapply(pseudoabs_pts, function(x) rbind(x,presence))



#--------------------------------------------
#--Visualize presence-pseudoabsence dataset--
#--------------------------------------------
pseudoabs_sf <- st_as_sf(pseudoabs_pts$X1, coords = c("x", "y"), crs=st_crs(rmiclimpreds))
euocc_sf <- st_as_sf(eu_occ, coords = c("decimalLongitude", "decimalLatitude"), crs=4326 )%>%
  select('geometry', 'species')
pseudoabs_sf$species<-0
eu_presabs_sf<-rbind(pseudoabs_sf, euocc_sf)
m<-mapview(pseudoSamplingArea_reverse, 
           col.regions = colorRampPalette(c("blue", "orange")),
           alpha=1, 
           na.color = "transparent", 
           layer.name = "Pseudosampling area") +
  mapview(eu_presabs_sf, zcol = "species", 
          col.regions = c("red", "yellow"),
          layer.name = "Species distribution")

m
mapshot(m, url = file.path(EuropeanModelOutput, "presence-pseudoabsence_map.html"))



#--------------------------------------------
#--Remove highly correlated predictors from training data --
#--------------------------------------------
# convert eu data to dataframe
eu_presabs.pts.df<-lapply(eu_presabs.pts,function(x) as.data.frame(x))

# find attributes that are highly corrected
highlyCorrelated_climate <-lapply(eu_presabs.pts.df, function(df) as.data.frame(cor(df[, 1:13], use = "complete.obs")))

#Calculate the mean correlation over the 10 datsets and identify highly correlated variables
mean_correlation_matrix <- Reduce("+", highlyCorrelated_climate) / length(highlyCorrelated_climate)
drop_climate<-findCorrelation(as.matrix(mean_correlation_matrix), cutoff=0.7,exact=TRUE,names=TRUE)

#Only keep layers that are not highly correlated
rmiclimpreds_uncor <- subset(rmiclimpreds, !(names(rmiclimpreds) %in% drop_climate))


#--------------------------------------------
#- Add habitat and anthropogenic predictors -
#--------------------------------------------
#combine uncorrelated climate variable selected earlier with habitat layers
fullstack<-c(rmiclimpreds_uncor,habitat_stack) 

# create a rasterstack for specified country (Belgium in example case)
fullstack_crop<-crop(fullstack,country_ext)
fullstack_be<-mask(fullstack_crop,country_vector)


#-----------------------------------------------------------
#- Extract predictor values for presences and pseudoabsences
#-----------------------------------------------------------
# Why first extract for environmental predictors and not immediately do this step?
occ.full.data <-lapply(eu_presabs.coord, function(x) extract(fullstack,x, ID=FALSE))


#--------------------------------------------
#--- Remove highly correlated predictors ----
#--------------------------------------------
highlyCorrelated_full <-lapply(names(occ.full.data),function(x)
  findCorrelation(cor(occ.full.data[[x]],use = 'complete.obs'), cutoff=0.7,exact=TRUE,names=TRUE))

highlyCorrelated_vec<-unlist(highlyCorrelated_full)
eupreds1<-as.data.frame(highlyCorrelated_vec)
kable(eupreds1) %>%
  kable_styling(bootstrap_options = c("striped"))

# Remove highly correlated predictors from dataset holding occurrences
occ.full.data<-sapply(names(occ.full.data),function (x) occ.full.data[[x]][,!(colnames(occ.full.data[[x]]) %in% highlyCorrelated_vec)],simplify=FALSE)

# Remove highly correlated predictors from rasterlayers
keep_layers <- !(names(fullstack) %in% highlyCorrelated_vec)
fullstack <- subset(fullstack, keep_layers)


#--------------------------------------------
#--- Remove near-zero variance predictors ---
#--------------------------------------------
# identify low variance predictors
nzv_preds<-lapply(names(occ.full.data),function(x) caret::nearZeroVar(occ.full.data[[x]],names=TRUE))
nzv_preds.vec<-unique(unlist(nzv_preds))
nzv_preds.vec

# remove near zero variance predictors. They don't contribute to the model.
occ.full.data<-sapply(names(occ.full.data),function (x) occ.full.data[[x]][,!(colnames(occ.full.data[[x]]) %in% nzv_preds.vec)],simplify=FALSE)


#--------------------------------------------
#-------- Prepare data for modelling --------
#--------------------------------------------
#Convert to dataframe
occ.full.data.df<-lapply(occ.full.data, function(x) as.data.frame(x))

#Add column with occurrence data (occ)
occ.full.data.df<- sapply(names(occ.full.data.df), function (x) cbind(occ.full.data.df[[x]],occ=eu_presabs.pts.df[[x]]$occ, deparse.level=0),simplify=FALSE)

#Recode factor levels of column 'occ' to absent (0) and present(1), and set present as the reference level
occ.full.data.factor<-sapply(names(occ.full.data.df), function (x) factorVars(occ.full.data.df[[x]], "occ"),simplify=FALSE)



#--------------------------------------------
#- Run models with climate and habitat data -
#--------------------------------------------
control <- trainControl(method="cv",
                        number=4,
                        savePredictions="final", 
                        preProc=c("center","scale"),
                        classProbs=TRUE)

mylist<-list(
  #glm =caretModelSpec(method = "glm",maxit=100),
  gbm= caretModelSpec(method = "gbm"),
  rf = caretModelSpec(method = "rf", importance = TRUE),
  earth= caretModelSpec(method = "earth"))

# set.seed(167)
eu_models<-sapply(names(occ.full.data.factor), function(x) model_train_habitat <- caretList(
  occ~., 
  data= occ.full.data.factor[[x]],
  trControl=control,
  tuneList=mylist), 
  simplify=FALSE)


#--------------------------------------------
#---- Display model evaluation statistics----
#--------------------------------------------
EU_ModelResults1<-sapply(names(eu_models), function(x) resamples(eu_models[[x]]),simplify=FALSE)
Results.summary<-sapply(names(EU_ModelResults1), function(x) summary(EU_ModelResults1[[x]]),simplify=FALSE)
Results.summary

#show_euModel_correlation
Model.cor<-sapply(names(eu_models), function(x) modelCor(resamples(eu_models[[x]])),simplify=FALSE)
Model.cor


#--------------------------------------------
#---------- Create ensemble model -----------
#--------------------------------------------
set.seed(458)
lm_ens_hab<-sapply(names(eu_models), function (x) caretEnsemble(eu_models[[x]], 
                                                                trControl=trainControl(method="cv", 
                                                                                       number=10,
                                                                                       savePredictions= "final",
                                                                                       classProbs = TRUE)),
                   simplify=FALSE)


#--------------------------------------------
#- Evaluate each ensemble model's performance -
#-------------------------------------------- 
#based on results from CV, 
#identify threshold where sensitivity=specifity
thresholds<-sapply(names(lm_ens_hab), function(x) findThresh(lm_ens_hab[[x]]$ens_model$pred),simplify=FALSE)

#Using thresholds identified for each model in the previous step, assess performance of each model
# accuracy measures
thresholds.df<-sapply(names(thresholds), function(x) accuracyStats(lm_ens_hab[[x]]$ens_model$pred,thresholds[[x]]$predicted),simplify=FALSE)
thresholds.comb<-do.call(rbind,thresholds.df)
kable(thresholds.comb,digits=2)


#--------------------------------------------
#---------Select best ensemble model --------
#--------------------------------------------
# select best model based on highest PCC, and, in case there are multiple rows with the same PCC, the highest AUC
bestmodelname <- thresholds.comb %>%
  filter(PCC == max(PCC)) %>%      
  slice_max(AUC, n = 1)%>%
  rownames()

bestModel<-lm_ens_hab[[bestmodelname]]


#--------------------------------------------
#-Create European predictions using best model -
#--------------------------------------------
system.time({
  ens_pred_hab_eu1<-terra::predict(fullstack,lm_ens_hab[[bestmodelname]],type="prob", na.rm = TRUE)
}) 


#--------------------------------------------
#- Create sf df with occurrences for plotting -
#--------------------------------------------
euocc1<-st_as_sf(as.data.frame(euocc), coords=c("X","Y"),crs=st_crs(rmiclimpreds))


#--------------------------------------------
#-------- ¨Plot predictions for Europe  -----
#--------------------------------------------
#Plot
brks <- seq(0, 1, by=0.1)
nb <- length(brks) - 1
viridis_palette <- viridis(nb)

eu_plot<-ggplot() + 
  geom_spatraster(data = ens_pred_hab_eu1) +
  scale_fill_gradientn(colors = viridis_palette, 
                       breaks = brks, 
                       labels = brks, 
                       na.value = NA) +
  #geom_sf(data = euocc1, color = "black", fill = "red", 
  #size = 1.5, shape = 21) +
  theme_bw() +
  labs(fill = "Suitability")+
  coord_sf(xlim = c(2254476, 6005897), 
           ylim = c(1363659, 5469923))

eu_plot

#Create an empty plot to fill PDF
empty_plot <- ggplot() + 
  theme_void() + 
  theme(plot.background = element_blank()) 

#Create final plot
plot_final<-eu_plot /empty_plot 


#--------------------------------------------
#- Export Eu predictions as raster and PDF --

#---------------Export raster-------------
writeRaster(ens_pred_hab_eu1,
            filename=file.path(paste(EuropeanModelOutput,"/hist_EU.tif",sep="")),
            overwrite=TRUE)

#---------------Export PDF----------------
#Define the file paths
plot_png_path <- file.path(paste(EuropeanModelOutput,"/hist_EU.png",sep=""))
plot_pdf_path <- file.path(paste(EuropeanModelOutput,"/hist_EU.pdf",sep=""))

# Save each plot as a PDF file
ggsave(filename = "hist_EU.png", plot = plot_final, 
       device = "png", width =8.27 , height = 11.69, path= EuropeanModelOutput)

# Read the PNG image back in
img <- image_read(plot_png_path)

# Start a PDF device for output
pdf(plot_pdf_path, width = 8.27, height = 11.69)

# Create a layout for title and image
grid.newpage()

# Add title at the top of the PDF
grid.text(
  label = bquote(italic(.("Vespa velutina")) ~ "( " * .(taxonkey) * ")"),
  x = 0.5, y = 0.95, just = "center", gp = gpar(fontsize = 12, fontface = "bold")
)

# Add the PNG image below the title
grid.raster(img, width = unit(0.9, "npc"), height = unit(0.9, "npc"), y = 0.47)

# Close the PDF device
while (dev.cur() > 1) dev.off()

# Remove the PNG file from the local directory
file.remove(plot_png_path)


#--------------------------------------------
#-Create country predictions using best model -
#--------------------------------------------
# creates  country level rasters using the European level models
system.time({
  ens_pred_hab_be<-terra::predict(fullstack_be,lm_ens_hab[[bestmodelname]],type="prob", na.rm=TRUE)
})


#--------------------------------------------
#- Create sf df with country occurrences for plotting -
#--------------------------------------------
country<-st_transform(country, crs=st_crs(euocc1))
be_occ <- euocc1 %>%
  st_intersection(country)


#--------------------------------------------
#-------- ¨Plot predictions for country -----
#--------------------------------------------
brks <- seq(0, 1, by=0.1)
nb <- length(brks) - 1
viridis_palette <- viridis(nb)

country_plot<-ggplot() + 
  geom_spatraster(data = ens_pred_hab_be) +
  scale_fill_gradientn(colors = viridis_palette, 
                       breaks = brks, 
                       labels = brks, 
                       na.value = NA) +
  #geom_sf(data = be_occ, color = "black", fill = "red", 
  #size = 1.5, shape = 21) +
  theme_bw() +
  labs(fill = "Suitability")
country_plot

#Create an empty plot to fill PDF
empty_plot <- ggplot() + 
  theme_void() + 
  theme(plot.background = element_blank()) 

#Create final plot
plot_final<-country_plot /empty_plot 


#--------------------------------------------
#- Export Eu predictions as raster and PDF --
#-------------------------------------------
#---------------Export raster-------------
writeRaster(ens_pred_hab_be,
            filename=file.path(EuropeanModelOutput,paste("hist_",country_name,".tif",sep="")),
            overwrite=TRUE)

#---------------Export PDF----------------
#Define the file paths
plot_png_path <- file.path(paste(EuropeanModelOutput,"/hist_",country_name,".png",sep=""))
plot_pdf_path <- file.path(paste(EuropeanModelOutput,"/hist_",country_name,".pdf",sep=""))

# Save each plot as a PDF file
ggsave(filename = paste0("/hist_",country_name,".png"), plot = plot_final, 
       device = "png", width =8.27 , height = 11.69, path= EuropeanModelOutput)

# Read the PNG image back in
img <- image_read(plot_png_path)

# Start a PDF device for output
pdf(plot_pdf_path, width = 8.27, height = 11.69)

# Create a layout for title and image
grid.newpage()

# Add title at the top of the PDF
grid.text(
  label = bquote(italic(.("Vespa velutina")) ~ "( " * .(taxonkey) * ")"),
  x = 0.5, y = 0.95, just = "center", gp = gpar(fontsize = 12, fontface = "bold")
)

# Add the PNG image below the title
grid.raster(img, width = unit(0.9, "npc"), height = unit(0.9, "npc"), y = 0.47)

# Close the PDF device
while (dev.cur() > 1) dev.off()

# Remove the PNG file from the local directory
file.remove(plot_png_path)


#--------------------------------------------
#- Save best model, european occurrences, and layers for Belgium -
#--------------------------------------------

eumodel <-list(species = species,
               taxonkey = taxonkey,
               occurrences = euocc1,
               ensemble_model=bestModel,
               fullstack_be=terra::wrap(fullstack_be)
)

qsave(eumodel, paste0(EuropeanModelOutput,"/EuropeanModelOutput.qs"))

print(paste("European model has been created for", species))


#--------------------------------------------
#------------ Open output qs file -----------
#--------------------------------------------
projectname<-"Spatial_thinning10km_nativeonly"
projectname<-"Spatial_thinning50EU10nonEU"
projectname<-"Gradual_pseudoabsences"
project_path <- file.path("./data/projects",projectname)
EuropeanModelOutput<-paste0(project_path,"/EuropeanModelOutput")

eumodels<-qread(paste0(EuropeanModelOutput,"/EuropeanModelOutput.qs"))
eumodels$ensemble_model

#global.ens.thresh<-findThresh(eumodels$ensemble_model$pred)
#ensemble_accurracy<-accuracyStats(global_stack$ens_model$pred,global.ens.thresh$predicted)

rm(list=ls())

