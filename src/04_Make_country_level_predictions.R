#--------------------------------------------
#-----------To do: specify project-----------
#--------------------------------------------
#specify project name
projectname<-"Buffered_occurrences_FINAL"


#--------------------------------------------
#-----------  Load packages  ----------------
#--------------------------------------------
packages <- c("viridis", "dplyr", "grid", "here", "qs","terra", "sf", "ggplot2","RColorBrewer","magick","patchwork",
              "ape", "geoR", "raster", "pdp", "purrr", "sp", "ggspatial", "tidyterra"
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


#--------------------------------------------
#----------  To do: specify country  --------
#--------------------------------------------
#If you'd like to predict for another country, change the shapefile
country_name<-"Belgium"
country<-sf::st_read(here("./data/external/GIS/Belgium/belgium_boundary.shp"))
country_ext<-terra::ext(country) 
country_vector <- terra::vect(country) #Convert to a SpatVector, used for masking


#--------------------------------------------
#----------- Load taxa info  ----------------
#--------------------------------------------
taxa_info<-read.csv2(paste0(project_path,"/taxa_info.csv"))
accepted_taxonkeys<-taxa_info%>%
  pull(speciesKey)%>%
  unique()



#-------------------------------------------------
#---------- Load habitat raster data -------------
#-------------------------------------------------
habitat<-list.files((here("./data/external/habitat")),pattern='tif',full.names = T)
habitat_stack<-rast(habitat[c(1:5,7)]) #Distance to water (layer 6) has another extent


#--------------------------------------------
#--------Source helper functions-------------
#--------------------------------------------
source("./src/helper_functions.R")


#--------------------------------------------
#-----------  Start loop   ----------------
#--------------------------------------------
key<-1311477

#Extract species name
species<-taxa_info%>%
  filter(accepted_taxonkeys==key)%>%
  pull(acceptedScientificName)%>%
  unique()

#Define taxonkey
taxonkey<- key

#Read in globalmodels object that was stored as part of  script 03_fit_European_model
eumodel<-qread( paste0("./data/projects/",projectname,"/EuropeanModelOutput/EuropeanModelOutput.qs"))

#Read in different data objects stored in globalmodels
euocc<-eumodel$occurrences #occurrences in point geometry
bestModel<-unwrap(eumodel$ensemble_model) #global_ensemble_model
fullstack_be<-unwrap(eumodel$fullstack_be)


### Subset Belgium occurrences 
#occ.eu is in WGS84, convert to same projection as country level shapefile (which is the same proj used for model outputs)
suppressWarnings(
  occ.country <- euocc%>%
    st_transform(crs=st_crs(country))%>%
    st_intersection(country)%>%
    dplyr::select(geometry)%>%
    mutate(decimalLongitude = sf::st_coordinates(.)[,1],
           decimalLatitude = sf::st_coordinates(.)[,2])
)

#--------------------------------------------
#-------- Plot country occurrences ---------
#--------------------------------------------
ggplot()+ 
  geom_sf(data = country,  colour = "black", fill = NA)+
  geom_point(data=occ.country, aes(x=decimalLongitude, y= decimalLatitude),  fill="green", shape = 22, colour = "black", size=3)+
  labs(x="Longitude", y="Latitude")+
  theme_bw()


#--------------------------------------------
#-Create country predictions using best model -
#--------------------------------------------
# creates  country level rasters using the European level models
system.time({
  ens_pred_hab_be<-terra::predict(fullstack_be,bestModel,type="prob", na.rm=TRUE)
})


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
  geom_sf(data = occ.country, color = "black", fill = "red", 
          size = 1.5, shape = 21) +
  theme_bw() +
  labs(fill = "Suitability")
country_plot

#Create an empty plot to fill PDF
empty_plot <- ggplot() + 
  theme_void() + 
  theme(plot.background = element_blank()) 

#Create final plot
plot_final<-country_plot /empty_plot 
plot_final

#-------------------------------------------------
#- Export country predictions as raster and PDF --
#-------------------------------------------------
#---------Specify folder paths------------
CountryPredictions <- file.path("./data/projects", projectname, "CountryPredictions")
if (!dir.exists(CountryPredictions)) {
  dir.create(CountryPredictions, recursive = TRUE)
}
#---------------Export raster-------------
writeRaster(ens_pred_hab_be,
            filename=file.path(CountryPredictions,paste("hist_",country_name,".tif",sep="")),
            overwrite=TRUE)

#---------------Export PDF----------------
#Define the file paths
plot_png_path <- file.path(CountryPredictions,paste("hist_",country_name,".png",sep=""))
plot_pdf_path <- file.path(CountryPredictions,paste("hist_",country_name,".pdf",sep=""))

# Save each plot as a PDF file
ggsave(filename = paste0("hist_",country_name,".png"), plot = plot_final, 
       device = "png", width =8.27 , height = 11.69, path= CountryPredictions)

# Read the PNG image back in
img <- image_read(plot_png_path)

# Start a PDF device for output
pdf(plot_pdf_path, width = 8.27, height = 11.69)

# Create a layout for title and image
grid.newpage()

# Add title at the top of the PDF
grid.text(
  label = bquote(italic("Vespa velutina")),
  x = 0.5, y = 0.95, just = "center", gp = gpar(fontsize = 12, fontface = "bold")
)

# Add the PNG image below the title
grid.raster(img, width = unit(0.9, "npc"), height = unit(0.9, "npc"), y = 0.47)

# Close the PDF device
while (dev.cur() > 1) dev.off()

# Remove the PNG file from the local directory
file.remove(plot_png_path)


#-------------------------------------------------
#- Clip habitat raster stack to extent of country --
#-------------------------------------------------
habitat_only_stack<-terra::crop(habitat_stack,country)
habitat_only_stack_be<-terra::mask(habitat_only_stack,country)


#-----------------------------------------------------------
#-Create individual RCP climate raster stacks for country --
#-----------------------------------------------------------
be26 <- list.files((here("./data/external/climate/byEEA_finalRCP/belgium_rcps/rcp26")),pattern='tif',full.names = T)
belgium_stack26 <- rast(be26)

be45 <- list.files((here("./data/external/climate/byEEA_finalRCP/belgium_rcps/rcp45")),pattern='tif',full.names = T)
belgium_stack45 <- rast(be45)

be85 <- list.files((here("./data/external/climate/byEEA_finalRCP/belgium_rcps/rcp85")),pattern='tif',full.names = T)
belgium_stack85 <- rast(be85)


#--------------------------------------------------------------------
#-Combine habitat stacks with climate stacks for each RCP scenario --
#--------------------------------------------------------------------
fullstack26_list <- list(belgium_stack26,habitat_only_stack_be)
fullstack26 <- rast(fullstack26_list) 

fullstack45_list <- list(belgium_stack45,habitat_only_stack_be)
fullstack45 <- rast(fullstack45_list) 

fullstack85_list <- list(belgium_stack85,habitat_only_stack_be)
fullstack85 <- rast(fullstack85_list) 

country_layers<-list(
  "historical"=list("layers"=fullstack_be,
                    "scenario"= "hist",
                    "scenario_title"="historical"),
  "rcp26"=list("layers"=fullstack26,
               "scenario"="rcp26",
               "scenario_title"="RCP 2.6"), 
  "rcp45"=list("layers"=fullstack45,
               "scenario"="rcp45",
               "scenario_title"="RCP 4.5"),
  "rcp85"=list("layers"=fullstack85,
               "scenario"="rcp85",
               "scenario_title"="RCP 8.5")
)

### Create and export RCP risk maps for each RCP scenario
ens_pred_hist <- raster::predict(fullstack_be, bestModel, type = "prob", na.rm = TRUE)
ens_pred_hab26<-raster::predict(fullstack26,bestModel,type="prob", na.rm=TRUE)
writeRaster(ens_pred_hab26, filename=file.path(CountryPredictions,paste("be_",taxonkey, "_rcp26.tif",sep="")), overwrite=TRUE) 
exportPDF(ens_pred_hab26,taxonkey,"Vespa velutina", "", "rcp26.pdf")
ens_pred_hab45<-raster::predict(fullstack45,bestModel,type="prob", na.rm=TRUE)
writeRaster(ens_pred_hab45, filename=file.path(CountryPredictions,paste("be_",taxonkey, "_rcp45.tif",sep="")), overwrite=TRUE) 
exportPDF(ens_pred_hab45,taxonkey,"Vespa velutina", "", "rcp45.pdf")
ens_pred_hab85<-raster::predict(fullstack85,bestModel,type="prob", na.rm=TRUE)
writeRaster(ens_pred_hab85, filename=file.path(CountryPredictions,paste("be_",taxonkey, "_rcp85.tif",sep="")), overwrite=TRUE) 
exportPDF(ens_pred_hab85,taxonkey,"Vespa velutina", "", "rcp85.pdf")



### Create and export RCP risk maps for each RCP scenario

par(mfrow=c(2,2), mar= c(2,3,0.8,0.8))
plot(ens_pred_hist,breaks=brks, lab.breaks=brks)
plot(ens_pred_hab26,breaks=brks, lab.breaks=brks)
plot(ens_pred_hab45,breaks=brks, lab.breaks=brks)
plot(ens_pred_hab85,breaks=brks, lab.breaks=brks)



### Create and export "difference maps": the difference between predicted risk by each RCP scenario and historical climate
hist26_diff_hab <- ens_pred_hab26 - ens_pred_hist
writeRaster(hist26_diff_hab,filename=file.path(CountryPredictions,paste("be_",taxonkey, "_rcp26_diff.tif",sep="")) , overwrite=TRUE) 
exportPDF(hist26_diff_hab,taxonkey,"Vespa velutina","rcp26_diff.pdf","TRUE")


hist45_diff_hab<-ens_pred_hab45 - ens_pred_hist
writeRaster(hist45_diff_hab,filename=file.path(CountryPredictions,paste("be_",taxonkey, "_rcp45_diff.tif",sep="")),overwrite=TRUE) 
exportPDF(hist45_diff_hab,taxonkey,"Vespa velutina","rcp45_diff.pdf","TRUE")


hist85_diff_hab<-ens_pred_hab85 - ens_pred_hist
writeRaster(hist85_diff_hab, filename=file.path(CountryPredictions,paste("be_",taxonkey, "_rcp_85_diff.tif",sep="")), overwrite=TRUE) 
exportPDF(hist85_diff_hab,taxonkey, "Vespa velutina","rcp85_diff.pdf","TRUE")

par(mfrow=c(2,2), mar= c(2,3,0.8,0.8))
plot(hist26_diff_hab)
plot(hist45_diff_hab)
plot(hist85_diff_hab)




### Check spatial autocorrelation of residuals to assess whether occurrence data should be thinned
#### derive residuals from best model
predEns1<-bestModel$ens_model$pred
obs.numeric<-ifelse(predEns1$obs == "absent",0,1)


#### standardize residuals
hab.res<-stdres(obs.numeric,predEns1$present)

# specify corresponding model number from eu_presabs.coord datafile to join data with xy locations. If best model is "X1", join with eu_presabs.coord$X1
res.best.coords1<-cbind(coordinates(eu_presabs.coord$X1),occ.full.data.factor$X1)
removedNAs.coords<-na.omit(res.best.coords1)
res.best.coords<-cbind(removedNAs.coords,hab.res)
res.best.geo<-as.geodata(res.best.coords,coords.col=1:2,data.col = 3)
summary(res.best.geo) #note distance is in meters


### Check Morans I.

#If Moran's I is very low (<0.10), or not significant, do not need to thin occurrences.
res.best.df<-as.data.frame(res.best.coords)
occ.dists <- as.matrix(dist(cbind(res.best.df[1], res.best.df[2])))
occ.dists.inv <- 1/occ.dists
diag(occ.dists.inv) <- 0
Moran.I(res.best.df$hab.res,occ.dists.inv,scaled=TRUE,alternative="greater")




### Quantify confidence of predicted values using class conformal prediction


# quantify confidence for country level predictions based on historical climate and under RCP scenarios of climate change

set.seed(1609)
pvalsdf_hist<-classConformalPrediction(bestModel,ens_pred_hist)
set.seed(447)
pvalsdf_rcp26<-classConformalPrediction(bestModel,ens_pred_hab26)
set.seed(568)
pvalsdf_rcp45<-classConformalPrediction(bestModel,ens_pred_hab45)
set.seed(988)
pvalsdf_rcp85<-classConformalPrediction(bestModel,ens_pred_hab85)

# option to export confidence and pvals as csv 
# write.csv(pvalsdf_hist,file=paste(genOutput,"confidence_",taxonkey, "_hist.csv",sep=""))


### Create confidence maps
brks <- seq(0, 1, by=0.1) 
nb <- length(brks)-1 
pal <- colorRampPalette(rev(brewer.pal(4, 'Spectral')))
cols<-pal(nb)


par(mfrow=c(2,2), mar= c(2,3,0.8,0.8))
hist.conf.map<-confidenceMaps(pvalsdf_hist,taxonkey,"Vespa velutina",maptype="hist_conf")
rcp26.conf.map<-confidenceMaps(pvalsdf_rcp26,taxonkey,"Vespa velutina",maptype="rcp26_conf")
rcp45.conf.map<-confidenceMaps(pvalsdf_rcp45,taxonkey,"Vespa velutina",maptype="rcp45_conf")
rcp85.conf.map<-confidenceMaps(pvalsdf_rcp85,taxonkey,"Vespa velutina",maptype="rcp85_conf")



### Mask areas of below a set confidence level  

# Cutoff for "high" confidence can be modified below. Cutoff should be a value between 0 and 1. Values that are less than the cutoff are shown in gray.
cutoff<-0.70

conf.brks <- seq(0,1, by=0.1) 
nb <- length(conf.brks) 
pal <- colorRampPalette(rev(brewer.pal(4, 'Spectral')))
cols<-pal(nb)

par(mfrow=c(2,2), mar= c(2,3,0.9,0.8))
m1<-hist.conf.map < cutoff
m1_spat <- rast(m1)
hist_masked <- mask(ens_pred_hist, m1_spat, maskvalue = TRUE)
plot(hist_masked,breaks=conf.brks, col=cols,lab.breaks=conf.brks)

m2<-rcp26.conf.map < cutoff
m2_spat<-rast(m2)
rcp26_masked<-mask(ens_pred_hab26,m2_spat,maskvalue=TRUE)
plot(rcp26_masked,breaks=conf.brks, col=cols,lab.breaks=conf.brks)

m3<-rcp45.conf.map < cutoff
m3_spat<-rast(m3)
rcp45_masked<-mask(ens_pred_hab45,m3_spat,maskvalue=TRUE)
plot(rcp45_masked,breaks=conf.brks, col=cols,lab.breaks=conf.brks)

m4<-rcp85.conf.map < cutoff
m4_spat<-rast(m4)
rcp85_masked<-mask(ens_pred_hab85,m4_spat,maskvalue=TRUE)
plot(rcp85_masked,breaks=conf.brks, col=cols,lab.breaks=conf.brks)

### confidence map of best model at EU level
brks <- seq(0, 1, by=0.1) 
nb <- length(brks)-1 
pal <- colorRampPalette(rev(brewer.pal(4, 'Spectral')))
set.seed(792)  
pvalsdf_hist_eu<-classConformalPrediction(bestModel,ens_pred_hist)
hist.conf.map.eu<-confidenceMaps(pvalsdf_hist_eu,taxonkey,"Vespa velutina",maptype="hist_conf_eu")


#Here responce curves moved to script 5

###  Evaluate the performance of each the EU level ensemble models using independent data set from the future 
#####################################################################

# read in and prepare independent data
#2011-2021
eval.data<-read.csv("C:/Users/amyjs/Documents/projects/xps15/xps15/wiSDM/data/external/0001753-230828120925497/0001753-230828120925497.csv",header=TRUE,sep ="\t",quote="")

#enter value for max coordinate uncertainty in meters.

eval.data.occ<-eval.data %>%
  filter(is.na(coordinateUncertaintyInMeters)| coordinateUncertaintyInMeters< 1000) 

eval.data.occ$lon_dplaces<-sapply(na.omit(eval.data.occ$decimalLongitude), function(x) decimalplaces(x))
eval.data.occ$lat_dplaces<-sapply(eval.data.occ$decimalLatitude, function(x) decimalplaces(x))
eval.data.occ[eval.data.occ$lon_dplaces < 4& eval.data.occ$lat_dplaces < 4 , ]<-NA
eval.data.occ<-eval.data.occ[ which(!is.na(eval.data.occ$lon_dplaces)),]
eval.data.occ<-within(eval.data.occ,rm("lon_dplaces","lat_dplaces"))

eval.data.occ<-eval.data.occ[c("decimalLongitude", "decimalLatitude")]
coordinates(eval.data.occ)<- c("decimalLongitude", "decimalLatitude")
proj4string(eval.data.occ)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")#specify here the existing coord.sys of the data
eval.data.occ.proj<-spTransform(eval.data.occ,rmiproj)

########################################################################
#Convert predicted probabilities of EU level risk maps into binary risk maps (present/absence) using thresholds from earlier step

# Eu level
binary_eu_rasters<-sapply(names(thresholds), function(x) raster::reclassify(ens_pred_hab_eu1[[x]],c(0,thresholds[[x]]$predicted,0, thresholds[[x]]$predicted,1,1)),simplify=FALSE)
testeval.eu.bin.rast<-sapply(names(binary_eu_rasters), function(x) eu_eval(binary_eu_rasters[[x]],eval.data.occ.proj),simplify=FALSE)
testeval.eu.bin.rast

binary_be_rasters<-sapply(names(thresholds), function(x) raster::reclassify(ens_pred_hab_be[[x]],c(0,thresholds[[x]]$predicted,0, thresholds[[x]]$predicted,1,1)),simplify=FALSE)
testeval.be.bin.rast<-sapply(names(binary_eu_rasters), function(x) eu_eval(binary_be_rasters[[x]],eval.data.occ.proj),simplify=FALSE)
testeval.be.bin.rast