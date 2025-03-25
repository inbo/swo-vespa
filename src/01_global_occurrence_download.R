#-----------------------------------------------------------------------------
#----To do: specify project, species, download key and thinning scenario------
#-----------------------------------------------------------------------------
#specify project name
projectname<-"Buffered_occurrences_without_glm"

# specify the scientific name of the species to be modelled
species<-c("Vespa velutina")

#put this to FALSE if not downloaded yet
#gbif_download_key<-FALSE
gbif_download_key<-"0023851-250127130748423"

#--------------------------------------------
#-----------------Load packages--------------
#--------------------------------------------
packages <- c("rgbif", "dplyr", "purrr", "assertthat", "readr", "here", "qs")

for(package in packages) {
  print(package)
  if( ! package %in% rownames(installed.packages()) ) { install.packages( package ) }
  library(package, character.only = TRUE)
}


#--------------------------------------------
#---------Load helper functions--------------
#--------------------------------------------
source("./src/helper_functions.R")


#--------------------------------------------
#--Create a project folder and a raw folder--
#--------------------------------------------
# Define the folder paths
project_path <- file.path("./data/projects",projectname)
raw_path<-file.path("./data/raw")

# Check and create each folder if necessary
create_folder(project_path, projectname)
create_folder(raw_path, "raw")

#--------------------------------------------
#-----------Retrieve GBIF taxonkeys----------
#--------------------------------------------
# Match species names with the GBIF backbone, retrieve taxon keys from GBIF when a match is found
taxon_df <- as.data.frame(species)

mapped_taxa <- purrr::map_dfr(
  taxon_df$species,
  ~ {
    tryCatch(
      {
        data <- rgbif::name_backbone(name = .x)
        if (length(data) == 0) {
          stop("No match with the GBIF backbone found")
        }
        data
      },
      error = function(e) {
        NULL
      }
    )
  }
)


  

#Make sure that only species info is stored as it is possible that genus information is captured when the species part of the name is not clear
mapped_taxa<-mapped_taxa %>%
  dplyr::filter(rank =="SPECIES")

#Make sure that all species were mapped to the GBIF backbone, if not an error will appear indicating which species are missing
assertthat::assert_that(nrow(mapped_taxa)==length(species),
                        msg=paste0("The following species could not be found in the GBIF backbone taxonomy: "
                                   ,species[!sapply(species, function(x) any(grepl(x,mapped_taxa$scientificName)))])
)

not_accepted <- mapped_taxa %>%
  dplyr::filter(status !="ACCEPTED")

if (nrow(not_accepted)!=0) {
  warning(paste0("The following species do not have an accepted taxonomic status in the GBIF backbone: ",paste(unique(not_accepted$scientificName), collapse=", "),". Their corresponding accepted species names will be used for downloading occurrence data.")
  )
} else {
  paste0("All species are accepted taxa in the GBIF backbone 🎉")
}

#Extract taxonkeys of each species, for synonyms the acceptedUsageKey is stored
accepted_taxonkeys<-mapped_taxa %>%
  dplyr::filter(status =="ACCEPTED")%>%
  pull(usageKey)

if(nrow(not_accepted!=0)){
  synonym_taxonkeys<-mapped_taxa %>%
    dplyr::filter(status !="ACCEPTED")%>%
    pull(acceptedUsageKey)

  accepted_taxonkeys<-c(accepted_taxonkeys, synonym_taxonkeys)
} 

#Keep unique accepted taxonkeys
accepted_taxonkeys<-unique(accepted_taxonkeys)


#--------------------------------------------
#-----------Define download settings---------
#--------------------------------------------
if (gbif_download_key==FALSE){
  
#All basis of record types, except `FOSSIL SPECIMEN` and `LIVING SPECIMEN`, which can have misleading location information (e.g. location of captive animal).
basis_of_record <- c(
  "OBSERVATION", 
  "HUMAN_OBSERVATION",
  "UNKNOWN", 
  "MACHINE_OBSERVATION",
  "OCCURRENCE"
)

#Time period
year_begin <- 1971
year_end <-2025

#Only georeferenced points
hasCoordinate <- TRUE


#--------------------------------------------
#---------------Perform download-------------
#--------------------------------------------
#Note that GBIF credentials are required
gbif_download_key <- occ_download(
  pred_in("taxonKey", accepted_taxonkeys),
  pred_in("basisOfRecord", basis_of_record),
  pred_gte("year", year_begin),
  pred_lte("year", year_end),
  pred("hasCoordinate", hasCoordinate),
  user = rstudioapi::askForPassword("GBIF username"),
  pwd = rstudioapi::askForPassword("GBIF password"),
  email = rstudioapi::askForPassword("Email address for notification")
)

occ_download_wait(gbif_download_key)#Check download status
} else {
  gbif_download_key<-"0023851-250127130748423"
}

#--------------------------------------------
#--------------Retrieve download-------------
#--------------------------------------------

occ_download_get(gbif_download_key, path = here("data","raw"), overwrite=TRUE)
metadata <- occ_download_meta(key = gbif_download_key)
gbif_download_key<-metadata$key

#extract_GBIF_occurrence
raw.path<- here("data", "raw", gbif_download_key)
unzip(paste0(raw.path,".zip"),exdir=raw.path, overwrite=TRUE)
gbif_original<-as.data.frame(data.table::fread(paste0(raw.path,"./occurrence.txt"),header=TRUE))


#--------------------------------------------
#-------- Clean global occurrences----------
#--------------------------------------------
#remove unverified records
identificationVerificationStatus_to_discard <- c( "unverified",
                                                  "unvalidated",
                                                  "not validated",
                                                  "under validation",
                                                  "not able to validate",
                                                  "control could not be conclusive due to insufficient knowledge",
                                                  "1",
                                                  "uncertain",
                                                  "unconfirmed",
                                                  "Douteux",
                                                  "Invalide",
                                                  "Non r\u00E9alisable",
                                                  "verification needed" ,
                                                  "Probable",
                                                  "unconfirmed - not reviewed",
                                                  "validation requested")

#enter value for max coordinate uncertainty in meters, default = 1000
gbif<-gbif_original %>%
  dplyr::filter(speciesKey%in%accepted_taxonkeys) %>%   
  dplyr::filter(is.na(coordinateUncertaintyInMeters)| coordinateUncertaintyInMeters<= 1000) %>%
  dplyr::filter(!str_to_lower(identificationVerificationStatus) %in% identificationVerificationStatus_to_discard)

#Remove coordinates that for both lon and lat values, have less than 4 decimal places
gbif$lon_dplaces<-sapply(gbif$decimalLongitude, function(x) decimalplaces(x))
gbif$lat_dplaces<-sapply(gbif$decimalLatitude, function(x) decimalplaces(x))
gbif[gbif$lon_dplaces < 4 & gbif$lat_dplaces < 4 , ]<-NA
gbif<-gbif[ which(!is.na(gbif$lon_dplaces)),]
gbif<-within(gbif,rm("lon_dplaces","lat_dplaces")) # n= 1758

# Only keep coordinates that are not flagged as potentially problematic
gbif<-gbif%>%
  cc_cen(buffer=100) %>% # remove points within a buffer of 100m around country centroids, default 1km
  cc_cap(buffer=100) %>% # remove capitals centroids (buffer 100m), default 10km
  cc_inst(buffer=100) %>% # remove zoo and herbaria records buffer of 100 m around biodiversity institutes, default 100m
  cc_gbif(buffer=100)%>% #remove around GBIF headquarters in Copenhagen (buffer 100m), default 100m
  cc_zero() #Remove around the 0/0 point (buffer 0.5 degrees)


#Delete occurrences in America, Sweden, Norway, Turkey and one in Spain
gbif <- dplyr::filter(gbif, 
                      !countryCode %in% c("SE", "NO", "JP", "US", "CO"), 
                      stateProvince!="Andalucía",
                      decimalLongitude<20 | decimalLongitude>40)

gbif<-dplyr::select(gbif, c(speciesKey,species, decimalLatitude, decimalLongitude, kingdom, phylum, class))

print(paste((nrow(gbif_original) - nrow(gbif)), "occurrences have been removed.", nrow(gbif), "remain."))


#--------------------------------------------
#------ Load extra data iAsset and manual-----
#--------------------------------------------


manual <- as.data.frame(data.table::fread(here("data", "raw","manual.csv"), header = TRUE))
manual <- manual[!is.na(manual$latitude), ]
iAsset <- as.data.frame(data.table::fread(here("data", "raw","iAsset.csv"), header = TRUE))


#--------------------------------------------
#----------Merge data with gbif--------------
#--------------------------------------------

manual_for_gbif <- data.frame(
  speciesKey = rep(unique(gbif$speciesKey), nrow(manual)),  
  species = rep(unique(gbif$species), nrow(manual)),
  decimalLatitude = manual$latitude,
  decimalLongitude = manual$longitude,
  kingdom = rep(unique(gbif$kingdom), nrow(manual)),  
  phylum = rep(unique(gbif$phylum), nrow(manual)),  
  class = rep(unique(gbif$class), nrow(manual))  
)

iAsset_for_gbif <- data.frame(
  speciesKey = rep(unique(gbif$speciesKey), nrow(iAsset)),  
  species = rep(unique(gbif$species), nrow(iAsset)),
  decimalLatitude = iAsset$latitude,
  decimalLongitude = iAsset$longitude,
  kingdom = rep(unique(gbif$kingdom), nrow(iAsset)),  
  phylum = rep(unique(gbif$phylum), nrow(iAsset)),  
  class = rep(unique(gbif$class), nrow(iAsset))
)

global <- rbind(gbif, manual_for_gbif, iAsset_for_gbif)

#delete foute coordinaat helemaal links boven wereldkaart
global <- dplyr::filter(global, decimalLongitude>-50)

#--------------------------------------------
#-----------Visualize global data------------
#--------------------------------------------

global.occ_sf <- st_as_sf(global, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
global.occ_sf <- global.occ_sf %>%
  mutate(
    decimalLongitude = st_coordinates(.)[, 1],
    decimalLatitude = st_coordinates(.)[, 2])
fig_velutina_global<-mapview(global.occ_sf,
                             layer.name = "Species distribution")
fig_velutina_global
mapshot(fig_velutina_global, url = file.path(project_path, "occurrences_map_global.html"))



#--------------------------------------------
#----------Visualize European data-----------
#--------------------------------------------
Europe <- global%>%
  filter(decimalLongitude<20)

europe.occ_sf <- st_as_sf(Europe, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
europe.occ_sf <- europe.occ_sf %>%
  mutate(
    decimalLongitude = st_coordinates(.)[, 1],
    decimalLatitude = st_coordinates(.)[, 2])
fig_velutina_europe<-mapview(europe.occ_sf,
           layer.name = "Species distribution")
fig_velutina_europe
mapshot(fig_velutina_europe, url = file.path(project_path, "occurrences_map_Europe.html"))

#--------------------------------------------
#--------Visualize non-European data---------
#--------------------------------------------

nonEurope <- global%>%
  filter(decimalLongitude>40)
noneurope.occ_sf <- st_as_sf(nonEurope, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
noneurope.occ_sf <- noneurope.occ_sf %>%
  mutate(
    decimalLongitude = st_coordinates(.)[, 1],
    decimalLatitude = st_coordinates(.)[, 2])

fig_velutina_noneurope<-mapview(noneurope.occ_sf,
                             layer.name = "Species distribution")
fig_velutina_noneurope
mapshot(fig_velutina_noneurope, url = file.path(project_path, "occurrences_map_nonEurope.html"))

print(paste("Before thinning, there are", length(europe.occ_sf$decimalLongitude), "occurrences in invaded range (Europe) and", length(noneurope.occ_sf$decimalLongitude), "in the native range (Asia)"))


#--------------------------------------------
#----------- Spatial thinning----------------
#--------------------------------------------


# Scenario 1: 50 km in invaded, 10 km in native
scenario<-1
thinned_50km_Europe<- thin_points(
  data = Europe, # Dataframe with coordinates
  long_col = "decimalLongitude", # Longitude column name
  lat_col = "decimalLatitude", # Latitude column name
  method = "grid",  # Method for thinning
  thin_dist = 50,  # Thinning distance in km,
  trials = 1, # Number of reps
  all_trials = TRUE # Return all trials
)
thinned_10km_nonEurope<- thin_points(
  data = nonEurope, # Dataframe with coordinates
  long_col = "decimalLongitude", # Longitude column name
  lat_col = "decimalLatitude", # Latitude column name
  method = "grid",  # Method for thinning
  thin_dist = 10,  # Thinning distance in km,
  trials = 1, # Number of reps
  all_trials = TRUE # Return all trials
)

print(paste("After thinning scenario 1, there are", length(thinned_50km_Europe[[1]]$decimalLongitude), "occurrences in invaded range (Europe) and", length(thinned_10km_nonEurope[[1]]$decimalLongitude), "in the native range (Asia)"))

europe50.occ_sf <- st_as_sf(thinned_50km_Europe[[1]], coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
europe50.occ_sf <- europe50.occ_sf %>%
  mutate(
    decimalLongitude = st_coordinates(.)[, 1],
    decimalLatitude = st_coordinates(.)[, 2])
plot_thinned_europe_50 <- ggplot() +
  annotation_map_tile(type = "osm") +
  geom_sf(data = europe50.occ_sf, aes(color = "red"), size = 0.5) +
  theme_minimal() +
  labs(title = "Thinned Locations Europe (50km)", color = "Hornet Data") +
  theme(legend.position = "none")
plot_thinned_europe_50

noneurope10.occ_sf <- st_as_sf(thinned_10km_nonEurope[[1]], coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
noneurope10.occ_sf <- noneurope10.occ_sf %>%
  mutate(
    decimalLongitude = st_coordinates(.)[, 1],
    decimalLatitude = st_coordinates(.)[, 2])
plot_thinned_noneurope_10 <- ggplot() +
  annotation_map_tile(type = "osm") +
  geom_sf(data = noneurope10.occ_sf, aes(color = "red"), size = 0.5) +
  theme_minimal() +
  labs(title = "Thinned Locations non-Europe (10km)", color = "Hornet Data") +
  theme(legend.position = "none")
plot_thinned_noneurope_10


#Scenario 2: Thinning 50km in both areas
#scenario <- 2
thinned_50km_nonEurope<- thin_points(
  data = nonEurope, # Dataframe with coordinates
  long_col = "decimalLongitude", # Longitude column name
  lat_col = "decimalLatitude", # Latitude column name
  method = "grid",  # Method for thinning
  thin_dist = 50,  # Thinning distance in km,
  trials = 1, # Number of reps
  all_trials = TRUE # Return all trials
)

print(paste("After thinning scenario 2, there are", length(thinned_50km_Europe[[1]]$decimalLongitude), "occurrences in invaded range (Europe) and", length(thinned_50km_nonEurope[[1]]$decimalLongitude), "in the native range (Asia)"))

noneurope50.occ_sf <- st_as_sf(thinned_50km_nonEurope[[1]], coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
noneurope50.occ_sf <- noneurope50.occ_sf %>%
  mutate(
    decimalLongitude = st_coordinates(.)[, 1],
    decimalLatitude = st_coordinates(.)[, 2])
plot_thinned_noneurope_50 <- ggplot() +
  annotation_map_tile(type = "osm") +
  geom_sf(data = noneurope50.occ_sf, aes(color = "red"), size = 0.5) +
  theme_minimal() +
  labs(title = "Thinned Locations non-Europe (50km)", color = "Hornet Data") +
  theme(legend.position = "none")
plot_thinned_noneurope_50

if (scenario==1){
  thinned_global<-rbind(thinned_50km_Europe[[1]], thinned_10km_nonEurope[[1]])
  thinned_Europe <-thinned_50km_Europe[[1]]
  thinned_nonEurope <-thinned_10km_nonEurope[[1]]
} else {
  thinned_global<-rbind(thinned_50km_Europe[[1]], thinned_50km_nonEurope[[1]])
  thinned_Europe <-thinned_50km_Europe[[1]]
  thinned_nonEurope <-thinned_50km_nonEurope[[1]]
  
}
  
  


#--------------------------------------------
#------------------Save data-----------------
#--------------------------------------------
#Create dataset taxa_info containing scientific name, canonical name, taxonkeys, gbif download key,...
taxa_info<-data.frame(speciesKey=unique(global$speciesKey),
                      acceptedScientificName=unique(global$species),
                      year_begin=metadata[["request"]][["predicate"]][["predicates"]][[3]][["value"]],
                      year_end=metadata[["request"]][["predicate"]][["predicates"]][[4]][["value"]],
                      gbif_download_key = gbif_download_key,
                      gbif_download_created = format(strptime(metadata$created, "%Y-%m-%dT%H:%M:%S"), "%Y-%m-%d %H:%M:%S"),
                      projectname = projectname)

#Save occurrence data as .qs file and taxa info as .csv
qsave(global, paste0(project_path,"/global_occurrences.qs"))
qsave(Europe, paste0(project_path,"/Europe_occurrences.qs"))
qsave(nonEurope, paste0(project_path,"/nonEurope_occurrences.qs"))

#TO DO: Choose here the right thinned data
qsave(thinned_global, paste0(project_path,"/thinned_global_occurrences.qs"))
qsave(thinned_Europe, paste0(project_path,"/thinned_Europe_occurrences.qs"))
qsave(thinned_nonEurope, paste0(project_path,"/thinned_nonEurope_occurrences.qs"))

write.csv2(taxa_info, paste0(project_path,"/taxa_info.csv"), row.names=FALSE)


#--------------------------------------------
#---- Clean up environment and local disk----
#--------------------------------------------
# Remove the zipped folder
suppressWarnings(file.remove(paste0(raw.path, ".zip"), full.names = TRUE))

# Remove the unzipped folder 
unlink(raw.path, recursive = TRUE)

rm(list=ls())



