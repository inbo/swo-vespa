######################################
###Creating bias layer for insects ###
######################################

#https://scottrinnan.wordpress.com/2015/08/31/how-to-construct-a-bias-file-with-r-for-use-in-maxent-modeling/


#--------------------------------------------
#-----------------Load packages--------------
#--------------------------------------------
packages <- c("rgbif", "dplyr", "purrr", "assertthat", "readr", "here", "qs", "dismo", "raster", "MASS", "magrittr", "maptools", "data.table")

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
#specify project name
projectname<-"Test_Vespa_velutina_12_12"

# Define the folder paths
project_path <- file.path("./data/projects",project)
raw_path<-file.path("./data/raw")

# Check and create each folder if necessary
create_folder(project_path, project)
create_folder(raw_path, "raw")



#--------------------------------------------
#-----------Retrieve GBIF taxonkeys----------
#--------------------------------------------

# specify the scientific name of the species to be modelled
class<-c("Insecta")

# Match species names with the GBIF backbone, retrieve taxon keys from GBIF when a match is found
taxon_df <- as.data.frame(class)

mapped_taxa <- purrr::map_dfr(
  taxon_df$class,
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
  dplyr::filter(rank =="CLASS")

#Make sure that all species were mapped to the GBIF backbone, if not an error will appear indicating which species are missing
assertthat::assert_that(nrow(mapped_taxa)==length(class),
                        msg=paste0("The following species could not be found in the GBIF backbone taxonomy: "
                                   ,class[!sapply(class, function(x) any(grepl(x,mapped_taxa$scientificName)))])
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
#All basis of record types, except `FOSSIL SPECIMEN` and `LIVING SPECIMEN`, which can have misleading location information (e.g. location of captive animal).
basis_of_record <- c(
  "OBSERVATION", 
  "HUMAN_OBSERVATION",
  "MATERIAL_SAMPLE",
  "PRESERVED_SPECIMEN", 
  "UNKNOWN", 
  "MACHINE_OBSERVATION",
  "OCCURRENCE"
)

#Time period
year_begin <- 2004
year_end <-2024

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



#--------------------------------------------
#--------------Retrieve download-------------
#--------------------------------------------

gbif_download_key<-"0037565-241126133413365"
occ_download_get(gbif_download_key, path = here("data","raw"), overwrite=TRUE)
metadata <- occ_download_meta(key = gbif_download_key)
gbif_download_key<-metadata$key

#extract_GBIF_occurrence
raw.path<- here("data", "raw", gbif_download_key)
unzip(paste0(raw.path,".zip"),exdir=raw.path, overwrite=TRUE)
occurrence_file<-as.data.frame(data.table::fread(paste0(raw.path,"/occurrence.txt"),header=TRUE)) #deze twee lijnen samenbrengen
locations <- fread(paste0(raw.path,"/occurrence.txt"),header=TRUE, select = c("decimalLongitude", "decimalLatitude"))            #deze twee lijnen samenbrengen

# Clean data by removing rows with NA values in the required columns
locations <- na.omit(locations)
# Remove rows with missing values (NA) in either of those columns
locations <- na.omit(locations)


#--------------------------------------------
#------------Load climate rasters------------
#--------------------------------------------

# Climate rasters are used as reference for their grid size. Same grid size will be used for bias layers
climdat <- brick("./data/external/climate/trias_CHELSA/CHELSA_annPrecip_12.tif")  
occur.ras <- rasterize(locations, climdat, fun = "length")
plot(occur.ras)



presences <- which(values(occur.ras) == 1)
pres.locs <- coordinates(occur.ras)[presences, ]

dens <- kde2d(pres.locs[,1], pres.locs[,2], n = c(nrow(occur.ras), ncol(occur.ras)))
dens.ras <- raster(dens)
plot(dens.ras)


writeRaster(dens.ras, "./data/external/bias_grids/final/trias/insects_1km.tif")

occurrences <- read.csv(occurdat[1])
mod1 <- maxent(climdat, occurrences, args = "biasfile=dens.ras")