
#-----------------------------------------------------------------------------------
#This function calculates the number of decimal places in any given numeric value 
# eg., 15.21 has 2 decimal places, 15.2569 has 4 decimal places, 15.25690 also has 4, as 0 in the end doesn't count
#-----------------------------------------------------------------------------------
decimalplaces <- function(x) {
  if (abs(x - round(x)) > .Machine$double.eps^0.5) {
    # Remove trailing zeros and split at the decimal point
    split_result <- strsplit(sub('0+$', '', as.character(x)), ".", fixed = TRUE)[[1]]
    # Check if there are any decimals
    if (length(split_result) > 1) {
      nchar(split_result[[2]]) # Count characters in the decimal part
    } else {
      return(0) # No decimal part
    }
  } else {
    return(0) # No decimals for whole numbers
  }
}

#-----------------------------------------------------------------------------------
#Divide a numerical value by 10
#-----------------------------------------------------------------------------------
divide10<-function(x){
  value<-x/10
  return(value)
}

#-----------------------------------------------------------------------------------
#Divide occurrence column with either y=0 (absences) or y=1 (presences)
#-----------------------------------------------------------------------------------
add.occ<-function(x,y){
  occ<-rep(y,nrow(x))
  cbind(x,occ)
}

#-----------------------------------------------------------------------------------
#Function to return threshold where sens=spec from caret results 
#-----------------------------------------------------------------------------------
findThresh<-function(df){
  df<-df[c("rowIndex","obs","present")]
  df<-df %>%
    dplyr::mutate(observed= ifelse(obs == "present",1,0)) %>%
    dplyr::select(rowIndex,observed,predicted=present)
  result<-PresenceAbsence::optimal.thresholds(df,opt.methods = 2)
  return(result)
}


#-----------------------------------------------------------------------------------
#Recalculate accuracy for a given model with the threshold that has been optimized
#-----------------------------------------------------------------------------------
accuracyStats<-function(df,y){
  df<-df[c("rowIndex","obs","present")]
  df<-df %>%
    dplyr::mutate(observed= ifelse(obs == "present",1,0)) %>%
    dplyr::select(rowIndex,observed,predicted=present)
  result<-PresenceAbsence::presence.absence.accuracy(df,threshold = y,st.dev=FALSE)
  return(result)
}


#-----------------------------------------------------------------------------------
# Model predictions for a large raster in a more efficient way using parallellization 
#-----------------------------------------------------------------------------------
predict_large_raster<-function(rasterstack, model, type) {
  
  # Ensure that connections are closed even in case of an error
  on.exit({
    plan(strategy = "sequential")  # Ensure that the parallel plan is returned to sequential
    gc()  # Trigger garbage collection
    closeAllConnections()  # Close any open file connections
  }, add = TRUE)
  
  gc() #Free up memory
  
  ncores<-min(4, availableCores()-2)  #Set up number of cores
  
  if(class(rasterstack)!="SpatRaster"){
    raster_terra<-rast(rasterstack)  #Convert raster to terra raster format if not already
  }else{
    raster_terra<-rasterstack
  }
  
  chunk_size <- ceiling(nrow(raster_terra) / ncores)   # Define chunk size
  
  # Create a list of row indices for each chunk
  chunk_indices <- split(seq_len(nrow(raster_terra)), ceiling(seq_along(seq_len(nrow(raster_terra))) / chunk_size))
  
  # Extract raster chunks and put raster chunks in list
  r_list<- vector(mode = "list", length = ncores)
  for (core in 1:ncores) {
    r_list[[core]]<- wrap(raster_terra[min(chunk_indices[[core]]):max(chunk_indices[[core]]), ,drop=FALSE])
  } #SpatRasters need to be wrapped before sending out to different cores
  
  # Save model to disk if it’s large
  saveRDS(model, "model.rds")
  options(future.globals.maxSize = 4.5 * 1024^3)
  plan(strategy = "multisession", workers=ncores) #Set up parallel
  
  out_list <- future_lapply(r_list,  function(chunk) {
    model <- readRDS("model.rds")  # Load model from disk
    unwrapped_raster <- unwrap(chunk)  # Unwrap raster for processing
    predicted_raster <- predict(unwrapped_raster, model, type = type, na.rm = TRUE)
    rm(unwrapped_raster)
    wrap(predicted_raster)  # Wrap the raster again
  }, future.seed = TRUE)
  
  
  plan(strategy = "sequential")   #Close parallel processing
  file.remove("model.rds")
  rm(r_list) #Remove large objects we don't need anymore
  out_list<- lapply(out_list, unwrap) #unwrap chunks
  gc() # Clean up memory after processing
  model_parallel<- do.call(terra::merge, out_list)  # Merge the chunks 
  rm(out_list) #Remove large objects we don't need anymore
  gc()  #Final garbage collect
  options(future.globals.maxSize = 500 * 1024^2)  # Reset to 500 MB
  return(model_parallel)
}


#-----------------------------------------------------------------------------------
# PDF export function
#-----------------------------------------------------------------------------------
exportPDF<-function(rst,taxonkey,taxonName,nameextension,is.diff="FALSE"){
  filename=file.path(PDF_folder,paste("be_",taxonkey, "_",nameextension,sep=""))
  pdf(file=filename,width=10,height=8,paper="a4r")
  par(bty="n")#to turn off box around plot
  ifelse(is.diff=="TRUE", brks<-seq(-1, 1, by=0.2), brks <- seq(0, 1, by=0.1)) 
  nb <- length(brks)-1 
  pal <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))
  cols<-pal(nb)
  maintitle<-paste(taxonName,taxonkey,"_",nameextension, sep= " ")
  plot(rst, breaks=brks, col=cols,main=maintitle, lab.breaks=brks,axes=FALSE)
  dev.off() 
} 


#-----------------------------------------------------------------------------------
# Export PNG function
#-----------------------------------------------------------------------------------
exportPNG<-function(rst,taxonkey,taxonName,nameextension,is.diff="FALSE"){
  filename=file.path(pdfOutput,paste("be_",taxonkey, "_",nameextension,sep=""))
  png(file=filename)
  par(bty="n")#to turn off box around plot
  ifelse(is.diff=="TRUE", brks<-seq(-1, 1, by=0.25), brks <- seq(0, 1, by=0.1)) 
  nb <- length(brks)-1 
  pal <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))
  cols<-pal(nb)
  maintitle<-paste(taxonName,taxonkey,"_",nameextension, sep= " ")
  plot(rst, breaks=brks, col=cols,main=maintitle, lab.breaks=brks,axes=FALSE)
  dev.off() 
} 


#-----------------------------------------------------------------------------------
# Generate pseudoabsences
#-----------------------------------------------------------------------------------
generate_pseudoabs <- function(index = NULL,mask, alternative_mask, n, p) {
  tryf_values <- c(50,100, 150)  # tryf values to attempt in each stage
  current_raster <- mask  # Start with the initial raster layer
  
  # Attempt to generate points
  for (tryf in tryf_values) {
    # Generate random points
    suppressWarnings(pseudoabs <- as.data.frame(
      randomPoints(
        current_raster, 
        n, 
        p, 
        ext = NULL, 
        extf = 1.1, 
        excludep = TRUE, 
        prob = FALSE, 
        cellnumbers = FALSE, 
        tryf = tryf, 
        warn = 2, 
        lonlatCorrection = TRUE
      )
    )
    )
    # Check if the number of pseudoabsences reaches required amount
    if (nrow(pseudoabs) == n) {
      # If index is provided, include it in the message (only for lists)
      if (!is.null(index)) {
        message(paste0(n, " out of ", n, " pseudoabsences generated while accounting for observer bias in set ", index))
      } else {
        message(paste0(n, " out of ", n, " pseudoabsences generated while accounting for observer bias."))
      }
      return(pseudoabs)  # Return dataset if the required amount of pseudoabsences are generated
    }
  }
  
  # If unsuccessful with biasgrid ecoregions raster, switch to the full ecoregions raster and retry
  current_raster <- alternative_mask
  
  for (tryf in tryf_values) {
    pseudoabs <- as.data.frame(
      randomPoints(
        current_raster, 
        n, 
        p, 
        ext = NULL, 
        extf = 1.1, 
        excludep = TRUE, 
        prob = FALSE, 
        cellnumbers = FALSE, 
        tryf = tryf, 
        warn = 2, 
        lonlatCorrection = TRUE
      )
    )
    
    # Check if the number of rows meets the desired count
    if (nrow(pseudoabs) == n) {
      # If index is provided, include it in the warning (only for lists)
      if (!is.null(index)) {
        warning(paste0(n, " out of ", n, " pseudoabsences generated without accounting for observer bias in set ", index))
      } else {
        warning(paste0(n, " out of ", n, " pseudoabsences generated without accounting for observer bias."))
      }
      return(pseudoabs)  # Return dataset if enough pseudoabsences were generated
    }
  }
  
  # If all attempts fail, return the last generated dataframe with fewer pseudoabsences than requested
  # If index is provided, include it in the warning
  if (!is.null(index)) {
    warning(paste0("Could not generate the required number of pseudoabsences: ", n, " out of ", n, " pseudoabsences generated without accounting for observer bias in set ", index))
  } else {
    warning(paste0("Could not generate the required number of pseudoabsences: ", n, " out of ", n, " pseudoabsences generated without accounting for observer bias."))
  }
  
  return(pseudoabs)  # Return the pseudoabs data, even if incomplete
}


#-----------------------------------------------------------------------------------
# Recode factor levels to absent (0) and present(1), and set present as the reference level
#-----------------------------------------------------------------------------------
factorVars<-function(df,var){
  df[,c(var)]<-as.factor(df[,c(var)])
  levels(df[,c(var)])<-c("absent","present")
  df[,c(var)]<-relevel(df[,c(var)], ref = "present")
  return(df)
}


#-----------------------------------------------------------------------------------
#----------------Create folders when they don't exist yet---------------------------
#-----------------------------------------------------------------------------------
create_folder <- function(path, name) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
    message(paste0("Folder '", name, "' created at path: '", path, "' 🎉"))
  } else {
    message(paste0("Folder '", name, "' already exists at path: '", path, "' 🎉"))
  }
}

#-----------------------------------------------------------------------------------
#-------------------------------Functions script 4----------------------------------
#-----------------------------------------------------------------------------------

stdres<-function(obs.numeric, yhat){
  num<-obs.numeric-yhat
  denom<-sqrt(yhat*(1-yhat))
  return(num/denom)
}


# functions needed for conformal prediction function


GetLength<-function(x,y){
  length(x[which(x<= y)])
}



# Code for Mondrian conformal prediction functions



CPconf<-function(pA,pB,confidence){
  if(pA > confidence && pB< confidence){
    predClass<-"classA"
  }else if(pA < confidence && pB> confidence){
    predClass<-"classB"
  }else if(pA< confidence && pB< confidence){
    predClass<-"noClass"
  }else{
    predClass<-"bothClasses"
    
    return(predClass)
  }}


#function to calculate confidence of each prediction

get.confidence<-function(pvalA,pvalB){
  secondHighest<-ifelse(pvalA>pvalB,pvalB,pvalA)
  conf<-(1-secondHighest)
  return(conf)
}

forcedCp<-function(pvalA,pvalB){
  ifelse(pvalA>pvalB,"presence","absence")
}

extractVals<-function(predras){
  vals <-  raster::values(predras)
  coord <-  raster::xyFromCell(predras,1:ncell(predras))
  raster_fitted <- cbind(coord,vals)
  raster_fitted.df<-as.data.frame(raster_fitted)
  raster_fitted.df1<-na.omit(raster_fitted.df)
  raster_fitted.df1$presence<-raster_fitted.df1$lyr1
  raster_fitted.df1$absence<- (1-raster_fitted.df1$presence)
  return(raster_fitted.df1)
}


classConformalPrediction<-function(x,y){
  ens_results<- get("x")
  ens_calib<-ens_results$ens_model$pred
  calibPresence<-ens_calib %>%
    filter(obs=='present')%>%
    select(present)
  calibPresence<-unname(unlist(calibPresence[c("present")]))
  calibAbsence<-ens_calib %>%
    filter(obs=='absent')%>%
    select(absent)
  calibAbsence<-unname(unlist(calibAbsence[c("absent")]))
  predicted.values<-extractVals(y)
  
  
  testPresence<-predicted.values$presence
  testAbsence<-predicted.values$absence
  
  #derive p.Values for class A
  smallrA<-lapply(testPresence,function(x) GetLength(calibPresence,x))
  smallrA_1<- unlist (smallrA)+1
  nCalibSet<-length(calibPresence)+1
  pvalA<-smallrA_1+1/nCalibSet
  
  # derive p.Values for Class B
  smallrB<-lapply(testAbsence,function(x) GetLength(calibAbsence,x))
  smallrB_1<- unlist (smallrB)+1
  nCalibSetB<-length(calibAbsence)
  pvalB<-smallrB_1/nCalibSetB
  
  pvalsdf<-as.data.frame(cbind(pvalA,pvalB,0.20))
  #raster_cp_20<-mapply(CPconf,pvalsdf$pvalA,pvalsdf$pvalB,pvalsdf[3])
  #table(raster_cp_20)
  
  pvalsdf$conf<-get.confidence(pvalsdf$pvalA,pvalsdf$pvalB)
  pvalsdf_1<-cbind(pvalsdf,predicted.values)
}


# Confidence maps

confidenceMaps<-function(x,taxonkey,taxonName,maptype){
  pvals_dataframe<-get("x")
  data.xyz <- pvals_dataframe[c("x","y","conf")]
  rst <- rasterFromXYZ(data.xyz)
  crs(rst)<-CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs") 
  plot(rst,breaks=brks, col=cols,lab.breaks=brks)
  writeRaster(rst, filename=file.path(raster_folder,paste("be_",taxonkey, "_",maptype,".tif",sep="")),overwrite=TRUE)
  exportPDF(rst,taxonkey,taxonName=taxonName,nameextension= paste(maptype,".pdf",sep=""))
  return(rst)
}

### Generate and export response curves in order of variable importance

partial_gbm<-function(x){
  m.gbm<-pdp::partial(bestModel$models$gbm$finalModel,pred.var=paste(x),train = bestModel.train,type="classification",
                      prob=TRUE,n.trees= bestModel$models$gbm$finalModel$n.trees, which.class = 1,grid.resolution=nrow(bestModel.train))
}

partial_glm<-function(x){
  m.glm<-pdp::partial(bestModel$models$glm$finalModel,pred.var=paste(x),train = bestModel.train,type="classification",
                      prob=TRUE,which.class = 1,grid.resolution=nrow(bestModel.train))
}

partial_rf<-function(x){
  pdp::partial(bestModel$models$rf$finalModel,pred.var=paste(x),train = bestModel.train,type="classification",
               prob=TRUE,which.class = 1,grid.resolution=nrow(bestModel.train))
}

partial_mars<-function(x){
  m.mars<-pdp::partial(bestModel$models$earth$finalModel,pred.var=paste(x),train = bestModel.train,type="classification",
                       prob=TRUE,which.class = 2,grid.resolution=nrow(bestModel.train)) # class=2 because in earth pkg, absense is the first class
}


responseCurves<-function(x,y) {
  colors <- c("GLM" = "gray", "GBM"="red","RF"="blueviolet","MARS"= "hotpink") 
  ggplot(all_dfs,(aes(x=.data[[x]],y=.data[[y]]))) +
    geom_line(aes(color = data), size =1.2, position=position_dodge(width=0.2))+
    theme_bw()+
    labs(y="Partial probability", x= gsub("//..*","",x),color="Legend") +
    scale_color_manual(values = colors)
} 


eu_eval<-function (ras,y){
  indep.bil<-raster::extract(ras,y,method="bilinear")
  indep.bil.df<-as.data.frame(indep.bil)
  indep.bil.df<-indep.bil.df %>%
    mutate(predicted= ifelse(indep.bil >= 0.5,"present","absent")) 
  indep.bil.df$observed<-rep("present",nrow(indep.bil.df))
  indep.bil.df$predicted<-as.factor(indep.bil.df$predicted)
  indep.bil.df$observed<-as.factor(indep.bil.df$observed)
  xtab<-table(indep.bil.df$predicted,indep.bil.df$observed)
  return(xtab)
}
