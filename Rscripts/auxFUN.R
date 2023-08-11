
# -----------------------------------------------------------------------------
#
# Elaborado por Inder Tecuapetla, May 31, 2023
#
# Funciones auxiliares
# 
# Hecho para SELPER/CEOS Working Group Chapter D Training Group
#
# -----------------------------------------------------------------------------


# --- 

LoadToEnvironment <- function(RData, env = new.env()){
  load(RData, env)
  return(env) 
}

# ---

get_timeSeries_byClicking <- function(toPlot, df){

  nRow <- length(unlist(toPlot)) / 2

  mat_toPlot <- matrix(as.numeric(unlist(toPlot)), nrow = nRow)

  dX <- matrix(NA, nrow = nrow(df))

  dY <- matrix(NA, nrow = nrow(df))

  aproxX <- numeric(nRow)

  aproxY <- numeric(nRow)

  dX <- sapply(1:nRow, function(s) abs(df[,1] - mat_toPlot[s,1]))

  aproxX <- sapply(1:nRow, function(s) df[which.min(dX[,s]),1] )

  dY <- sapply(1:nRow, function(s) abs(df[,2] - mat_toPlot[s,2]))

  aproxY <- sapply(1:nRow, function(s) df[which.min(dY[,s]),2] )

  toExtract <- matrix(NA, nrow = nRow, ncol = 2)

  toExtract[,1] <- aproxX
  toExtract[,2] <- aproxY

  pixels <- matrix(NA, nrow = nRow, ncol = ncol(df)-2)

  for(i in 1:nRow){
    pixels[i,] <- df[(df[,1] == toExtract[i,1]) & (df[,2] == toExtract[i,2])][-c(1:2)]
  }

  list(ts = pixels, coord = toExtract)
}

# ---

ts_at_breaks <- function(ts, at){
  sapply(1:length(at), function(s) ts[at[s]] )
}

# getYear <- function(start=2003, end=2016, bp, freq=23){
#   period <- start:end
#   totalDays <- c(0, freq * 1:length(start:end))
#   
#   if( length(bp) == 1 ){
#     year <- period[sum( totalDays - bp < 0 )]
#   } else {
#     year <- unlist( lapply(1:length(bp), function(s) period[sum( totalDays - bp[s] < 0 )]  ) )
#   }
#   
#   year
# }

# ---

get_dNBR_dNDVI <- function(ndvi, nbr, breaks, before = 23,
                            after = 1, scale = TRUE, scaleFactor = 1e-4){
  
  validIndices <- FALSE
  validBreaks <- NULL
  dNBR <- NA
  dNDVI <- NA
  
  if(length(breaks) != 0){
    
    validBreaks <- breaks
    
    validIndices <- TRUE
    
    nbrPrevious <- ts_at_breaks(ts=nbr, at=breaks-before) 
    nbrPosterior <- ts_at_breaks(ts=nbr, at=breaks+after) 
    ndviPrevious <- ts_at_breaks(ts=ndvi, at=breaks-before)
    ndviPosterior <- ts_at_breaks(ts=ndvi, at=breaks+after)
    
    dNBR <- nbrPrevious - nbrPosterior
    dNDVI <- ndviPrevious - nbrPosterior
    
    if(scale){
      dNBR <- dNBR * scaleFactor
      dNDVI <- dNDVI * scaleFactor
    }
    
  }
  

list(dNBR=dNBR, dNDVI=dNDVI, validIndice=validIndices, validBreaks=breaks)
  
}

# ---

plot_ndvi_nbr_cps <- function(ndvi, nbr, ndvi_bfast){
  n_iter <- length(ndvi_bfast$output)
  
  AUX <- as.matrix(ndvi_bfast$output[[n_iter]]$ci.Vt$confint)
  
  draw_params <- c(5.1, 4.5, 0.5, 2.1)
  par(mfrow=c(2,1), cex.lab = 1.5, cex.axis = 1.25, 
      mar = draw_params)
  
  plot(as.numeric(ndvi$ts), col = "gray", ylab = "NDVI", xlab = "Years",
       lwd = 4, type = "l", axes=F)
  axis(1, at = seq(1,345, by=46), labels = seq(2003,2017, by=2))
  lines(1:length(ndvi$ts), as.numeric(ndvi_bfast$output[[n_iter]]$Tt), 
        col = "royalblue", lwd = 4)
  yRan <- range(ndvi$ts)
  yAxes <- seq(yRan[1], yRan[2], length=5)
  axis(2, at = yAxes, labels = round(yAxes, 2))
  TEMP <- unlist(sapply(1:nrow(AUX), function(s) AUX[s,1]:AUX[s,3]))
  for(i in 1:nrow(AUX)){
    axis(1, at = c(AUX[i,1], AUX[i,3]), labels = c("", ""), col = "red",
         lwd = 3, lwd.ticks = 4)
  }
  
  plot(as.numeric(nbr$ts), col = "gray", ylab = "NBR", xlab = "Years",
       lwd = 4, type = "l", axes=F)
  axis(1, at = seq(1,345,by=46), labels = seq(2003,2017,by=2))
  points(AUX[,2], nbr$ts[AUX[,2]], pch = 8, lwd = 3)
  for(i in 1:nrow(AUX)){
    points(AUX[i,2]-23, nbr$ts[AUX[i,2]-23], pch = 24, col = "blue", lwd = 3)
    points(AUX[i,2]+1, nbr$ts[AUX[i,2]+1], pch = 25, col = "blue", lwd = 3)
  }
  legend("bottomleft", pch = c(24,8,25), pt.cex = rep(2,3),
         col = c("blue","black","blue"), pt.lwd = c(2,2,2),
         legend = c("Before", "Change", "After"),
         # horiz = TRUE, 
         bty = "n")
  yRan <- range(nbr$ts)
  yAxes <- seq(yRan[1], yRan[2], length=5)
  axis(2, at = yAxes, labels = round(yAxes, 2))
  
}

# ---

vegCondition <- function(x){
  ifelse(x < -0.25, "High-regrowth", 
         ifelse(x < -0.1, "Low-regrowth", 
                ifelse(x < 0.1, "No burned",
                       ifelse(x < 0.27, "Low-severity",
                              ifelse(x < 0.66, "Moderate-severity", "High-severity")))))
  
}

# ---

get_SEVmap <- function(path, nameLayer){
  
  # path=tifDIR[5]
  # nameLayer = "2005"
  
  r <- raster(path) #raster(paste0(pathSeverityMaps, "/severityMap_0p15_", YEAR, ".tif"))
  
  r2 <- reclassify(r, c(-0.1,0.1,NA))
  
  df_r2 <- rasterToPoints(r2)
  
  datos <- matrix(NA, nrow = nrow(df_r2), ncol = 1)
  LONG <- matrix(NA, nrow = nrow(df_r2), ncol = 1)
  LAT <- matrix(NA, nrow = nrow(df_r2), ncol = 1)
  DIAMETER <- matrix(NA, nrow = nrow(df_r2), ncol = 1)
  
  MIN <- min(df_r2[,3])
  MAX <- max(df_r2[,3])
  for(i in 1:nrow(datos)){
    datos[[i]] <- df_r2[i,3]
    LONG[[i]] <- paste0("X: ", df_r2[i,1])
    LAT[[i]] <- paste0("Y: ", df_r2[i,2])
    DIAMETER[[i]] <- ((df_r2[i,3] - MIN)/(MAX-MIN))^20
  }
  #
  TEST <- data.frame(x = df_r2[,1], y = df_r2[,2], 
                     dNBR = datos, X = LONG, Y = LAT, 
                     diam = DIAMETER)
  TEST$X <- unlist(lapply(TEST$X, as.character))
  TEST$Y <- unlist(lapply(TEST$Y, as.character))
  coordinates(TEST) <- ~ x + y
  proj4string(TEST) <- "+proj=utm +zone=13 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
  
  myColorPal <- c("darkgreen", "palegreen4", "#CFE6B8", #"darkseagreen2", 
                  "gold", "orange", "firebrick1")
  
  myPal <- colorRampPalette(colors = myColorPal, alpha = TRUE)
  
  mp_severity <- mapview(TEST, legend = TRUE, 
                         layer.name = nameLayer,
                         col.regions = myPal(256),
                         na.color = "transparent", alpha = 0,
                         zcol = "dNBR", cex = "diam",
                         at = c(-1,-0.25,-0.1,0.1,0.27,0.66,1))
  
  mp_severity@map$x$calls[[11]]$args[[1]]$labels <- c("High regrowth",
                                                      "Low regrowth",
                                                      "Unburned",
                                                      "Low severity",
                                                      "Moderate severity",
                                                      "High severity")
  
  mp_severity@map$x$calls[[11]]$args[[1]]$na_label <- c("")
  
  mp_severity@map
}

# ---

# --- interpolation hybrid method

get_pixel_matrix <- function(x,lenPeriod=23){
  output <- matrix(nrow=length(x)/lenPeriod, ncol=lenPeriod)
  
  for(i in seq_len(nrow(output))){
    output[i,] <- x[((i-1) * lenPeriod + 1):(i * lenPeriod)]
  }
  output
}

# # --- this function fill gaps in first 3 years of the series stored in m
# # m is a matrix nrow=nYears, ncol=23
# fill_gap <- function(m){
#   apply(m[-1,1:3],MARGIN=2,FUN=median)
# }
# # ---

# NOTE: for this to work, length(x) must be a multiple of lenPeriod
climatology <- function(x, lenPeriod){
  MAT <- get_pixel_matrix(x=x, lenPeriod=lenPeriod)
  
  BOXPLOT <- boxplot(MAT, plot=FALSE)
  
  list(matrix=MAT, boxplot=BOXPLOT)
}

# =============================================================================
# --- Functions used for break-points estimation ---
# =============================================================================

#' Applies bfast01 to data to get a breakpoint estimate
#' 
#' start     numeric specifying starting year of data
#' end       numeric specifying ending year of data
#' frequency numeric givin number of observations per year
#' data      numeric vector, the data to abalyze
#' 
#' Details. This functions assumes that data has been regularly sampled. No missing
#' value is allowed - a priori. It is assumed that the first observation has been 
#' acquired on the first day of the start year.
#' 
#' Value
#' 
#' A list with objects
#' 
#' bPs          numeric vector with the estimated breakpoints
#' type         numeric vector specifying which type of breakpoint has been detected
#' significance numeric
#' stability    numeric
#' 
#' For further details see bfast01
#' 
getBreak <- function(data, start=2000, end=2018, frequency=23, bw=0.15){
  output <- NA
  breakType <- NA
  significance <- NA
  stability <- NA
  
  # data = pixel_full_type1; start=2000; end=2009; bw=0.2
  # frequency=23
  
  dataTS <- ts(data, start=c(start, 1), end=c(end, frequency), 
               frequency=frequency)
  
  getBFAST <- bfast01(data=dataTS, bandwidth=bw)
  
  # if(getBFAST$breaks == 1){
  output <- getBFAST$breakpoints
  temp <- bfast01classify(getBFAST)
  breakType <- as.numeric(temp[1])
  significance <- as.numeric(temp[2])
  stability <- as.numeric(temp[3])
  # }
  
  list(bPs=output, type=breakType, significance=significance, 
       stability=stability)
}

#' Calculates the year on which a breakpoint has occurred as a function of
#' start and end years of a time series, as well as the estimated breakpoint
#' and the number of observations per year.
#' 
#' start numeric giving starting year of analyzed period
#' end   numeric giving ending year of analyzed period
#' bp    numeric, a breakpoint as estimated by bfast or bfast01
#' freq  numeric giving the number of observations pear year
#' 
#' Value
#' 
#' A numeric 
#' 
# getYear <- function(start=2000, end=2018, bp, freq=23){
#   period <- start:end
#   totalDays <- c(0, freq * 1:length(start:end))
#   year <- period[sum( totalDays - bp < 0 )]
#   year  
# }

getYear <- function(start=2003, end=2016, bp, freq=23){
  period <- start:end
  totalDays <- c(0, freq * 1:length(start:end))
  
  if( length(bp) == 1 ){
    year <- period[sum( totalDays - bp < 0 )]
  } else {
    year <- unlist( lapply(1:length(bp), function(s) period[sum( totalDays - bp[s] < 0 )]  ) )
  }
  
  year
}

# ---

getMapview <- function(mapRaster, colPal, nameLayer,
                       typeQuery = c("mousemove", "click"),
                       label = TRUE){
  
  # mapRaster = RASTER
  # colPal = brewer.pal(256, color)
  # nameLayer = "mean"
  # typeQuery = "click"
  
  palTest <- colorRampPalette( colors = colPal, alpha = TRUE )
  palTest_rev <- colorRampPalette( colors = rev(colPal) )
  
  typeQuery <- match.arg(typeQuery)
  
  m <- mapview(mapRaster, na.color = "transparent",
               col.regions = palTest(256),
               map.types = c("Esri.WorldImagery",
                             "CartoDB.Positron",
                             "OpenStreetMap",
                             "OpenTopoMap"),
               query.type = typeQuery,
               label = label,
               homebutton = FALSE,
               layer.name = nameLayer, legend = TRUE)
  
  m_rev <- mapview(mapRaster, na.color = "transparent",
                   col.regions = palTest_rev(256),
                   map.types = c("Esri.WorldImagery",
                                 "OpenTopoMap"),
                   homebutton = FALSE,
                   layer.name = nameLayer, legend = TRUE)
  
  m@map[[1]]$calls[[10]]$args[[1]]$labels <- rev(m@map[[1]]$calls[[10]]$args[[1]]$labels)
  
  m@map[[1]]$calls[[10]]$args[[1]]$title = paste0("slope-", nameLayer) #"Slope"
  
  m@map[[1]]$calls[[10]]$args[[1]]$colors <- m_rev@map[[1]]$calls[[8]]$args[[1]]$colors
  
  m
}

plot_central_pnorm <- function(q){
  x <- seq(-5, 5, by = 0.1)
  y <- dnorm(x = x, mean = 0, sd = 1)
  
  xz <- x
  xz[x <= q] <- NA
  xz[x > -q] <- NA
  
  yz <- y
  yz[x <= q] <- NA
  yz[x > -q] <- NA
  
  yRan <- range(y, yz, na.rm = T)
  
  plot(x, y, type = "l", col = "red", lwd = 3, 
       ylab = "f(z)", xlab = "", ylim = yRan,
       xlim = xRan, cex.main = 0.8)
  par(new = TRUE)
  plot(x, yz, type = "h",  ylim = yRan, #xaxt = "n", yaxt = "n", 
       xlab = "z", ylab = "", 
       main = paste0("Cuantil: ", abs(round(q,2)), ". Área bajo curva: ", 1-2*pnorm(q)))
  abline(v=q, col="blue", lty=2)
  abline(v=-q, col="blue", lty=2)
}

# --- 

spRast_valueCoords <- function(spRaster, na_rm=FALSE){
  
  spPoints <- as.points(spRaster, na.rm=na_rm)
  
  spValues <- extract(spRaster, spPoints)
  
  DIM <- dim(spValues)
  
  spRasterToPoints <- as.matrix(spValues[1:DIM[1],2:DIM[2]])
  
  spCoords <- crds(spRaster, na.rm=na_rm)
  
  list(values=spRasterToPoints, coords=spCoords)  
}

gapfill_climatology <- function(y, box=c("lower", "median", "upper"),
                                gapType=c(NA, NaN), lenPeriod=23){
  
  if(is.na(gapType)){
    toFill_ind <- (1:length(y))[is.na(y)]
  }
  
  if(is.nan(gapType)){
    toFill_ind <- (1:length(y))[is.nan(y)]
  }
  
  output <- y
  
  clima <- climatology(x=y, lenPeriod=lenPeriod)
  
  box <- match.arg(box)
  
  quant <- ifelse(box=="lower", 2, ifelse(box=="median", 3, 4))
  
  output[toFill_ind] <- clima$boxplot$stats[quant,toFill_ind%%lenPeriod]
  
  list(filled=output, original=y)
}



trendAnalysis <- function(x, startYear, endYear, frequency, productName){
  
  pixel_MK <- mk.test(x)
  
  # pixel_MK$p.value
  
  pixel_SenTheil <- sens.slope(x)
  
  # pixel_SenTheil$estimates
  
  b_hat <- as.numeric(pixel_SenTheil$estimates)
  a_hat <- median( x - b_hat * 1:length(x) )
  
  lineaTheilSen <- ts(a_hat + b_hat * 1:length(x), 
                      start = c(startYear, 1),
                      end = c(endYear, 12), frequency = frequency)
  
  x_ts <- ts(x, start=c(startYear,1), end=c(endYear,frequency), 
             frequency = frequency)
  
  
  par(mfrow=c(1,1), mar = c(2,2,1,2), adj=0)
  plot(x_ts , type="l", col = "gray", ylab = productName)
  lines(lineaTheilSen, lwd = 5, col = "lightcoral")
  legend("topright", legend = c("raw data", "linear trend"),
         col = c("gray", "lightcoral"), lty = rep(1,2), lwd = c(1,5), bty = "n")
  

}


