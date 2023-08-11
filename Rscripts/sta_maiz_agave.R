
# --- Diplomado Geomática, IG, UNAM, 2023
# --- Módulo X: Análisis de series de tiempo de imágenes satelitales con R

# --- Elaborado: Agosto 11, 2023
# --- Impartido: Agosto 11, 2023

# --- En este script hacemos:
# --- Análisis de tendencias estacionales
# --- Consultar https://conabio.shinyapps.io/viSTA_esp/ para una explicacion
# --- interactiva de sta

# --- DATASET: LP_NDMI_S2_2019_2022.tif pegado en la nube del diplomado


library(terra)
library(rgdal)
library(raster)
library(mapview)
library(bfast)
library(Kendall)
library(trend)
library(sta)

library(doParallel)
library(foreach)

source( paste0( getwd(), "/Rscripts/auxFUN.R" ) )

# ---

tifFILES <- list.files(path = paste0( getwd(), "/data/LaPiedad_Mich" ),
                       pattern = ".tif",
                       full.names = TRUE)

shpFILES <- list.files(path = paste0( getwd(), "/data/muni_2018" ),
                       pattern = ".shp",
                       full.names = TRUE)

LPsentinel <- stack(tifFILES)
shp_LP <- shapefile(shpFILES)

shp_LP_utm <- spTransform(shp_LP[855,], crs(LPsentinel))

LPsentinel_shp <- raster_intersect_sp(LPsentinel, shp_LP_utm)

LPsentinel_shp_rTp <- rasterToPoints(LPsentinel_shp)

# ---

plot(subset(sentinel_LP_shp,25))

# Agricola
# 20.365365, -102.024281

# xy <- SpatialPoints(cbind(-102.024281, 20.365365))
xy <- SpatialPoints(cbind(-102.0216491, 20.365379))
proj4string(xy) <- CRS('+init=epsg:4326')
xy_UTM <- spTransform(xy, CRS('+init=epsg:32613'))

plot(xy_UTM, add=TRUE, pch="M")

# Agave
# 20.363361, -102.034041

xyT <- SpatialPoints(cbind(-102.034041, 20.363361))
proj4string(xyT) <- CRS('+init=epsg:4326')
xyT_UTM <- spTransform(xyT, CRS('+init=epsg:32613'))

plot(xyT_UTM, add=TRUE, pch="A")

# --- 

xy_list <- list(x=extent(xy_UTM)[1], y=extent(xy_UTM)[3])

pixel_homemade <- get_timeSeries_byClicking(toPlot = xy_list, 
                                            df = LPsentinel_shp_rTp)

pixel_filled <- gapfill_climatology(y=as.numeric(pixel_homemade$ts),
                                    box = "median", gapType = NaN)

pixel_sta <- sta(pixel_filled$filled,
                 numFreq = 3,
                 startYear = 2019,
                 endYear = 2022,
                 freq = 12)

plot(pixel_sta)

LPsentinel_shp_rTp[137,]
# ---

# ---------------------------------------------------
# --- Datos faltantes: Imputación a gran escala --- #
# ---------------------------------------------------

# --- CODIGO EN PARALELO

# filled_layers guardará las imputaciones
filled_layers <- matrix(nrow=nrow(LPsentinel_shp_rTp), 
                        ncol=ncol(LPsentinel_shp_rTp))
filled_layers[,1:2] <- LPsentinel_shp_rTp[,1:2]

numCores <- detectCores()

# --- Asegurarse de crear folder /RData/progressReports

# progress report file (to check out on the process)
progressReportFile <- paste0(getwd(), "/RData/progressReports/LP_clima_Imputation.txt" )
file.create(path=progressReportFile, showWarnings=FALSE)

write("===CLIMATOLOGY imputation began at===",
      file=progressReportFile, append=TRUE)
write(as.character(Sys.time()[1]), file=progressReportFile,
      append=TRUE)

kluster <- parallel::makeCluster(numCores-1, outfile="")
registerDoParallel(kluster)

output <- foreach(i=1:nrow(LPsentinel_shp_rTp), .combine="rbind", .packages = "imputeTS") %dopar% { # nrow(sp_ndvi_rTp)
  
  pixel <- LPsentinel_shp_rTp[i, 3:ncol(LPsentinel_shp_rTp)]
  
  numNaN <- sum( is.nan(pixel) )
  
  if(numNaN == 1){
    pixel_filled <- gapfill_climatology(y=as.numeric(pixel),
                                        box = "median", gapType = NaN)$filled
  } else {
    pixel_filled <- na_interpolation(x=as.numeric(pixel))
  }
  
  s <- pixel_filled
  
  if(i %% 100 ==0){
    texto <- paste0("Working on ROW: ", i)
    write(texto, file=progressReportFile, append=TRUE)
  }
  
  return(s)
}
stopCluster(kluster)

write( as.character(Sys.time()[1]), file=progressReportFile, append=TRUE)
write( "===CLIMATOLOGY imputation ended here===", file=progressReportFile, 
       append=TRUE)
# ---

# --- Guardando las imputaciones como objetos matrix
filled_layers[,3:50] <- output[,1:48]

# --- Asegurarse de crear folder /RData/LaPiedad_imputation

save(filled_layers, file=paste0(getwd(),"/RData/LaPiedad_imputation/LP_NDMI_S2_2019_2022.RData"))
# ---

# --- Leer docu sta

LP_sta <- sta(data=filled_layers,
              freq = 12, numFreq = 3, 
              startYear = 2019, endYear = 2022,
              save=TRUE, 
              dirToSaveSTA = paste0( getwd(), "/RData/LaPiedad_sta" ))

MASTER <- getMaster(x=subset(LPsentinel_shp,1))

plot(x=LP_sta, master=MASTER, significance=0.1)

LP_sta_matrix <- LoadToEnvironment(paste0( getwd(), "/RData/LaPiedad_sta/sta_matrix_output.RData" ))

ls(LP_sta_matrix)

# -----------------------
# --- RASTERIZACION --- #
# -----------------------

media <- LP_sta_matrix$mean
media_slope <- media[,c(1,2,3)]
media_pval <- media[,c(1,2,4)]

map_media_slope <- matrixToRaster(media_slope, 
                                  projection = crs(LPsentinel_shp))

map_media_pval <- matrixToRaster(media_pval, 
                                 projection = crs(LPsentinel_shp))

# ---

annual <- LP_sta_matrix$annual

annual_slope <- annual[,c(1,2,3)]
annual_pval <- annual[,c(1,2,4)]

map_annual_slope <- matrixToRaster(annual_slope, 
                                   projection = crs(LPsentinel_shp))

map_annual_pval <- matrixToRaster(annual_pval, 
                                  projection = crs(LPsentinel_shp))

# ---

s_annual <- LP_sta_matrix$semiannual
s_annual_slope <- s_annual[,c(1,2,3)]
s_annual_pval <- s_annual[,c(1,2,4)]

map_s_annual_slope <- matrixToRaster(s_annual_slope, 
                                     projection = crs(LPsentinel_shp))

map_s_annual_pval <- matrixToRaster(s_annual_pval, 
                                    projection = crs(LPsentinel_shp))

# -----------------------
# --- VISUALIZACION --- #
# -----------------------

map_media <- map_media_slope
map_media[ map_media_pval > 0.09 ] <- NA

plot(map_media)

mp_media <- mapview(map_media, na.color="transparent")

mp_media

mp_M <- mapview(xy_UTM)

mp_A <- mapview(xyT_UTM)

mp_media + mp_M + mp_A

map_annual <- map_annual_slope
map_annual[ map_annual_pval > 0.09 ] <- NA

plot(map_annual)

mp_annual <- mapview(map_annual, na.color="transparent")

mp_annual + mp_M + mp_A

map_s_annual <- map_s_annual_slope
map_s_annual[ map_s_annual_pval > 0.09 ] <- NA

plot(map_s_annual)

mp_s_annual <- mapview(map_s_annual, na.color="transparent")

mp_s_annual + mp_M + mp_A

# ---

rgb_sta <- stack(map_media, map_annual, map_s_annual)

viewRGB(rgb_sta, na.color="transparent")
