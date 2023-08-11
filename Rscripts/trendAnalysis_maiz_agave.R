# --- Diplomado Geomática, IG, UNAM, 2023
# --- Módulo X: Análisis de series de tiempo de imágenes satelitales con R

# --- Elaborado: Agosto 10, 2023
# --- Impartido: Agosto 11, 2023

# --- En este script hacemos:
# --- Análisis de tendencias y de cambio abrupto para todo el dataset

# --- DATASET: LP_NDMI_S2_2019_2022.tif pegado en la nube del diplomado

library(terra)
library(mapview)
library(bfast)
library(Kendall)
library(trend)
library(foreach)
library(doParallel)
library(raster)
library(geoTS)

source( paste0( getwd(), "/Rscripts/auxFUN.R" ) )

# ---

tifFILES <- list.files(path = paste0( getwd(), "/data/LaPiedad_Mich" ),
                       pattern = ".tif",
                       full.names = TRUE)

sentinel_LP <- rast(tifFILES)

shpFILES <- list.files(path = paste0( getwd(), "/data/muni_2018" ),
                       pattern = ".shp",
                       full.names = TRUE)

# --- leyendo shapefile con terra

shp_LP <- vect(shpFILES)

shp_LP@ptr

str(shp_LP@ptr$df)

shp_LP@ptr$df$names[450,]

shp_LP[450,]$NOM_MUN

poly <- (1:length(shp_LP))[shp_LP$NOM_MUN == "La Piedad"]

shp_LP_polygon <- shp_LP[poly,]

# --- cambio de proyeccion con terra

shp_LP_lon_lat <- project(x=shp_LP_polygon, y=sentinel_LP)

# --- haciendo un recorte con terra

sentinel_LP_shp <- terra::crop(sentinel_LP, shp_LP_lon_lat)

plot(subset(sentinel_LP_shp, 20))

sentinel_LP_shp <- terra::mask(sentinel_LP_shp, shp_LP_lon_lat)

plot(subset(sentinel_LP_shp, 20))

# --- extrayendo los "numeritos" con version propia de rasterToPoint para terra

sentinel_LP_shp_rTp <- spRast_valueCoords(spRaster=sentinel_LP_shp, na_rm=TRUE)

# ---

mat_pvalue <- matrix(nrow=nrow(sentinel_LP_shp_rTp$values), ncol=3)
mat_pvalue[,1:2] <- sentinel_LP_shp_rTp$coords

mat_slope <- matrix(nrow=nrow(sentinel_LP_shp_rTp$values), ncol=3)
mat_slope[,1:2] <- sentinel_LP_shp_rTp$coords

numCores <- detectCores()

# --- progress report file (to check out on the process)
progressReportFile <- paste0(getwd(), "/RData/progressReports/LP_trendAnalysis.txt" )
file.create(path=progressReportFile, showWarnings=FALSE)

write("===TREND analysis began at===",
      file=progressReportFile, append=TRUE)
write(as.character(Sys.time()[1]), file=progressReportFile,
      append=TRUE)

kluster <- parallel::makeCluster(numCores-1, outfile="")
registerDoParallel(kluster)

output <- foreach(i=1:nrow(sentinel_LP_shp_rTp$values), 
                  .combine="rbind",
                  .packages=c("trend") ) %dopar% { # nrow(sp_ndvi_rTp)
                    
                    pixel <- as.numeric(sentinel_LP_shp_rTp$values[i, ]) * 1e-4
                    
                    # pixel_interpol <- na_interpolation(pixel * 1e-4)
                    pixel_filled <- gapfill_climatology(y=pixel, 
                                                        box="median", gapType = NaN, 
                                                        lenPeriod = 12)
                    
                    pixel_MannKendall <- mk.test(pixel_filled$filled)$p.value
                    
                    pixel_SenTheil <- sens.slope(pixel_filled$filled)$estimates
                    
                    s <- c(as.numeric(pixel_MannKendall), as.numeric(pixel_SenTheil))
                    
                    if(i %% 100 ==0){
                      texto <- paste0("Working on ROW: ", i)
                      write(texto, file=progressReportFile, append=TRUE)
                    }
                    
                    return(s)
                  }
stopCluster(kluster)

write( as.character(Sys.time()[1]), file=progressReportFile, append=TRUE)
write( "===TREND analysis ended here===", file=progressReportFile, append=TRUE)
# ---

# --- Guardando los estadísticos como objetos matrix
mat_pvalue[,3] <- output[,1]
mat_slope[,3] <- output[,2]

# --- Asegurarse de crear /LaPiedad_trendAnalysis en /RData

save(mat_pvalue, file=paste0(getwd(),"/RData/LaPiedad_trendAnalysis/pvalue.RData"))
save(mat_slope, file=paste0(getwd(),"/RData/LaPiedad_trendAnalysis/slope.RData"))
# ---

# -----------------------
# --- RASTERIZACION --- #
# -----------------------

PROJECTION <- raster::projection(raster(tifFILES[1]))

# map_pvalue <- LoadToEnvironment(paste0(getwd(),"/RData/LaPiedad_trendAnalysis/pvalue.RData"))$df_pvalue
# map_slope <- LoadToEnvironment(paste0(getwd(),"/RData/LaPiedad_trendAnalysis/slope.RData"))$df_slope

map_pvalue <- matrixToRaster(matrix=mat_pvalue, projection=PROJECTION)
map_slope <- matrixToRaster(matrix=mat_slope, projection=PROJECTION)

# --- Raster de pixeles con tendencia significativa

pvalueMap <- map_pvalue
pvalueMap[pvalueMap > 0.05] <- NA

slopeMap <- map_slope
slopeMap[is.na(pvalueMap)] <- NA

plot(slopeMap, main="Pendiente significativa estimada al 95% de significancia")

# --- Inspección visual espacial

shp_raster <- shapefile(shpFILES[1])

mp_shp <- mapview(shp_raster[855,])

mp <- mapview(slopeMap)


mp + mp_shp





