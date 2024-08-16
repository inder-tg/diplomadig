# --- Diplomado Geomática, IG, UNAM, 2022
# --- Módulo VI Análisis de series de tiempo de imágenes satelitales con R

# --- Elaborado: Junio 2, 2022
# --- Impartido: Junio 25, 2022
# --- Modificado: Agosto 3, 2023, Agosto 10, 2024

# --- En este script presentamos un análisis de anomalías.
# --- Usamos código en paralelo para eficientar el cómputo a alta escala.

# --- DATASET: NDVI MOD13Q1 v061 en Cerro Mohinora, Chihuahua, 2000-2023 

# --- ADDicionalmente, este script requiere archivos
# --- MOD13Q1_061_250m_16_days_NDVI_interpol.tif

# --- NOTA: Este código no es completamente automático; de vez en vez se requerirá
# --- crear folders o descomentar líneas de código (por ejemplo en el uso de rutinas ligada al paquete raster),
# --- esas líneas están marcada con el texto "ACTION REQUIRED!!!"

# --- Preámbulo
library(raster)
library(terra)
library(rasterVis)
library(mapview)
library(RColorBrewer)
library(gtools)
library(foreach)
library(doParallel)
library(geoTS)

source( paste0( getwd(), "/Rscripts/auxFUN.R" ) )

# ---

listFILES_mohinora <- list.files(path=paste0(getwd(), "/data/mohinora"), 
                                 pattern=".tif$", 
                                 full.names=TRUE)

shpFILES_mohinora <- list.files(path=paste0(getwd(), "/data/shp_mohinora"), 
                                pattern=".shp$", 
                                full.names=TRUE)

# stack_NDVI_Mohinora <- stack( listFILES_mohinora[1] ) # ACTION REQUIRED!!!
stack_NDVI_Mohinora <- rast( listFILES_mohinora[1] )

# shp_mohinora <- shapefile( shpFILES_mohinora[1] ) # ACTION REQUIRED!!!
shp_mohinora <- read_sf( shpFILES_mohinora[1] )

# mohinora_NDVI_rTp_full <- rasterToPoints(stack_NDVI_Mohinora) # ACTION REQUIRED!!!
mohinora_NDVI_rTp_full <- spRast_valueCoords(stack_NDVI_Mohinora)

# ---

# -----------------------------------------------
# --- Análisis de cambio abrupto: anomalías --- #
# -----------------------------------------------

# XY <- list(x=-10698785, y=2893289)

# xy <- get_timeSeries_byClicking(c(XY$x, XY$y),
#                                 df=mohinora_NDVI_rTp_full) # ACTION REQUIRED!!!
xy <- get_timeSeries_byClicking(c(XY$x, XY$y),
                                df=mohinora_NDVI_rTp_full$coords)

# pixel_full <- mohinora_NDVI_rTp_full[xy$coord, 3:ncol(mohinora_NDVI_rTp_full)] # ACTION REQUIRED!!!
pixel_full <- mohinora_NDVI_rTp_full$values[xy$coord, ]

pixel_mat <- get_pixel_matrix(pixel_full * 1e-4)

# --- promedio por fecha de adquisicion
mu <- apply(pixel_mat, 2, mean)

# --- desviacion estándar por fecha de adquisicion
sigma <- apply(pixel_mat, 2, sd)


# --- anomalias estandarizadas
anomalia <- (pixel_mat - mu)/sigma

anomalia_ts <- ts(c(t(anomalia)), start = c(2000,1), end = c(2023,23),
                  frequency = 23)

plot(anomalia_ts, ylab="")
abline(h=1, col="green")
abline(h=-1, col="green")

abline(h=2, col="darkgreen")
abline(h=-2, col="darkgreen")

abline(h=3, col="darkolivegreen")
abline(h=-3, col="darkolivegreen")

abline(h=4, col="brown")
abline(h=-4, col="brown")

# ---
startYear <- 2004
endYear <- 2005
Title <- paste0("pixel de ", startYear, " a ", endYear)

plot(anomalia_ts, ylab="", xlim=c(startYear, endYear),
     main=Title)
abline(h=1, col="green")
abline(h=-1, col="green")

abline(h=2, col="darkgreen")
abline(h=-2, col="darkgreen")

abline(h=3, col="darkolivegreen")
abline(h=-3, col="darkolivegreen")

abline(h=4, col="brown")
abline(h=-4, col="brown")


# ------------------------------------------
# --- Anomalías: cómputo a gran escala --- #
# ------------------------------------------

# --- CODIGO EN PARALELO

# # ACTION REQUIRED!!!
# df_anomalias2000 <- matrix(nrow=nrow(mohinora_NDVI_rTp_full), ncol=25)
# df_anomalias2000[,1:2] <- mohinora_NDVI_rTp_full[,1:2]
# 
# df_anomalias2004 <- matrix(nrow=nrow(mohinora_NDVI_rTp_full), ncol=25)
# df_anomalias2004[,1:2] <- mohinora_NDVI_rTp_full[,1:2]
# 
# df_anomalias2008 <- matrix(nrow=nrow(mohinora_NDVI_rTp_full), ncol=25)
# df_anomalias2008[,1:2] <- mohinora_NDVI_rTp_full[,1:2]

df_anomalias2000 <- matrix(nrow=nrow(mohinora_NDVI_rTp_full$coords), ncol=25)
df_anomalias2000[,1:2] <- mohinora_NDVI_rTp_full$coords

df_anomalias2004 <- matrix(nrow=nrow(mohinora_NDVI_rTp_full$coords), ncol=25)
df_anomalias2004[,1:2] <- mohinora_NDVI_rTp_full$coords

df_anomalias2008 <- matrix(nrow=nrow(mohinora_NDVI_rTp_full$coords), ncol=25)
df_anomalias2008[,1:2] <- mohinora_NDVI_rTp_full$coords


numCores <- detectCores()

# ACTION REQUIRED!!!
# --- Asegurarse de crear /RData/progressReports
progressReportFile <- paste0(getwd(), "/RData/progressReports/mohinora_anomalies.txt" )
file.create(path=progressReportFile, showWarnings=FALSE)

write("===ANOMALIES computation began at===",
      file=progressReportFile, append=TRUE)
write(as.character(Sys.time()[1]), file=progressReportFile,
      append=TRUE)

kluster <- parallel::makeCluster(numCores-1, outfile="")
registerDoParallel(kluster)

output <- foreach(i=1:nrow(mohinora_NDVI_rTp_full$coords), .combine="rbind") %dopar% { 
                    
                    # pixel <- mohinora_NDVI_rTp_full[i, 3:ncol(mohinora_NDVI_rTp_full)] * 1e-4 # ACTION REQUIRED!!!
  
                    pixel <- mohinora_NDVI_rTp_full$values[i,] * 1e-4
                    
                    pixel_mat <- get_pixel_matrix(pixel)
                    
                    mu <- apply(pixel_mat, 2, mean)
                    
                    sigma <- apply(pixel_mat, 2, sd)
                    
                    anomalia <- (pixel_mat - mu)/sigma
                    
                    anomalia_2000 <- anomalia[1,]
                    
                    anomalia_2004 <- anomalia[5,]
                    
                    anomalia_2008 <- anomalia[9,]
                    
                    s <- c(anomalia_2000, anomalia_2004, anomalia_2008)
                    
                    if(i %% 100 ==0){
                      texto <- paste0("Working on ROW: ", i)
                      write(texto, file=progressReportFile, append=TRUE)
                    }
                    
                    return(s)
                  }
stopCluster(kluster)

write( as.character(Sys.time()[1]), file=progressReportFile, append=TRUE)
write( "===ANOMALIES analysis ended here===", file=progressReportFile, append=TRUE)
# ---

# --- Guardando las anomalías como objetos matrix
df_anomalias2000[,3:25] <- output[,1:23]
df_anomalias2004[,3:25] <- output[,24:46]
df_anomalias2008[,3:25] <- output[,47:69]

# --- Asegurarse de haber creado el folder /RData/mohinora_anomalies
save(df_anomalias2000, file=paste0(getwd(),"/RData/mohinora_anomalies/2000.RData"))
save(df_anomalias2004, file=paste0(getwd(),"/RData/mohinora_anomalies/2004.RData"))
save(df_anomalias2008, file=paste0(getwd(),"/RData/mohinora_anomalies/2008.RData"))
# ---

# -----------------------
# --- RASTERIZACION --- #
# -----------------------

PROJECTION <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs"

# df_anomalias2000 <- LoadToEnvironment(paste0(getwd(),"/RData/mohinora_anomalies/2000.RData"))$df_anomalias2000
# df_anomalias2004 <- LoadToEnvironment(paste0(getwd(),"/RData/mohinora_anomalies/2004.RData"))$df_anomalias2004
# df_anomalias2008 <- LoadToEnvironment(paste0(getwd(),"/RData/mohinora_anomalies/2008.RData"))$df_anomalias2008

map_anomalies2000 <- brick()
map_anomalies2004 <- brick()
map_anomalies2008 <- brick()

for( i in 1:23 ){
  temp2000 <- matrixToRaster(matrix=df_anomalias2000[,c(1,2,i+2)],
                                  projection=PROJECTION)
  
  temp2004 <- matrixToRaster(matrix=df_anomalias2004[,c(1,2,i+2)], 
                                  projection=PROJECTION)
  
  temp2008 <- matrixToRaster(matrix=df_anomalias2008[,c(1,2,i+2)], 
                                  projection=PROJECTION)
  
  map_anomalies2000 <- addLayer(map_anomalies2000, temp2000)
  
  map_anomalies2004 <- addLayer(map_anomalies2004, temp2004)

  map_anomalies2008 <- addLayer(map_anomalies2008, temp2008)
}

DoY <- seq(1, 365, by=16)
LABELS <- sapply(1:length(DoY), function(s) paste0("DoY-", DoY[s]) )

names(map_anomalies2000) <- LABELS

names(map_anomalies2004) <- LABELS

names(map_anomalies2008) <- LABELS

# #ACTION REQUIRED!!!
# Asegurarse de crear /TIF/mohinora_anomalies
raster::writeRaster(map_anomalies2000,
                    filename = paste0( getwd(), "/TIF/mohinora_anomalies/2000" ),
                    format="GTiff", datatype="FLT4S", overwrite=TRUE)

raster::writeRaster(map_anomalies2004,
                    filename = paste0( getwd(), "/TIF/mohinora_anomalies/2004" ),
                    format="GTiff", datatype="FLT4S", overwrite=TRUE)

raster::writeRaster(map_anomalies2008,
                    filename = paste0( getwd(), "/TIF/mohinora_anomalies/2008" ),
                    format="GTiff", datatype="FLT4S", overwrite=TRUE)

# -----------------------
# --- VISUALIZACION --- #
# -----------------------

myPal <- brewer.pal('RdYlGn', n=7)
myTheme <- rasterTheme(region = myPal)

levelplot(map_anomalies2000, main="Anomalías 2000", 
          par.settings = myTheme,
          at = seq(-6,6,by=1))

levelplot(map_anomalies2004, main="Anomalías 2004", 
          par.settings = myTheme,
          at = seq(-6,6,by=1))

levelplot(map_anomalies2008, main="Anomalías 2008",
          par.settings = myTheme,
          at = seq(-6,6,by=1))
# ---

