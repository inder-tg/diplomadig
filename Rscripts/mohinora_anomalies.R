# --- Diplomado Geomática, IG, UNAM, 2022
# --- Módulo VI Análisis de series de tiempo de imágenes satelitales con R

# --- Elaborado: Junio 2, 2022
# --- Impartido: Junio 25, 2022
# --- Modificado: Agosto 3, 2023

# --- En este script presentamos un análisis de anomalías.
# --- Usamos código en paralelo para eficientar el cómputo a alta escala.

# --- DATASET: NDVI MOD13Q1 en Cerro Mohinora, Chihuahua
# --- ADDicionalmente, este script requiere archivos
# --- MOD13Q1_NDVI_2000001.tif
# --- MOD13Q1_NDVI_2000017.tif
# --- MOD13Q1_NDVI_2000033.tif
# --- creados a través del archivo mohinora_imputation.R

# --- Preámbulo
library(raster)
library(rasterVis)
library(mapview)
library(RColorBrewer)
library(gtools)
library(foreach)
library(doParallel)

source( paste0( getwd(), "/Rscripts/auxFUN.R" ) )

# ---

listFILES_mohinora <- list.files(path=paste0(getwd(), "/data/mohinora"), 
                                 pattern=".tif", 
                                 full.names=TRUE)

stack_primeras3Imagenes <- stack(listFILES_mohinora[2:4])

mohinora_NDVI_rTp <- LoadToEnvironment( paste0( getwd(), "/RData/mohinora_imputation/NDVI_MOD13Q1_00_09_lineal.RData" ) )$mohinora_NDVI_rTp

primeras3_NDVI_rTp <- rasterToPoints(stack_primeras3Imagenes)

mohinora_NDVI_rTp_full <- matrix(nrow = nrow(mohinora_NDVI_rTp), ncol=232)

mohinora_NDVI_rTp_full[,1] <- primeras3_NDVI_rTp[,1]

mohinora_NDVI_rTp_full[,2] <- primeras3_NDVI_rTp[,2]

mohinora_NDVI_rTp_full[,3:5] <- primeras3_NDVI_rTp[,3:5]

mohinora_NDVI_rTp_full[,6:232] <- mohinora_NDVI_rTp[,3:229]

# ---

# SHP_anp <- list.files( path = paste0( getwd(), "/data/anp_2021" ),
#                        pattern = ".shp", 
#                        full.names = TRUE)
# 
# shp_anp <- shapefile( SHP_anp[1] )
# 
# shp_anp_sinu <- spTransform(shp_anp, crs(stack_primeras3Imagenes))

# ---


# -----------------------------------------------
# --- Análisis de cambio abrupto: anomalías --- #
# -----------------------------------------------

XY <- list(x=-10698785, y=2893289)

xy <- get_timeSeries_byClicking(c(XY$x, XY$y),
                                df=mohinora_NDVI_rTp_full)

pixel_full <- mohinora_NDVI_rTp_full[xy$coord, 3:ncol(mohinora_NDVI_rTp_full)]

pixel_mat <- get_pixel_matrix(pixel_full * 1e-4)

mu <- apply(pixel_mat, 2, mean)

sigma <- apply(pixel_mat, 2, sd)

anomalia <- (pixel_mat - mu)/sigma

anomalia_ts <- ts(c(t(anomalia)), start = c(2000,1), end = c(2009,23),
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

df_anomalias2000 <- matrix(nrow=nrow(mohinora_NDVI_rTp_full), ncol=25)
df_anomalias2000[,1:2] <- mohinora_NDVI_rTp_full[,1:2]

df_anomalias2004 <- matrix(nrow=nrow(mohinora_NDVI_rTp_full), ncol=25)
df_anomalias2004[,1:2] <- mohinora_NDVI_rTp_full[,1:2]

df_anomalias2008 <- matrix(nrow=nrow(mohinora_NDVI_rTp_full), ncol=25)
df_anomalias2008[,1:2] <- mohinora_NDVI_rTp_full[,1:2]

# df_slope <- matrix(nrow=nrow(mohinora_NDVI_rTp_full), ncol=3)
# df_slope[,1:2] <- mohinora_NDVI_rTp_full[,1:2]

numCores <- detectCores()

# --- progress report file (to check out on the process)
progressReportFile <- paste0(getwd(), "/RData/progressReports/mohinora_anomalies.txt" )
file.create(path=progressReportFile, showWarnings=FALSE)

write("===ANOMALIES computation began at===",
      file=progressReportFile, append=TRUE)
write(as.character(Sys.time()[1]), file=progressReportFile,
      append=TRUE)

kluster <- parallel::makeCluster(numCores-1, outfile="")
registerDoParallel(kluster)

# nrow(mohinora_NDVI_rTp_full)

output <- foreach(i=1:nrow(mohinora_NDVI_rTp_full), .combine="rbind") %dopar% { 
                    
                    pixel <- mohinora_NDVI_rTp_full[i, 3:ncol(mohinora_NDVI_rTp_full)] * 1e-4
                    
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

PROJECTION <- raster::projection(stack_primeras3Imagenes)

df_anomalias2000 <- LoadToEnvironment(paste0(getwd(),"/RData/mohinora_anomalies/2000.RData"))$df_anomalias2000
df_anomalias2004 <- LoadToEnvironment(paste0(getwd(),"/RData/mohinora_anomalies/2004.RData"))$df_anomalias2004
df_anomalias2008 <- LoadToEnvironment(paste0(getwd(),"/RData/mohinora_anomalies/2008.RData"))$df_anomalias2008

map_anomalies2000 <- brick()
map_anomalies2004 <- brick()
map_anomalies2008 <- brick()

for( i in 1:23 ){
  temp2000 <- matrixToRaster_test(matrix=df_anomalias2000[,c(1,2,i+2)],
                                  projection=PROJECTION)
  
  temp2004 <- matrixToRaster_test(matrix=df_anomalias2004[,c(1,2,i+2)], 
                                  projection=PROJECTION)
  
  temp2008 <- matrixToRaster_test(matrix=df_anomalias2008[,c(1,2,i+2)], 
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

writeRaster(map_anomalies2000,
            filename = paste0( getwd(), "/data/mohinora_anomalies/2000" ),
            format="GTiff", datatype="FLT4S", overwrite=TRUE)

writeRaster(map_anomalies2004,
            filename = paste0( getwd(), "/data/mohinora_anomalies/2004" ),
            format="GTiff", datatype="FLT4S", overwrite=TRUE)

writeRaster(map_anomalies2008,
            filename = paste0( getwd(), "/data/mohinora_anomalies/2008" ),
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

