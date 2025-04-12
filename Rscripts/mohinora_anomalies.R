# --- GASTIS 2024

# --- Elaborado Mar 22, 2024
# --- PRELIM version in mohinora_anomalies.R para 
# --- el Diplomado en Geomática 2023

# --- En este script presentamos un análisis de anomalías.
# --- Usamos código en paralelo para eficientar el cómputo 

# --- DATASET: NDVI MOD13Q1 v061 en Cerro Mohinora, Chihuahua, 2000-2024

# --- ADDicionalmente, este script requiere archivos
# --- MOD13Q1_061_250m_16_days_NDVI_interpol.tif
# --- creado con el archivo mohinora_interpolation.R

# --- Preámbulo
library(terra)
library(rasterVis)
library(mapview)
library(RColorBrewer)
library(gtools)
library(foreach)
library(doParallel)
library(geoTS)
library(sf)
# library(raster)

source( "Rscripts/auxFUN.R" )

# ---

SHPfiles <- list.files(path = paste0( getwd(), "/data/outputs" ),
                       pattern = ".shp$",
                       full.names = TRUE)

mohinora_shp <- read_sf(SHPfiles[1])


FILES_NDVI_imputation <- list.files(path = paste0( getwd(), "/data/outputs/mohinora_imputation" ),
                                    pattern = ".tif$",
                                    full.names = TRUE)

FILES_NDVI_interpolation <- list.files(path = paste0( getwd(), "/data/outputs/mohinora_interpolation" ),
                                       pattern = ".tif",
                                       full.names = TRUE)

FILES_NDVI <- c(FILES_NDVI_imputation, FILES_NDVI_interpolation)

mohinora_DATA <- rast(FILES_NDVI) # raster::stack(FILES_NDVI)

mohinora_DATA_rTp <- spRast_valuesCoords(mohinora_DATA)

# -----------------------------------------------
# --- Análisis de cambio abrupto: anomalías --- #
# -----------------------------------------------

plot(subset(mohinora_DATA,453))
lines( mohinora_SHP, lwd=4 )

# --- Para analizar la serie de tiempo de cualquier píxel 
# --- en la imagen sigue estos pasos:
# --- 1. Ejecuta la línea de abajo
XY <- locator()
# --- 2. Haz click (SOLO UNA VEZ) en algún píxel en la imagen
# --- 3. Presiona la tecla ESC de tu teclado
# --- 4. Continúa con el script a partir de la línea 70

pixel_full <- get_timeSeries_byClicking(c(XY$x, XY$y),
                                        df=mohinora_DATA_rTp$coords)

pixel_mat <- get_pixel_matrix(mohinora_DATA_rTp$values[pixel_full$coord,])

# --- promedio por fecha de adquisicion
mu <- apply(pixel_mat, 2, mean)

# --- desviacion estándar por fecha de adquisicion
sigma <- apply(pixel_mat, 2, sd)

# --- anomalias estandarizadas
anomalia <- (pixel_mat - mu)/sigma

# --- Recordando la densidad normal
# --- y el cálculo de algunas probabilidades

x <- seq(-5,5, by = 0.05)
y <- dnorm(x = x)

# --- APROX
lB <- c(-1,-2,-3,-4) # lB:= lower bound, límite inferior o límite izquierdo
uB <- c(1,2,3,4) # uB:= upper bouund, límite superior o límite derecho

# --- EXACT
# lB <- c(-qnorm(0.85),
#         -qnorm(0.975),
#         -qnorm(0.999),
#         -qnorm(0.999975)) # lB:= lower bound, límite inferior o límite izquierdo
# 
# uB <- c(qnorm(0.85),
#         qnorm(0.975),
#         qnorm(0.999),
#         qnorm(0.999975)) # uB:= upper bouund, límite superior o límite derecho

# --- 1sigma, aprox 85% (qnorm(0.85);) de 
# --- la probabilidad total
# --- está concentrada en esta región
xz1 <- x
xz1[x >= uB[1]] <- NA
xz1[x <= lB[1]] <- NA

yz1 <- y
yz1[x >= uB[1]] <- NA
yz1[x <= lB[1]] <- NA

# --- 2 sigma, aprox 97.5% (qnorm(0.975);) de 
# --- la probabilidad total
# --- está concentrada en esta región

xz2 <- x
xz2[x >= uB[2]] <- NA
xz2[x <= lB[2]] <- NA

yz2 <- y
yz2[x >= uB[2]] <- NA
yz2[x <= lB[2]] <- NA

# --- 3 sigma, aprox 99.9% (qnorm(0.999);) de 
# --- la probabilidad total
# --- está concentrada en esta región

xz3 <- x
xz3[x >= uB[3]] <- NA
xz3[x <= lB[3]] <- NA

yz3 <- y
yz3[x >= uB[3]] <- NA
yz3[x <= lB[3]] <- NA

# --- 4 sigma, aprox 99.9975% (qnorm(0.999975);) de 
# --- la probabilidad total
# --- está concentrada en esta región

xz4 <- x
xz4[x >= uB[4]] <- NA
xz4[x <= lB[4]] <- NA

yz4 <- y
yz4[x >= uB[4]] <- NA
yz4[x <= lB[4]] <- NA

# ----
# --- PRESTAR atención a los colores

myPal <- brewer.pal('RdYlGn', n=8)

yRan <- range(y,yz1,yz2,yz3,yz4,na.rm = T)

plot(x, y, type = "l", col = "red", lwd = 3,
     ylab = "", xlab = "", main = "",
     ylim=yRan)

# --- 1sigma region
a <- x[!is.na(yz1)]
b <- yz1[!is.na(yz1)]

polygon(x=c(a[21:39],rev(a[21:39])),
        y=c(b[21:39],rep(0,19)),
        col = myPal[5], border = NA)

polygon(x=c(a[1:21],rev(a[1:21])),
        y=c(b[1:21],rep(0,21)),
        col = myPal[4], border = NA)

abline(v=lB[1], col=myPal[4], lwd=2)
abline(v=uB[1], col=myPal[5], lwd=2)

# --- 2sigma region
a2 <- x[!is.na(yz2)]
b2 <- yz2[!is.na(yz2)]

polygon(x=c(a2[60:79],rev(a2[60:79])),
        y=c(b2[60:79],rep(0,20)),
        col = myPal[6], border = NA)

polygon(x=c(a2[1:20],rev(a2[1:20])),
        y=c(b2[1:20],rep(0,20)),
        col = myPal[3], border = NA)

abline(v=lB[2], col=myPal[3], lwd=2)
abline(v=uB[2], col=myPal[6], lwd=2)

# --- 3sigma region
a3 <- x[!is.na(yz3)]
b3 <- yz3[!is.na(yz3)]


polygon(x=c(a3[101:119],rev(a3[101:119])),
        y=c(b3[101:119],rep(0,19)),
        col = myPal[7], border = NA)

polygon(x=c(a3[1:22],rev(a3[1:22])),
        y=c(b3[1:22],rep(0,22)),
        col = myPal[2], border = NA)

abline(v=lB[3], col=myPal[2], lwd=2)
abline(v=uB[3], col=myPal[7], lwd=2)

# --- 4sigma region
abline(v=lB[4], col=myPal[1], lwd=2)
abline(v=uB[4], col=myPal[8], lwd=2)

# --- obj anomalia como un obj 'ts'
anomalia_ts <- ts(c(t(anomalia)), 
                  start = c(2000,1), 
                  end = c(2024,23),
                  frequency = 23)

# --- PRESTAR atención a los colores
# --- ligar las líneas de abajo con las
# --- regiones de probabilidades definidas a partir
# --- de la densidad normal

plot(anomalia_ts, ylab="")
abline(h=1, col=myPal[5], lwd=3) # una desviación estándar
abline(h=-1, col=myPal[4], lwd=3) 

abline(h=2, col=myPal[6], lwd=3) # 2 desviaciones estándar
abline(h=-2, col=myPal[3], lwd=3)

abline(h=3, col=myPal[7], lwd=3) # 3 desviaciones estándar
abline(h=-3, col=myPal[2], lwd=3)

abline(h=4, col=myPal[8], lwd=3) # 4 desviaciones estándar
abline(h=-4, col=myPal[1], lwd=3)

# ---

startYear <- 2004
endYear <- 2005
Title <- paste0("pixel de ", startYear, " a ", endYear)

plot(anomalia_ts, ylab="", 
     xlim=c(startYear, endYear),
     main=Title)
abline(h=1, col=myPal[5], lwd=3) # una desviación estándar
abline(h=-1, col=myPal[4], lwd=3) 

abline(h=2, col=myPal[6], lwd=3) # 2 desviaciones estándar
abline(h=-2, col=myPal[3], lwd=3)

abline(h=3, col=myPal[7], lwd=3) # 3 desviaciones estándar
abline(h=-3, col=myPal[2], lwd=3)

abline(h=4, col=myPal[8], lwd=3) # 4 desviaciones estándar
abline(h=-4, col=myPal[1], lwd=3)


# ------------------------------------------
# --- Anomalías: cómputo en paralelo --- #
# ------------------------------------------

# --- CODIGO EN PARALELO

df_anomalias2000 <- matrix(nrow=nrow(mohinora_DATA_rTp$values), 
                           ncol=25) # primeras 2 columnas x, y, restantes 23 los valores de las anomalias
df_anomalias2000[,1:2] <- mohinora_DATA_rTp$coords[,1:2]

df_anomalias2004 <- matrix(nrow=nrow(mohinora_DATA_rTp$values), 
                           ncol=25)
df_anomalias2004[,1:2] <- mohinora_DATA_rTp$coords[,1:2]

df_anomalias2008 <- matrix(nrow=nrow(mohinora_DATA_rTp$values), 
                           ncol=25)
df_anomalias2008[,1:2] <- mohinora_DATA_rTp$coords[,1:2]

numCores <- detectCores()

# --- Asegurarse de crear /RData/progressReports/mohinora
progressReportFile <- paste0(getwd(), "/RData/progressReports/mohinora/progress_anomalies.txt" )
file.create(path=progressReportFile, showWarnings=FALSE)

write("===ANOMALIES computation began at===",
      file=progressReportFile, append=TRUE)
write(as.character(Sys.time()[1]), file=progressReportFile,
      append=TRUE)

kluster <- parallel::makeCluster(numCores-1, outfile="")
registerDoParallel(kluster)


output <- foreach(i=1:nrow(mohinora_DATA_rTp$values), .combine="rbind") %dopar% { 
  
  pixel <- mohinora_DATA_rTp$values[i, ]
  
  pixel_mat <- get_pixel_matrix(pixel)
  
  mu <- apply(pixel_mat, 2, mean, na.rm=TRUE)
  
  sigma <- apply(pixel_mat, 2, sd, na.rm=TRUE)
  
  anomalia <- (pixel_mat - mu) / sigma
  
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

DIR_anomalies <- paste0( getwd(), "/RData/mohinora_anomalies" )

dir.create(DIR_anomalies, recursive = TRUE)

save(df_anomalias2000, file=paste0(DIR_anomalies, "/mohinora_anomalies/2000.RData"))
save(df_anomalias2004, file=paste0(DIR_anomalies, "/mohinora_anomalies/2004.RData"))
save(df_anomalias2008, file=paste0(DIR_anomalies, "/mohinora_anomalies/2008.RData"))
# ---

# -----------------------
# --- RASTERIZACION --- #
# -----------------------

PROJECTION <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs"

# --- descomentar las siguientes 3 lineas
# --- si no tienes creados los objetos df_anomalias2000,
# --- df_anomalias2004 y df_anomalias2008
# df_anomalias2000 <- LoadToEnvironment(paste0(DIR_anomalies, "/mohinora_anomalies/2000.RData"))$df_anomalias2000
# df_anomalias2004 <- LoadToEnvironment(paste0(DIR_anomalies, "/mohinora_anomalies/2004.RData"))$df_anomalias2004
# df_anomalias2008 <- LoadToEnvironment(paste0(DIR_anomalies, "/mohinora_anomalies/2008.RData"))$df_anomalias2008

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


# --- Asegurarse de haber creado el folder /outputs/mohinora_anomalies

DIR_outpus <- paste0( getwd(), "/data/outputs/mohinora_anomalies" )
dir.create(DIR_outpus, recursive = TRUE)

writeRaster(map_anomalies2000,
            filename = paste0( DIR_outpus, "/2000" ),
            format="GTiff", datatype="FLT4S", overwrite=TRUE)

writeRaster(map_anomalies2004,
            filename = paste0( DIR_outpus, "/2004" ),
            format="GTiff", datatype="FLT4S", overwrite=TRUE)

writeRaster(map_anomalies2008,
            filename = paste0( DIR_outpus, "/2008" ),
            format="GTiff", datatype="FLT4S", overwrite=TRUE)

# -----------------------
# --- VISUALIZACION --- #
# -----------------------
# --- PRESTAR atención a los colores
# --- interpretación de las áreas coloreadas
# --- es similar a la mostrada con las regiones
# --- de probabilidad normal

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
