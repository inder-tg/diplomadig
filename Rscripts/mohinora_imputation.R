
# --- Diplomado Geomática, IG, UNAM, 2022
# --- Módulo VI Análisis de series de tiempo de imágenes satelitales con R

# --- Elaborado: Junio 2, 2022
# --- Impartido: Junio 24, 2022
# --- Modificado: Julio 25, 2023; Agosto 10, 2024

# --- En este script presentamos un ejemplo para imputar estadísticamente (rellenar)
# --- las primeras 3 fechas del DATASET con base en la curva de climatología.
# --- Usamos código en paralelo para eficientar el cómputo a alta escala

# --- DATASET: NDVI MOD13Q1 en Cerro Mohinora, Chihuahua

# --- NOTA: Este código no es completamente automático; de vez en vez se requerirá
# --- crear folders o descomentar líneas de código (por ejemplo en el uso de rutinas ligada al paquete raster),
# --- esas líneas están marcada con el texto "ACTION REQUIRED!!!"

# --- Preámbulo
library(raster)
library(terra)
library(mapview)
library(sf)
library(RColorBrewer)
library(gtools)
library(geoTS)
library(foreach)
library(doParallel)

source( paste0( getwd(), "/Rscripts/auxFUN.R" ) )

# --- 

# --- Carga de datos

listFILES_mohinora <- list.files( path = paste0( getwd(), "/data/mohinora" ), 
                              pattern = ".tif", 
                              full.names = TRUE )

# SHP_anp <- list.files( path = paste0( getwd(), "/data/anp_2021" ),
#                        pattern = ".shp",
#                        full.names = TRUE)

SHP_anp <- list.files( path = paste0( getwd(), "/data/ANP" ),
                       pattern = ".shp",
                       full.names = TRUE)


# stack_NDVI_Mohinora <- stack( listFILES_mohinora[2] ) # ACTION REQUIRED!!!
stack_NDVI_Mohinora <- rast( listFILES_mohinora[2] )

# shp_anp <- shapefile( SHP_anp[1] ) # ACTION REQUIRED!!!
shp_anp <- read_sf( SHP_anp[1] )

# ------------------------------
# --- cambio de proyección --- #
# ------------------------------

# ACTION REQUIRED!!!
# plot( subset(stack_NDVI_Mohinora, 100) )
# plot( shp_anp[144,] )
# crs(subset(stack_NDVI_Mohinora, 100))
# crs(shp_anp)


plot( subset(stack_NDVI_Mohinora, 100) )
lines( shp_anp[144,] )

subset(stack_NDVI_Mohinora, 100)
shp_anp$geometry

# ACTION REQUIRED!!!
# crs shp original es distinto a crs de stack_NDVI_Mohinora
# spTranform ayuda a realizar la reproyección
# shp_anp_sinu <- spTransform(shp_anp, raster::crs(stack_NDVI_Mohinora))

plot( subset(stack_NDVI_Mohinora, 100) )
plot( shp_anp_sinu[144,], add=TRUE, lwd=4)

# mohinora_poligono_sinu <- st_transform(x=shp_anp[144,], crs=crs(stack_NDVI_Mohinora)) # ACTION REQUIRED!!!
shp_anp_sinu <- st_transform(shp_anp[144,], st_crs(stack_NDVI_Mohinora))

# ACTION REQUIRED!!!
# Asegurarse de crear /data/shp_mohinora
write_sf(shp_anp_sinu, paste0( getwd(), "/data/shp_mohinora/mohinora.shp" ))

# ----------------------------------
# --- Datos faltantes: Exploración #
# ----------------------------------

# --- Accediendo a los numeritos

# mohinora_NDVI_rTp <- rasterToPoints(stack_NDVI_Mohinora)
mohinora_NDVI_rTp <- spRast_valueCoords(stack_NDVI_Mohinora)

plot( subset(stack_NDVI_Mohinora, 100) )
# plot( shp_anp_sinu[144,], add=TRUE, lwd=4) # ACTION REQUIRED!!!
lines(shp_anp_sinu, lwd=4)

# -----------------------------------------------------------------------------
# --- Para analizar la serie de tiempo de cualquier píxel en la imagen sigue estos
# --- pasos:
# --- 1. Ejecuta la línea de abajo
XY <- locator()
# --- 2. Haz click (SOLO UNA VEZ) en algún píxel en la imagen
# --- 3. Presiona la tecla ESC de tu teclado
# --- 4. Continúa con el script a partir de la línea 107

XY <- list(x=-10698785, y=2893289)

# xy <- get_timeSeries_byClicking(c(XY$x, XY$y),
#                                 df=mohinora_NDVI_rTp) # ACTION REQUIRED!!!

xy <- get_timeSeries_byClicking(c(XY$x, XY$y),
                                df=mohinora_NDVI_rTp$coords)

# pixel <- mohinora_NDVI_rTp[xy$coord, 3:ncol(mohinora_NDVI_rTp)] # ACTION REQUIRED!!!
pixel <- mohinora_NDVI_rTp$values[xy$coord,]

pixel_ts <- ts(as.numeric(pixel), start = c(2000,1), 
               end = c(2023,23), frequency = 23)

plot(as.numeric(pixel), main="pixel original")

plot(pixel_ts, xlab="Años", ylab="NDVI", col="darkgreen", 
     main="pixel como objeto 'ts'")

# ------------------------------------------------------------------------------
# --- CONOCE tu DATASET!!
# --- Extremo cuidado al usar ts() para definir un objeto

pixel_ts[550:552]
as.numeric(pixel[1:3])

# ------------------------------------------------------------------------------

pixel_aug <- c(NA,NA,NA, as.numeric(pixel))

pixel_aug_ts <- ts(pixel_aug, start = c(2000,1), end = c(2023,23),
                   frequency = 23)

plot(pixel_aug_ts, xlab="Años", ylab="NDVI", col="darkgreen", 
     main="pixel aumentado como objeto 'ts'")

pixel_aug_ts[550:552]
as.numeric(pixel[547:549])

# ------------------------------------------------------------------------------
# --- Uso de la curva de climatología para imputar las primeras 3 fechas del pixel

clima <- climatology(x=pixel_aug, lenPeriod=23)

boxplot(clima$matrix)

pixel_aug[1:3] <- ceiling(apply(clima$matrix[-1,1:3], MARGIN=2, FUN=median))

pixel_ts_correct <- ts(pixel_aug, start = c(2000,1), end = c(2023,23), 
                       frequency = 23)


pixel_aug[1:3]
as.numeric(pixel_ts_correct[1:3])

plot(pixel_ts_correct, ylab="NDVI", col="darkgreen", 
     main="pixel imputado como objeto 'ts'")

# pixel_output <- c(mohinora_NDVI_rTp[xy$coord,1:2], pixel_aug) # ACTION REQUIRED!!!
pixel_output <- c(mohinora_NDVI_rTp$coords[xy$coord,], pixel_aug)
# ------------------------------------------------------------------------------

# ---------------------------------------------------
# --- Datos faltantes: Imputación a gran escala --- #
# ---------------------------------------------------

# --- CODIGO EN PARALELO

# # ACTION REQUIRED!!!
# df_layer1 <- matrix(nrow=nrow(mohinora_NDVI_rTp), ncol=3)
# df_layer1[,1:2] <- mohinora_NDVI_rTp[1:2,]
# 
# df_layer2 <- matrix(nrow=nrow(mohinora_NDVI_rTp), ncol=3)
# df_layer2[,1:2] <- mohinora_NDVI_rTp[1:2,]
# 
# df_layer3 <- matrix(nrow=nrow(mohinora_NDVI_rTp), ncol=3)
# df_layer3[,1:2] <- mohinora_NDVI_rTp[1:2,]

# df_layers guardará las imputaciones
df_layer1 <- matrix(nrow=nrow(mohinora_NDVI_rTp$coords), ncol=3)
df_layer1[,1:2] <- mohinora_NDVI_rTp$coords

df_layer2 <- matrix(nrow=nrow(mohinora_NDVI_rTp$coords), ncol=3)
df_layer2[,1:2] <- mohinora_NDVI_rTp$coords

df_layer3 <- matrix(nrow=nrow(mohinora_NDVI_rTp$coords), ncol=3)
df_layer3[,1:2] <- mohinora_NDVI_rTp$coords

# ---

numCores <- detectCores()

# ACTION REQUIRED!!!
# --- Asegurarse de crear folder /RData/progressReports

progressReportFile <- paste0(getwd(), "/RData/progressReports/mohinora_imputation.txt" )
file.create(path=progressReportFile, showWarnings=FALSE)

write("===CLIMATOLOGY imputation began at===",
      file=progressReportFile, append=TRUE)
write(as.character(Sys.time()[1]), file=progressReportFile,
      append=TRUE)

kluster <- parallel::makeCluster(numCores-1, outfile="")
registerDoParallel(kluster)

output <- foreach(i=1:nrow(mohinora_NDVI_rTp$coords), .combine="rbind") %dopar% { 
  
  # pixel <- mohinora_NDVI_rTp[i, 3:ncol(mohinora_NDVI_rTp)] # ACTION REQUIRED!!!
  pixel <- mohinora_NDVI_rTp$values[i, ]
  
  clima <- climatology(x=pixel, lenPeriod=23)
  
  s <- ceiling(apply(clima$matrix[-1,1:3], MARGIN=2, FUN=median))
  
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
df_layer1[,3] <- output[,1]
df_layer2[,3] <- output[,2]
df_layer3[,3] <- output[,3]


# --- Asegurarse de crear folder /RData/mohinora_imputation
save(df_layer1, file=paste0(getwd(),"/RData/mohinora_imputation/MOD13Q1.A2000001.RData"))
save(df_layer2, file=paste0(getwd(),"/RData/mohinora_imputation/MOD13Q1.A2000017.RData"))
save(df_layer3, file=paste0(getwd(),"/RData/mohinora_imputation/MOD13Q1.A2000033.RData"))
# ---

# -----------------------
# --- RASTERIZACION --- #
# -----------------------

# ACTION REQUIRED!!!
# df_layer1 <- LoadToEnvironment(paste0(getwd(), "/RData/mohinora_imputation/MOD13Q1.A2000001.RData"))$df_layer1
# df_layer2 <- LoadToEnvironment(paste0(getwd(), "/RData/mohinora_imputation/MOD13Q1.A2000017.RData"))$df_layer2
# df_layer3 <- LoadToEnvironment(paste0(getwd(), "/RData/mohinora_imputation/MOD13Q1.A2000033.RData"))$df_layer3

PROJECTION <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs"

layer1 <- matrixToRaster(matrix=df_layer1, projection=PROJECTION)
layer2 <- matrixToRaster(matrix=df_layer2, projection=PROJECTION)
layer3 <- matrixToRaster(matrix=df_layer3, projection=PROJECTION)


# --- Guardando las imputaciones como objetos rasterLayer

# Asegurarse de haber creado el folder /TIF/mohinora
writeRaster(layer1,
            filename = paste0(getwd(), "/TIF/mohinora/MOD13Q1_NDVI_2000001.tif"),
            datatype="INT2S", overwrite=TRUE)

writeRaster(layer2,
            filename = paste0(getwd(), "/TIF/mohinora/MOD13Q1_NDVI_20000017.tif"),
            datatype="INT2S", overwrite=TRUE)

writeRaster(layer3,
            filename = paste0(getwd(), "/TIF/mohinora/MOD13Q1_NDVI_2000033.tif"),
            datatype="INT2S", overwrite=TRUE)


# -----------------------
# --- VISUALIZACION --- #
# -----------------------

# imputedTIFs <- list.files(path=paste0( getwd(), "/TIF/mohinora" ),
#                           pattern = ".tif$",
#                           full.names = TRUE)
# 
# layer1 <- rast(imputedTIFs[1])

mp <- mapview(layer1)

shp_mohinora_mp <- mapview(shp_anp_sinu)

mp + shp_mohinora_mp
