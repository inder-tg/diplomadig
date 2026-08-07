
# --- Elaborado: Feb 29, 2024
# --- Actualizado: Mar 5, 2025, Abril 4, 2025

# --- En este script presentamos un ejemplo para imputar estadísticamente (rellenar)
# --- las primeras 3 fechas del DATASET con base en la curva de climatología.
# --- Usamos código en paralelo para eficientar el cómputo a alta escala

# --- DATASET: NDVI MOD13Q1 v061 2000-2024 en Cerro Mohinora, Chihuahua

# --- Actualizado: Ago 6, 2026
# --- Se han actualizado algunas líneas de código empleando funciones actuales 
# --- para alinearnos a las expectativas del Diplomado en Geomática Edición XIX
# --- NOTA: Crear subdirectorios necesarios (ej. /data/mohinora)

# --- Preámbulo
library(terra)
library(mapview)
library(gtools)
library(geoTS)
library(foreach)
library(doParallel)
library(raster)

source("Rscripts/auxFUN.R")

# --- 

# DIR <- paste0( getwd(), "/data/mohinora" )
DIR <- list.dirs( here("data", "mohinora") )
lst_DIR <- setNames( as.list(DIR), basename(DIR) )

# --- Carga de datos

NDVIfiles <- list.files( path = lst_DIR$`250m_16_days_NDVI_QA`, # paste0( DIR, "/250m_16_days_NDVI_QA" ), 
                         pattern = ".tif",
                         full.names = TRUE )

mohinora_DATA <- rast(NDVIfiles)

SHPfiles <- list.files(path = here("data", "outputs"), # paste0( getwd(), "/data/outputs" ),
                       pattern = ".shp$",
                       full.names = TRUE)

mohinora_shp <- read_sf(SHPfiles[1])

# OJO: lst_dirs fue definido en mohinora_tmap.R
usvFILES <- list.files(path = lst_dirs$mohinora_usv7,
                       full.names = TRUE,
                       pattern = ".shp$")

mohinora_USV <- read_sf(usvFILES)

# ----------------------------------
# --- Datos faltantes: Exploración #
# ----------------------------------

# --- Accediendo a los numeritos

mohinora_DATA_rTp <- spRast_valuesCoords(mohinora_DATA) # rasterToPoints(mohinora_DATA) #spRast_valuesCoords(mohinora_DATA)

# Graficar un subset de tu objeto base
plot(subset(mohinora_DATA, 10))

# Añadir los polígonos con colores
lines(mohinora_USV, col = usv_COLORS, lwd = 6)
# plot(mohinora_USV, col = usv_COLORS, lwd = 6, add=TRUE)

# Añadir la leyenda
legend("topleft", # posición en el gráfico
       legend = usv_NAMES,         # nombres de las categorías
       col = usv_COLORS,           # colores de borde
       lwd = 6,                    # grosor de línea en la leyenda
       pt.cex = 1,                  # tamaño del texto
       bty = "n",                  # sin caja alrededor
       inset = c(0.015, 0.05))

# ACÁ ME QUEDÉ ANOCHE


# -----------------------------------------------------------------------------
# --- Para analizar la serie de tiempo de cualquier píxel en la imagen sigue estos
# --- pasos:
# --- 1. Ejecuta la línea de abajo
XY <- locator()
# --- 2. Haz click (SOLO UNA VEZ) en algún píxel en la imagen
# --- 3. Presiona la tecla ESC de tu teclado
# --- 4. Continúa con el script a partir de la línea 107

# XY <- list(x=-10698697, y=2894003)

xy <- get_timeSeries_byClicking(c(XY$x, XY$y),
                                df=mohinora_DATA_rTp$coords)

pixel <- mohinora_DATA_rTp$values[xy$coord, ]

# pixel <- mohinora_DATA_rTp$values[295, ]

pixel_ts <- ts(pixel, start = c(2000,1), end = c(2024,23),
               frequency = 23)

plot(pixel, main="pixel original")

plot(pixel_ts, xlab="Años", ylab="NDVI", col="darkgreen", 
     main="pixel como objeto 'ts'")

# ------------------------------------------------------------------------------
# --- CONOCE tu DATASET!!
# --- Extremo cuidado al usar ts() para definir un objeto

pixel_ts[573:575] # últimas 3 entradas del vector pixel_ts 
as.numeric(pixel[570:572]) # últimas 3 entradas del vector pixel
as.numeric(pixel[1:3]) # primeras 3 entradas del vector pixel
# ------------------------------------------------------------------------------

pixel_aug <- c(NA,NA,NA, as.numeric(pixel))

pixel_aug_ts <- ts(pixel_aug, start = c(2000,1), end = c(2024,23),
                   frequency = 23)

plot(pixel_aug_ts, xlab="Años", ylab="NDVI", col="darkgreen", 
     main="pixel aumentado como objeto 'ts'")

pixel_aug_ts[573:575]
as.numeric(pixel[570:572])

# -----------------------------------------------------------------------------
# --- Uso de la curva de climatología para imputar las primeras 3 fechas del pixel

clima <- climatology(x=pixel_aug, lenPeriod=23)

boxplot(clima$matrix)

pixel_aug[1:3] <- ceiling(apply(clima$matrix[-1,1:3], MARGIN=2, 
                                FUN=median))

pixel_ts_correct <- ts(pixel_aug, start = c(2000,1), end = c(2024,23), 
                       frequency = 23)

pixel_aug[1:3]
as.numeric(pixel_ts_correct[1:3])

plot(pixel_aug_ts, xlab="Años", ylab="NDVI", col="darkgreen", 
     main="pixel aumentado como objeto 'ts'")

plot(pixel_ts_correct, ylab="NDVI", col="darkgreen", 
     main="pixel imputado visto como objeto 'ts'")

pixel_output <- c(mohinora_DATA_rTp$coords[xy$coord,], pixel_aug)
# ------------------------------------------------------------------------------

# ---------------------------------------------------
# --- Datos faltantes: Imputación a gran escala --- #
# ---------------------------------------------------

# --- CODIGO EN PARALELO

# df_layers guardará las imputaciones
df_layer1 <- matrix(nrow=nrow(mohinora_DATA_rTp$values), ncol=3)
df_layer1[,1:2] <- mohinora_DATA_rTp$coords

df_layer2 <- matrix(nrow=nrow(mohinora_DATA_rTp$values), ncol=3)
df_layer2[,1:2] <- mohinora_DATA_rTp$coords

df_layer3 <- matrix(nrow=nrow(mohinora_DATA_rTp$values), ncol=3)
df_layer3[,1:2] <- mohinora_DATA_rTp$coords

# progress report file (to check out on the process)

DIR_progress <- paste0( getwd(), "/RData/progressReports" )

if( !dir.exists(DIR_progress) ){
  dir.create(DIR_RData, recursive = TRUE)
}

progressReportFile <- paste0( DIR_progress, "/progress_imputation.txt" )
file.create(path=progressReportFile, showWarnings=FALSE)

write("===CLIMATOLOGY imputation began at===",
      file=progressReportFile, append=TRUE)
write(as.character(Sys.time()[1]), file=progressReportFile,
      append=TRUE)

numCores <- detectCores()

kluster <- parallel::makeCluster(numCores-1, outfile="")
registerDoParallel(kluster)

output <- foreach(i=1:nrow(mohinora_DATA_rTp$values), .combine="rbind") %dopar% {
  
  pixel <- mohinora_DATA_rTp$values[i, ]
  
  pixel_aug <- c(NA,NA,NA, as.numeric(pixel))
  
  clima <- climatology( x=pixel_aug, lenPeriod=23 )
  
  s <- ceiling( apply( clima$matrix[-1,1:3], MARGIN=2, 
                     FUN = median, na.rm = TRUE ) )
  
  if(i %% 100 ==0){
    texto <- paste0("Working on ROW: ", i)
    write(texto, file=progressReportFile, append=TRUE)
  }
  
  return(s)
}
stopCluster(kluster)

write( as.character(Sys.time()[1]), file=progressReportFile, append=TRUE)
write( "===CLIMATOLOGY imputation ended here===", 
       file=progressReportFile, append=TRUE)
# ---

# --- Guardando las imputaciones como objetos matrix
df_layer1[,3] <- output[,1]
df_layer2[,3] <- output[,2]
df_layer3[,3] <- output[,3]

# --- asegurarse de crear /RData/mohinora_imputation

DIR_imputation <- paste0( getwd(), "/RData/mohinora_imputation" )

dir.create(DIR_imputation, recursive = TRUE)

save(df_layer1, file=paste0(DIR_imputation, "/MOD13Q1.A2000001.RData"))
save(df_layer2, file=paste0(DIR_imputation, "/MOD13Q1.A2000017.RData"))
save(df_layer3, file=paste0(DIR_imputation, "/MOD13Q1.A2000033.RData"))

# ---

# -----------------------
# --- RASTERIZACION --- #
# -----------------------

# --- usar las siguientes 3 líneas si se ha empezado una nueva sesión
# --- de trabajo

# df_layer1 <- LoadToEnvironment(paste0(DIR_RData, "/mohinora_imputation/MOD13Q1.A2000001.RData"))$df_layer1
# df_layer2 <- LoadToEnvironment(paste0(DIR_RData, "/mohinora_imputation/MOD13Q1.A2000017.RData"))$df_layer2
# df_layer3 <- LoadToEnvironment(paste0(DIR_RData, "/mohinora_imputation/MOD13Q1.A2000033.RData"))$df_layer3

PROJECTION <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs"

layer1 <- matrixToRaster(matrix=df_layer1, projection=PROJECTION)
layer2 <- matrixToRaster(matrix=df_layer2, projection=PROJECTION)
layer3 <- matrixToRaster(matrix=df_layer3, projection=PROJECTION)

# --- Asegurarse de crear /data/outputs/mohinora_imputation
# --- Guardando las imputaciones en archivos GeoTiff

DIR_outputs <- paste0( getwd(), "/data/outputs/mohinora_imputation" )

dir.create(DIR_outputs, recursive = TRUE)

baseName <- "h08v06.061.2024022950859.250m_16_days_NDVI.tif"

raster::writeRaster(layer1,
                    filename = paste0(DIR_outputs, "/MOD13Q1.A2000001.",
                                      baseName),
                    datatype="INT2S", overwrite=TRUE)

raster::writeRaster(layer2,
                    filename = paste0(DIR_outputs, "/MOD13Q1.A2000017.",
                                      baseName),
                    datatype="INT2S", overwrite=TRUE)

raster::writeRaster(layer3,
                    filename = paste0(DIR_outputs, "/MOD13Q1.A2000033.",
                                      baseName),
                    datatype="INT2S", overwrite=TRUE)

# -----------------------
# --- VISUALIZACION --- #
# -----------------------

mp <- mapview(layer1)

shp_mohinora_mp <- mapview(mohinora_shp)

mp + shp_mohinora_mp
