# --- Diplomado Geomática, IG, UNAM, 2022
# --- Módulo VI Análisis de series de tiempo de imágenes satelitales con R

# --- Elaborado: Junio 2, 2022
# --- Impartido: Julio 1, 2022
# --- Actualizado: Agosto 4, 2023; Agosto 16, 2024; Abril 12, 2025

# --- Más detalles técnicos sobre este análisis se pueden encontrar en 
# --- "Clasificación de tendencias de NDVI en la península de Yucatán,
# --- México, de 2014 a 2020" de Tecuapetla et al. 2022, documento compartido en 
# --- la carpeta PDF de este proyecto.

# --- DATASET: NDVI MOD13Q1 v061 en Cerro Mohinora, Chihuahua, 2000-2024 

# --- NOTA: Este código no es completamente automático; de vez en vez se requerirá
# --- crear folders o descomentar líneas de código (por ejemplo en el uso de rutinas ligada al paquete raster),
# --- esas líneas están marcada con el texto "ACTION REQUIRED!!!"

# --- Preámbulo
library(terra)
library(raster)
library(mapview)
library(RColorBrewer)
library(gtools)
library(foreach)
library(doParallel)
library(bfast)
library(tmap)
library(geoTS)
library(sf)
library(tidyverse)
library(kableExtra)

source( "Rscripts/auxFUN.R" )

# --- DATA Loading

SHPfiles <- list.files(path = paste0( getwd(), "/data/outputs" ),
                       pattern = ".shp$",
                       full.names = TRUE)

mohinora_shp <- read_sf(SHPfiles[1])

SHPfiles_USV7 <- list.files(path = paste0( getwd(), "/data/mohinora_usv7" ),
                            pattern = ".shp$",
                            full.names = TRUE) # cambiar 7 por 5

mohinora_shp_usv7 <- read_sf(SHPfiles_USV7[1])

FILES_NDVI_imputation <- list.files(path = paste0( getwd(), "/data/outputs/mohinora_imputation" ), 
                                    pattern = ".tif$",
                                    full.names = TRUE)

FILES_NDVI_interpol <- list.files(path = paste0( getwd(), "/data/outputs/mohinora_interpolation" ),
                                  pattern = ".tif$",
                                  full.names = TRUE)

FILES_NDVI <- c(FILES_NDVI_imputation, FILES_NDVI_interpol)

mohinora_DATA <- rast(FILES_NDVI) # raster::stack(FILES_NDVI)

mohinora_DATA_rTp <- spRast_valuesCoords(mohinora_DATA) # rasterToPoints(mohinora_DATA)

# -------------------------------------------------
# --- Análisis de cambio abrupto: Exploración --- #
# -------------------------------------------------

plot(subset(mohinora_DATA,453))
lines( mohinora_shp, lwd=4 )

# --- Para analizar la serie de tiempo de cualquier píxel en la imagen sigue estos
# --- pasos:
# --- 1. Ejecuta la línea de abajo
XY <- locator()
# --- 2. Haz click (SOLO UNA VEZ) en algún píxel en la imagen
# --- 3. Presiona la tecla ESC de tu teclado
# --- 4. Continúa con el script a partir de la línea 76

xy <- get_timeSeries_byClicking(c(XY$x, XY$y), 
                                df=mohinora_DATA_rTp$coords)

# # ACTION REQUIRED!!! (SOLO PARA EL INSTRUCTOR)
# --- Clase de tendencia: 5; Año: 2002

# XY <- list(x=-10698785, y=2893289)
# 
# xy <- get_timeSeries_byClicking(c(XY$x, XY$y),
#                                 df=mohinora_DATA_interpol_rTp$coords)

# ---- otro buen ejemplo en coord 2564

pixel <- mohinora_DATA_rTp$values[xy$coord, ] * 1e-4

# --- objeto ts
pixel_ts <- ts(pixel, 
               start = c(2000,1), 
               end = c(2024,23),
               frequency = 23)
# ---

pixel_bfast01 <- bfast01( data=pixel_ts )

plot(pixel_bfast01)

pixel_bfast01$breakpoints

bfast01classify(pixel_bfast01)

getYear(start=2000, end=2024, bp=pixel_bfast01$breakpoints, freq=23)

# --- TYPES
# --- 1 - monotonic increase
# --- 2 - monotonic decrease
# --- 3 - monotonic increase with positive break
# --- 4 - monotonic decrease with negative break
# --- 5 - interruption: increase with negative break
# --- 6 - interruption: decrease with positive break
# --- 7 - reversal: increase to decrease
# --- 8 - reversal: decrease to increase

# ---

# --- Clase de tendencia: 6; Año: 2006

XY <- list(x=--10698418, y=2887497)

xy <- get_timeSeries_byClicking(c(XY$x, XY$y),
                                df=mohinora_DATA_rTp$coords)

pixel <- mohinora_DATA_rTp$values[xy$coord, ] * 1e-4

pixel_ts <- ts(as.numeric(pixel), 
               start = c(2000,1), end = c(2024, 23),
               frequency = 23)

pixel_bfast01 <- bfast01(pixel_ts, bandwidth = 0.2)

plot(pixel_bfast01)

pixel_bfast01$breakpoints

bfast01classify(pixel_bfast01)

getYear(start=2000, end=2024, bp=pixel_bfast01$breakpoints, freq=23)

# --- Clase de tendencia: 8

XY <- list(x=-10687435, y=2894415)

xy <- get_timeSeries_byClicking(c(XY$x, XY$y),
                                df=mohinora_DATA_rTp$coords)

pixel <- mohinora_DATA_rTp$values[xy$coord, ] * 1e-4

pixel_ts <- ts(as.numeric(pixel), 
               start = c(2000,1), end = c(2024, 23),
               frequency = 23)

pixel_bfast01 <- bfast01(pixel_ts, bandwidth = 0.2)

plot(pixel_bfast01)

pixel_bfast01$breakpoints

bfast01classify(pixel_bfast01)

getYear(start=2000, end=2024, bp=pixel_bfast01$breakpoints, freq=23)

# -----------------------------------------------
# --- Cambio abrupto: cómputo a gran escala --- #
# -----------------------------------------------

# --- CODIGO EN PARALELO

# --- Vamos a guardar TYPE, SIGN, STABLE, YEARS, CP

TYPE <- matrix(nrow=nrow(mohinora_DATA_rTp$coords), ncol=3)
TYPE[,1:2] <- mohinora_DATA_rTp$coords[,1:2]

SIGN <- matrix(nrow=nrow(mohinora_DATA_rTp$coords), ncol=3)
SIGN[,1:2] <- mohinora_DATA_rTp$coords[,1:2]

STABLE <- matrix(nrow=nrow(mohinora_DATA_rTp$coords), ncol=3)
STABLE[,1:2] <- mohinora_DATA_rTp$coords[,1:2]

YEARS <- matrix(nrow=nrow(mohinora_DATA_rTp$coords), ncol=3)
YEARS[,1:2] <- mohinora_DATA_rTp$coords[,1:2]

CP <- matrix(nrow=nrow(mohinora_DATA_rTp$coords), ncol=3)
CP[,1:2] <- mohinora_DATA_rTp$coords[,1:2]

# ---  

progressReportFile <- paste0( getwd(), "/RData/progressReports/mohinora/progress_cps.txt" )
file.create(path=progressReportFile, showWarnings=FALSE)

write("===BFAST01 analysis began at===",
      file=progressReportFile, append=TRUE)
write(as.character(Sys.time()[1]), file=progressReportFile,
      append=TRUE)

numCores <- detectCores()

kluster <- parallel::makeCluster(10, outfile="")
registerDoParallel(kluster)

output <- foreach(i=1:nrow(mohinora_DATA_rTp$coords), .combine="rbind",
                  .packages=c("bfast") ) %dopar% { 
                    
                    if(i %% 100 ==0){
                      texto <- paste0("Working on ROW: ", i)
                      write(texto, file=progressReportFile, append=TRUE)
                    }
                    
                    pixel <- mohinora_DATA_rTp$values[i,] * 1e-4
                    
                    pixel_ts <- ts(pixel, 
                                   start = c(2000,1), 
                                   end = c(2024, 23),
                                   frequency = 23)
                    
                    pixel_bfast01 <- bfast01( data=pixel_ts )
                    
                    # TYPE, SIGN, STABLE, YEARS, CP
                    
                    TEMP <- bfast01classify(pixel_bfast01)
                    
                    YEAR <- getYear(start=2000, end=2024, 
                                    bp=pixel_bfast01$breakpoints, freq=23)
                    
                    s <- c(TEMP$flag_type, TEMP$flag_significance, 
                           TEMP$flag_pct_stable, YEAR, 
                           pixel_bfast01$breakpoints)
                    
                    return(s)
                  }
stopCluster(kluster)

write( as.character(Sys.time()[1]), file=progressReportFile, append=TRUE )
write( "===BFAST01 analysis ended here===", file=progressReportFile, append=TRUE )

# ---

TYPE[,3] <- output[,1]
SIGN[,3] <- output[,2]
STABLE[,3] <- output[,3]
YEARS[,3] <- output[,4]
CP[,3] <- output[,5]

# --- asegurarse de haber creado /RData/mohinora_cps

RData_cps <- paste0( getwd(), "/RData/mohinora_cps" )

dir.create(RData_cps)

save(TYPE, file=paste0(RData_cps, "/TYPE.RData"))
save(SIGN, file=paste0(RData_cps, "/SIGN.RData"))
save(STABLE, file=paste0(RData_cps, "/STABLE.RData"))
save(YEARS, file=paste0(RData_cps, "/YEARS.RData"))
save(CP, file=paste0(RData_cps, "/CP.RData"))

# -----------------------
# --- RASTERIZACION --- #
# -----------------------

PROJECTION <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs"

map_TYPE <- matrixToRaster(matrix=TYPE, projection=PROJECTION)
map_SIGN <- matrixToRaster(matrix=SIGN, projection=PROJECTION)
map_YEARS <- matrixToRaster(matrix=YEARS, projection=PROJECTION)

# --- asegurarse de haber creado /data/outputs/mohinora_cps

outputs_cps <- paste0( getwd(), "/data/outputs/mohinora_cps" )

dir.create(outputs_cps)

raster::writeRaster(map_TYPE,
                    filename = paste0( outputs_cps, "/map_TYPE" ),
                    format="GTiff", datatype="INT2S", overwrite=TRUE)

raster::writeRaster(map_YEARS,
                    filename = paste0( outputs_cps, "/map_YEARS" ),
                    format="GTiff", datatype="INT2S", overwrite=TRUE)

raster::writeRaster(map_SIGN,
                    filename = paste0( outputs_cps, "/map_SIGN" ),
                    format="GTiff", datatype="INT2S", overwrite=TRUE)

# --- just the ANP Cerro Mohinora

TYPEmap <- rast(paste0( outputs_cps, "/map_TYPE.tif" ))

SIGNmap <- rast(paste0( outputs_cps, "/map_SIGN.tif" ))

signTypeMap <- TYPEmap
signTypeMap[ SIGNmap != 0 ] <- NA

mohinora_cps_TYPE <- terra::crop(x=signTypeMap, 
                                 y=mohinora_shp,
                                 mask=TRUE)

# raster::writeRaster(mohinora_cps_TYPE,
#                     filename = paste0( getwd(), "/TIF/mohinora_cps/map_signTYPE.tif"),
#                     datatype="INT2S", overwrite=TRUE)

# -----------------------
# --- VISUALIZACION --- #
# -----------------------

# --- objetos auxiliares para el tmap

COLORES_update <- c("#50C878", "#C08F73",
                    "#006B3C", "#E32636",
                    "#87A96B", "#2F221E",
                    "#F2C185", "#66FF00") # Colores para TIPOS de tendencia

usv_COLORS <- c("#A1E5A5", "#E9D66B", "#00A877", "#66B032", "#83A4F0", "#FC8FAB", "#F500A1")# Colores para tipos de SUELO y VEGETACION
usv_NAMES <- c("Pino", "Pastizal", "Pino-Encino",
               "Ayarin", "Agro", "Arbustiva", "Arborea")

COLOR_USV <- c(usv_COLORS[1], 
               rep(usv_COLORS[3],2), 
               usv_COLORS[4], 
               rep(usv_COLORS[2],2),
               rep(usv_COLORS[5], length(7:12)), 
               rep(usv_COLORS[6], 9),
               usv_COLORS[7])

# --- definiendo un bbox a usar en el tmap

bbox_new <- st_bbox(mohinora_shp_usv7) # current bounding box

xrange <- bbox_new$xmax - bbox_new$xmin # rango x
yrange <- bbox_new$ymax - bbox_new$ymin # rango y

bbox_new[1] <- bbox_new[1] - (0.25 * xrange) # xmin - left
bbox_new[3] <- bbox_new[3] + (0.25 * xrange) # xmax - right
bbox_new[2] <- bbox_new[2] - (0.3 * yrange) # ymin - bottom
bbox_new[4] <- bbox_new[4] + (0.25 * yrange) # ymax - top

bbox_new <- bbox_new %>%  
  st_as_sfc() 

# --- definiedo un st multistring para usarlo en tm_lines()

visual_mohinora <- mohinora_shp_usv7 %>%
  sf::st_cast("MULTILINESTRING")

visual_mohinora$COLOR <- COLOR_USV

# ---

type_map = tm_shape(mohinora_cps_TYPE, 
                    bbox = bbox_new) +
  tm_raster(col.scale = tm_scale(), 
            # style = "cont", 
            palette = COLORES_update, 
            legend.show = TRUE,
            title="Trend types") +
  tm_shape(visual_mohinora) + tm_lines(col="COLOR", lwd=3) + 
  tm_title("", frame = FALSE, bg.color = NA) +
  tm_compass(type = "8star", position = c("right", "bottom")) +
  tm_scalebar(text.size = 0.65,
               position = c("right", "bottom")) +
  tm_add_legend("symbol", 
                labels=usv_NAMES, 
                fill=usv_COLORS,
                border.col = "grey40",
                size=1,
                shape=18,
                is.portrait = TRUE)

type_map

# ---

# -------------------------------
# --- RESUMEN dE RESULTADOS --- #
# --- Pendiente, Tercera Clas ---
# -------------------------------

getCorrectPercent <- function(x){
  x_na <- x %>%
    filter(!is.na(map_TYPE))
  
  types <- sort(unique(x_na$map_TYPE))
  
  y <- unlist(lapply(types, function(s) sum(x_na$fraction[x_na$map_TYPE == s]) ))
  
  cbind.data.frame(TYPE=types, PERCENT=round(y/sum(x_na$fraction) * 100, digits=3) )
}


trunk <- mohinora_USV7 %>%
  dplyr::select( DESCRIPCIO ) %>%
  st_drop_geometry() 

df_exact <- trunk %>%
  apply(MARGIN = 1,
        FUN = function(x) terra::extract(mohinora_cps_TYPE,
                                         dplyr::filter(mohinora_USV7,
                                                       DESCRIPCIO == x),
                                         exact = TRUE)) 

names(df_exact) <- trunk$DESCRIPCIO

df_exact <- df_exact[!duplicated(df_exact)]

df_exact_percent <- lapply(df_exact, function(s) getCorrectPercent(s)) 

# ---

percent_type <- matrix(nrow=7, ncol=8)

for(i in 1:nrow(percent_type)){
  percent_type[i,as.numeric(df_exact_percent[[i]]$TYPE)] <- as.numeric(df_exact_percent[[i]]$PERCENT) 
}

df_type <- data.frame("Pino" = percent_type[1,],
                      "PinoEncino"= percent_type[2,],
                      "Ayarin" = percent_type[3,],
                      "Pastizal" = percent_type[4,],
                      "Agro" = percent_type[5,],
                      "ArbustivaPino" = percent_type[6,],
                      "ArboreaPino" = percent_type[7,])

row.names(df_type) <- 1:8

df_type[is.na(df_type)] <- 0

t(df_type) %>%
  kbl(digits=2, 
      caption = "Porcentaje de área clasificada por tipo de tendencias .") %>%
  kable_minimal(full_width = FALSE, html_font = "Cambria",
                font_size=20)

# --- Hermoseando la tabla anterior

usv_COLORS_sorted <- c("#A1E5A5", "#00A877", "#66B032",
                       "#E9D66B", "#83A4F0", "#FC8FAB", 
                       "#F500A1")

usv_NAMES_sorted <- c("Pino", "Pino-Encino", "Ayarin",
                      "Pastizales", "Agro", "Arbustiva", "Arborea")

t(df_type) %>%
  kbl(digits=2, booktabs = TRUE,
      caption = "Porcentaje de distintos tipos de
      tendencia por tipo de uso de suelo y vegetación usando 'terra::exact'") %>%
  kable_minimal(full_width = FALSE, font_size=20) %>%
  # kable_styling(latex_options = "striped", full_width = FALSE,
  #               font_size = 10) %>%
  column_spec(1, color = usv_COLORS_sorted) %>%
  row_spec(row=0, bold = TRUE)
# color = COLORES_update)

# ----

