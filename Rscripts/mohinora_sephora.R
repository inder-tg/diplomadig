# --- Diplomado Geomática, IG, UNAM, 2024
# --- Módulo IX Percepción Remota || Análisis de series de tiempo de imágenes satelitales con R

# --- Elaborado: Agosto 21, 2024
# --- Impartido: Agosto 23, 2024

# --- Más detalles técnicos sobre este análisis se pueden encontrar en "Phenology curve estimation via a 
# --- mixed model representation of functional principal components: Characterizing time series of satellite-derived 
# --- vegetation indices" de Tecuapetla et al. 2024, disponible en https://arxiv.org/abs/2403.14451

# --- DATASET: NDVI MOD13Q1 v061 en Cerro Mohinora, Chihuahua, 2000-2023 

# --- ADDicionalmente, este script requiere archivos
# --- MOD13Q1_061_250m_16_days_NDVI_interpol.tif

# --- NOTA: Este código no es completamente automático; de vez en vez se requerirá
# --- crear folders o descomentar líneas de código (por ejemplo en el uso de rutinas ligada al paquete raster),
# --- esas líneas están marcada con el texto "ACTION REQUIRED!!!"

# --- Preámbulo
library(terra)
# library(raster)
# library(mapview)
# library(RColorBrewer)
# library(gtools)
# library(foreach)
# library(doParallel)
# library(bfast)
# library(tmap)
# library(geoTS)
library(sf)
# library(tidyverse)
# library(kableExtra)
library(sephora)

source( paste0( getwd(), "/Rscripts/auxFUN.R" ) )

# ---

# --- DATA Loading

dataDIR <- paste0( getwd(), "/data" )

listDIRS <- list.dirs(path = dataDIR)[-1]

listFILES_mohinora <- list.files(path=paste0(getwd(), "/TIF"), 
                                 pattern=".tif$", 
                                 full.names=TRUE)

listFILES_mohinora <- list.files(path=listDIRS[3], # cambiar 3 por 2 
                                 pattern=".tif$", 
                                 full.names=TRUE)

shpFILES_mohinora <- list.files(path=listDIRS[6], # cambiar 6 por 4 
                                pattern=".shp$", 
                                full.names=TRUE)

# mohinora_DATA_interpol <- stack( listFILES_mohinora[1] ) # ACTION REQUIRED!!!
mohinora_DATA_interpol <- rast( listFILES_mohinora[1] )

# mohinora_shp <- shapefile( shpFILES_mohinora[1] ) # ACTION REQUIRED!!!
mohinora_shp <- read_sf( shpFILES_mohinora[1] )

mohinora_USV7 <- read_sf( shpFILES_mohinora[2] )

# mohinora_DATA_interpol_rTp <- rasterToPoints(stack_NDVI_Mohinora) # ACTION REQUIRED!!!
mohinora_DATA_interpol_rTp <- spRast_valueCoords(mohinora_DATA_interpol)

# ------------------------------
# --- sephora: Exploración --- #
# ------------------------------

usv_COLORS <- c("#A1E5A5", "#E9D66B", "#00A877", 
                "#66B032", "#83A4F0", "#FC8FAB", "#F500A1")

usv_NAMES <- c("Pino", "Pastizal", "Pino-Encino",
               "Ayarin", "Agro", "Arbustiva", "Arborea")

COLOR_USV <- c(usv_COLORS[1], 
               rep(usv_COLORS[3],2), 
               usv_COLORS[4], 
               rep(usv_COLORS[2],2),
               rep(usv_COLORS[5], length(7:12)), 
               rep(usv_COLORS[6], 9),
               usv_COLORS[7])

# --- PLOT base

image(subset(mohinora_DATA_interpol, 453), col=rev(terrain.colors(255)))
legend("topright", legend = usv_NAMES, col = usv_COLORS, lwd=4)
lines(mohinora_USV7, col=COLOR_USV, lwd=6)

# --- Exploración sobre Pastizal

# --- SELECCIONA un píxel al interior del polígono "Pastizal"
# --- YA TE LA SABES:
# --- 1. Ejecuta la línea de abajo
XY <- locator()
# --- 2. Haz click (SOLO UNA VEZ) en algún píxel en la imagen
# --- 3. Presiona la tecla ESC de tu teclado
# --- 4. Continúa con el script a partir de la línea 76

xy <- get_timeSeries_byClicking(c(XY$x, XY$y), 
                                df=mohinora_DATA_interpol_rTp$coords)

pixel_pastizal <- mohinora_DATA_interpol_rTp$values[xy$coord, ]

# ---

BASIS <- drbasis(n=50, q=2)

output_pastizal <- phenopar(x=pixel_pastizal * 1e-4,
                            startYear=2000, endYear=2023, 
                            frequency=23, basis=BASIS,
                            distance="dtw_basic", 
                            clusterSize=15, numFreq=3, k=2)

# ---

plot(pixel_pastizal,
     startYear=2000, endYear=2023, frequency=23)

plot(c(t(output_pastizal$x_smooth)), 
     startYear=2000, endYear=2023, frequency=23)

plot(c(t(output_pastizal$x_smooth)),
     startYear=2000, endYear=2023, frequency=23,
     type="profiles", position_legend="bottom")

plot(output_pastizal, type="ms")

sephora_fpca(x_sephora = output_pastizal, x=pixel_pastizal)

# =============================================================================

# Fecha fenológica | Notation |  Definición
#       Green up   |      GU  |  máximo (global)  2da derivada
#   Start of season|      SoS |  máximo (global) 1era derivada 
#         Maturity |      Mat |  mínimo (global) 2da derivada

#       Senescence |      Sen |  mínimo (local) 2da derivada 
#   End of season  |      EoS |  mínimo (global) 1era. derivada
#        Dormancy  |      Dor |  máximo (local) 2da. derivada

# =============================================================================

plot(output_pastizal, type="derivatives")

str(output_pastizal)

output_pastizal$phenoparams

# --- PLOT base

image(subset(mohinora_DATA_interpol, 453), col=rev(terrain.colors(255)))
legend("topright", legend = usv_NAMES, col = usv_COLORS, lwd=4)
lines(mohinora_USV7, col=COLOR_USV, lwd=6)

# --- Exploración sobre Agricultura

# --- SELECCIONA un píxel al interior del polígono "Agro"
# --- YA TE LA SABES:
# --- 1. Ejecuta la línea de abajo
XY <- locator()
# --- 2. Haz click (SOLO UNA VEZ) en algún píxel en la imagen
# --- 3. Presiona la tecla ESC de tu teclado
# --- 4. Continúa con el script a partir de la línea 76

xy <- get_timeSeries_byClicking(c(XY$x, XY$y), 
                                df=mohinora_DATA_interpol_rTp$coords) #xy$coord = 210

pixel_agro <- mohinora_DATA_interpol_rTp$values[xy$coord, ]

# ---

BASIS <- drbasis(n=50, q=2)

output_agro <- phenopar(x=pixel_agro * 1e-4,
                            startYear=2000, endYear=2023, 
                            frequency=23, basis=BASIS,
                            distance="dtw_basic", 
                            clusterSize=15, numFreq=3, k=3)

# ---

plot(pixel_agro,
     startYear=2000, endYear=2023, frequency=23)

plot(c(t(output_agro$x_smooth)), 
     startYear=2000, endYear=2023, frequency=23)

plot(c(t(output_agro$x_smooth)),
     startYear=2000, endYear=2023, frequency=23,
     type="profiles", position_legend="bottom")

plot(output_agro, type="ms")

sephora_fpca(x_sephora = output_agro, x=pixel_agro)

plot(output_agro, type="derivatives")

str(output_agro)

output_agro$phenoparams

# Fecha fenológica | Notation |  Definición
#       Green up   |      GU  |  máximo (global)  2da derivada
#   Start of season|      SoS |  máximo (global) 1era derivada 
#         Maturity |      Mat |  mínimo (global) 2da derivada

#       Senescence |      Sen |  mínimo (local) 2da derivada 
#   End of season  |      EoS |  mínimo (global) 1era. derivada
#        Dormancy  |      Dor |  máximo (local) 2da. derivada




































