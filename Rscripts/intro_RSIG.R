
# --- Diplomado Geomática, IG, UNAM, 2023
# --- Módulo X Análisis de series de tiempo de imágenes satelitales con R
# --- Diplomado Geomática, IG, UNAM, 2024
# --- Módulo IX Percepción Remota: Análisis de series de tiempo de imágenes satelitales con R


# --- Elaborado: Junio 2, 2022
# --- Impartido: Junio 24, 2022
# --- Modificado: Julio 25, 2023; Agosto 9, 2024

# --- En este script presentamos ejemplos para usar los paquetes raster,
# --- mapview, RColorBrewer y discutimos algunas ideas para realizar análisis 
# --- exploratorio de imágenes satelitales.
# --- Basado en Cap. 2 de Remote sensing image analysis disponible en
# --- http://rspatial.org/analysis/rst/9-remotesensing.html

# --- DATASET: Escena Landsat 8 tomada 14 Junio, 2017; cubre el área entre Concord y Stockton

# --- NOTA: Abajo las líneas de código comentadas con el formato
# --- "# CODE # === son compatibles con el paquete raster; para tener un ejemplo
# --- de una línea de código comentada con este formato ver la línea 29

# --- Preámbulo

# --- Instalación de todos los paquetes a utilizar en este módulo
source( paste0( getwd(), "/Rscripts/auxPKG.R" ) )

# library(raster) # ===
library(terra)
library(mapview)
library(RColorBrewer)
library(gtools)

# --- Lectura de archivos
dirDATA_rspatial <- paste0( getwd(), "/data/rspatial" )

listFILES_rspatial <- mixedsort(list.files( path=dirDATA_rspatial,
                                            pattern=".tif", full.names=TRUE ))

# ------------------------------------------------------------------------------
# Escena Landsat 8 tomada 14 Junio, 2017; cubre el área entre Concord y Stockton
# LANDSAT <- stack(listFILES_rspatial) # ===
LANDSAT <- rast(listFILES_rspatial)
# ------------------------------------------------------------------------------

# --- Visualización

par(mfrow=c(2,2))

plot(subset(LANDSAT, 2), main = "Landsat Blue")
plot(subset(LANDSAT, 3), main = "Landsat Green")
plot(subset(LANDSAT, 4), main = "Landsat Red")
plot(subset(LANDSAT, 5), main = "Landsat NIR")

# LANDSAT_RGB <- stack(listFILES_rspatial[c(4,3,2)]) # ===
# LANDSAT_FCC <- stack(listFILES_rspatial[c(5,4,3)]) # ===
LANDSAT_RGB <- rast(listFILES_rspatial[c(4,3,2)])
LANDSAT_FCC <- rast(listFILES_rspatial[c(5,4,3)])


# Compuesto de color verdadero 
par(mfrow = c(1,1))
plotRGB(LANDSAT_RGB, r = 1, g = 2, b = 3, axes = TRUE, stretch = "lin",  
        main = "Landsat True Color Composite")

# Compuesto de color verdadero y falso
par(mfrow = c(1, 2))
plotRGB(LANDSAT_RGB, r = 1, g = 2, b = 3, axes = TRUE, stretch = "lin",  
        main = "Landsat True Color Composite")
plotRGB(LANDSAT_FCC, r = 1, g = 2, b = 3, axes = TRUE, stretch = "lin",  
        main = "Landsat False Color Composite")

# Mapa interactivo

mp <- mapview(subset(LANDSAT,5))

mp

# --- Algunas herramientas para el análisis exploratorio

# ------------------------------------------------------------------------------
# Generando perfiles espectrales a partir de la combinación de varias funciones:
# subset, names, extent, crop, sample, readRDS, spsample, over, extract, aggregate,
# 
# ------------------------------------------------------------------------------

# ---------------
# subset y name #
# ---------------

# Remover las 4 últimas bandas de "LANDSAT"
landsat <- subset(LANDSAT, 1:7)

# Renombrando las bands:
names(landsat)

names(landsat) <- c("ultra-blue", "blue", "green", "red", "NIR", "SWIR1", "SWIR2")

names(landsat)

# --------------------------
# Spatial subset or "crop" #
# --------------------------

# extent(landsat) # ===
ext(landsat)

# newExtent <- extent(624387, 635752, 4200047, 4210939) # ===
newExtent <- ext(624387, 635752, 4200047, 4210939)

landsatCrop <- crop(landsat, newExtent)

par(mfrow=c(1,1))
par(mar = c(4.5, 5, 1.5, 2))
plotRGB(landsat, 1,2,3, stretch = "lin")
plotRGB(landsatCrop, 3,2,1, stretch = "lin", add = T)

# saving results to disk
# test <- getwd()
# setwd(paste(rootDir, "/data/rspatial", sep = ""))
# writeRaster(landsatCrop, filename = "cropped-landsat.tif", format = "GTiff", 
#             overwrite = TRUE)
# setwd(test)

# --------------------------------------------
# Extraer valores de un subconjunto del stack:
# readRDS, sample, readRDS, spsample, over, extract
# --------------------------------------------

# cargar polígonos con info de uso y cobertura del suelo
poligono <- readRDS(paste0(dirDATA_rspatial, "/samples.rds"))

# tomemos una muestra aleatoria de 300 
puntos_muestra <- spsample(poligono, 300, type = "random")

# agregar info a los puntos muestreados
puntos_muestra$class <- over(puntos_muestra, poligono)$class

# extract values with points
# mat_landsat_puntos <- extract(landsat, puntos_muestra) # ===
mat_landsat_puntos <- extract(landsat, puntos_muestra@coords)

# imprime en consola primeros 6 renglones
head(mat_landsat_puntos)

# ----------------------
# Perfiles espectrales #
# ----------------------

# -----------------------------------------------------------------------------
# Perfil espectral: gráfica de todas las bandas en pixeles que representa
# alguna característica de la superficie terrestre
# -----------------------------------------------------------------------------

ms <- aggregate(mat_landsat_puntos, list(puntos_muestra$class), mean)

rownames(ms) <- ms[,1]
ms <- ms[,-1]
ms

# plot de perfiles espectrales

# paleta de colores para las clases de cobertura del suelo
mycolor <- c("darkred", "yellow", "burlywood", "cyan", "blue")

# convertir ms de data.frame a matrix (una conveniencia)
ms <- as.matrix(ms)

# un plot "vacío"
plot(0, ylim = c(0, 0.6), xlim = c(1,7), type = "n", xlab = "Bands", 
     ylab = "Reflectance")

# aggregamos bandas vs reflectancia
for (i in 1:nrow(ms)) {
  lines(ms[i,], lwd = 3, col = mycolor[i])
}

# agregamos un título
title(main = "Spectral Profile from Landsat (subset)", font.main = 2)

# agregamos una leyenda
legend("topleft", rownames(ms), cex = 0.8, col = mycolor, lwd = 3, bty = "n" )

# --------------------------------
# Cálculo de índices espectrales #
# --------------------------------

ndvi <- (subset(landsat,5) - subset(landsat,4))/(subset(landsat,5) + subset(landsat,4))

plot(ndvi)

# -----------------------------------------------------------------------------
# --- Definiendo una función
# Esta función calcula un índice espectral utilizando como parámetros "imagen" 
# y el número de dos bandas, a saber "k" e "i"
spectralIndex <- function(image, k, i){
  bandK <- image[[k]] # get k-th band from image
  bandi <- image[[i]] # get i-th band from image
  index <- (bandK - bandi)/(bandK + bandi)
  index
}
# -----------------------------------------------------------------------------

# Para "landsat" NIR = 5, red = 4
ndvi_fun <- spectralIndex(landsat, 5, 4)
plot(ndvi_fun, col = rev(terrain.colors(10)), main = "NDVI Landsat8")


# display.brewer.all()
myBlues <- brewer.pal(9, "Blues")

waterIndex <- spectralIndex(landsat, 1, 7)
plot(waterIndex, col = myBlues, main = "waterIndex Landsat8")

# ------------------------
# histograma de ndvi_fun #
# ------------------------

hist(ndvi)

hist(ndvi,
     main = "Distribution of NDVI values",
     xlab = "NDVI",
     ylab="Frequency",
     col = "wheat",
     xlim = c(-0.5, 1),
     breaks = 30,
     xaxt = 'n')
axis(side=1, at = seq(-0.5,1, 0.05), labels = seq(-0.5,1, 0.05))

# --------------
# Thresholding #
# --------------

# ---
# Suponiendo que valores de NDVI mayores que 0.4 representan "verdor saludable" 
# Dibujamos zonas con verdor saludable

ndviVeg <- ndvi_fun
ndviVeg[ndviVeg < 0.4] <- NA

plot(ndviVeg, main = "Verdor saludable")

# --- ndviVeg es equivalente a vegNDVI ---
# vegNDVI <- calc(ndvi_fun, function(x){x[x < 0.4]<-NA; x}) # ===
vegNDVI <- app(ndvi_fun, function(x){x[x < 0.4]<-NA; x})
plot(vegNDVI, main = "Verdor saludable con app")

# ?reclassify
# vegReclassify <- reclassify( ndvi_fun, cbind(-Inf, 0.4, NA) ) # ===
vegReclassify <- classify( ndvi_fun, cbind(-Inf, 0.4, NA) )

# compareRaster(vegNDVI, vegReclassify) # ===
compareGeom(vegNDVI, vegReclassify)

# ---
# El histograma de ndvi_fun tiene un "pico" alrededor del intervalo (0.25, 0.3)
# abajo mostramos 2 formas de resaltar las áreas con valores en  


# land <- reclassify(ndvi_fun, c(-Inf, 0.25, NA, 0.25, 0.3, 1, 0.3, Inf, NA)) # ===
(matClass <- matrix(c(-Inf, 0.25, NA, 0.25, 0.3, 1, 0.3, Inf, NA), ncol=3, byrow = TRUE))

land <- classify(ndvi_fun,  matClass)
plot(land, main = "Qué es esto?")

# Sobreponiedo LANDSAT_RGB y land
plotRGB(LANDSAT_RGB, r = 1, g = 2, b = 3, axes = TRUE, stretch = "lin", 
        main = "Landsat False Color Composite")
plot(land, add = TRUE, legend = FALSE)

# Diferentes clases para NDVI
# vegc <- reclassify(ndvi_fun, c(-Inf, 0.25, 1, 0.25, 0.3, 2, 0.3, 0.4, 3, 0.4, 0.5, 
#                            4, 0.5, Inf, 5)) # ===

matClass <- matrix(c(-Inf, 0.25, 1, 0.25, 0.3, 2, 0.3, 0.4, 3, 0.4, 0.5, 
                     4, 0.5, Inf, 5), byrow = TRUE, ncol=3)
vegc <- classify(ndvi_fun, matClass)
plot(vegc, col = rev(terrain.colors(4)), main = 'NDVI based thresholding')
# -----------------------------------------------------------------------------
