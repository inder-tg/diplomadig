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
library(sephora)
library(sf)
library(terra)
library(foreach)
library(doParallel)

source( paste0( getwd(), "/Rscripts/auxFUN.R" ) )

# ---

# --- DATA Loading

dataDIR <- paste0( getwd(), "/data" )

listDIRS <- list.dirs(path = dataDIR)[-1]

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

# --- Exploración sobre Pastizal (coord 2168 tiene un pixel de textbook)

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
                            clusterSize=15, numFreq=2, k=2)

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

plot(output_pastizal, type="derivatives")

str(output_pastizal)

output_pastizal$phenoparams

# =============================================================================

# Fecha fenológica | Notation |  Definición
#       Green up   |      GU  |  máximo (global)  2da derivada
#   Start of season|      SoS |  máximo (global) 1era derivada 
#         Maturity |      Mat |  mínimo (global) 2da derivada

#       Senescence |      Sen |  mínimo (local) 2da derivada 
#   End of season  |      EoS |  mínimo (global) 1era. derivada
#        Dormancy  |      Dor |  máximo (local) 2da. derivada

# =============================================================================

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
                            clusterSize=15, numFreq=3, k=2)

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

output_agro$phenoparams

# =============================================================================

# Fecha fenológica | Notation |  Definición
#       Green up   |      GU  |  máximo (global)  2da derivada
#   Start of season|      SoS |  máximo (global) 1era derivada 
#         Maturity |      Mat |  mínimo (global) 2da derivada

#       Senescence |      Sen |  mínimo (local) 2da derivada 
#   End of season  |      EoS |  mínimo (global) 1era. derivada
#        Dormancy  |      Dor |  máximo (local) 2da. derivada

# =============================================================================

# --- PLOT base

image(subset(mohinora_DATA_interpol, 453), col=rev(terrain.colors(255)))
legend("topright", legend = usv_NAMES, col = usv_COLORS, lwd=4)
lines(mohinora_USV7, col=COLOR_USV, lwd=6)

# --- Exploración sobre Ayarin (coord 2734 provee buen ejemplo)

# --- SELECCIONA un píxel al interior del polígono "Agro"
# --- YA TE LA SABES:
# --- 1. Ejecuta la línea de abajo
XY <- locator()
# --- 2. Haz click (SOLO UNA VEZ) en algún píxel en la imagen
# --- 3. Presiona la tecla ESC de tu teclado
# --- 4. Continúa con el script a partir de la línea 76

xy <- get_timeSeries_byClicking(c(XY$x, XY$y), 
                                df=mohinora_DATA_interpol_rTp$coords) #xy$coord = 210

pixel_ayarin <- mohinora_DATA_interpol_rTp$values[xy$coord, ]

# ---

BASIS <- drbasis(n=50, q=2)

output_ayarin <- phenopar(x=pixel_ayarin * 1e-4,
                          startYear=2000, endYear=2023, 
                          frequency=23, basis=BASIS,
                          distance="dtw_basic", 
                          clusterSize=15, numFreq=3, k=2)
# ---

plot(pixel_ayarin,
     startYear=2000, endYear=2023, frequency=23)

plot(c(t(output_ayarin$x_smooth)), 
     startYear=2000, endYear=2023, frequency=23)

plot(c(t(output_ayarin$x_smooth)),
     startYear=2000, endYear=2023, frequency=23,
     type="profiles", position_legend="bottom")

plot(output_ayarin, type="ms")

sephora_fpca(x_sephora = output_ayarin, x=pixel_ayarin)

plot(output_ayarin, type="derivatives")

output_ayarin$phenoparams

# =============================================================================

# Fecha fenológica | Notation |  Definición
#       Green up   |      GU  |  máximo (global)  2da derivada
#   Start of season|      SoS |  máximo (global) 1era derivada 
#         Maturity |      Mat |  mínimo (global) 2da derivada

#       Senescence |      Sen |  mínimo (local) 2da derivada 
#   End of season  |      EoS |  mínimo (global) 1era. derivada
#        Dormancy  |      Dor |  máximo (local) 2da. derivada

# =============================================================================

# --- Uso de argumento series, usar coord 2734

pixel_ayarin <- mohinora_DATA_interpol_rTp$values[2734, ]

BASIS <- drbasis(n=50, q=2)

output_ayarin <- phenopar(x=pixel_ayarin * 1e-4,
                          startYear=2000, endYear=2023, 
                          frequency=23, basis=BASIS,
                          distance="dtw_basic", 
                          clusterSize=15, numFreq=3, k=2)

test_ayarin <- pixel_ayarin

# --- CLUSTER 1

BASIS <- drbasis(n=50, q=2)

CLUSTER1 <- which(output_ayarin$clustering@cluster == 1)

cluster1_ayarin <- phenopar(x=test_ayarin * 1e-4,
                            startYear=2000, endYear=2023, 
                            frequency=23, basis=BASIS,
                            distance="dtw_basic", 
                            series = CLUSTER1,
                            clusterSize=15, numFreq=2, k=2)

CLUSTER2 <- which(output_ayarin$clustering@cluster == 2)

cluster2_ayarin <- phenopar(x=test_ayarin * 1e-4,
                            startYear=2000, endYear=2023, 
                            frequency=23, basis=BASIS,
                            distance="dtw_basic", 
                            series = CLUSTER2,
                            clusterSize=15, numFreq=2, k=2)

# ---

plot(test_ayarin,
     startYear=2000, endYear=2023, frequency=23)

plot(c(t(cluster1_ayarin$x_smooth)), 
     startYear=2000, endYear=2023, frequency=23)

plot(c(t(cluster1_ayarin$x_smooth)),
     startYear=2000, endYear=2023, frequency=23,
     type="profiles", position_legend="bottom")

plot(cluster1_ayarin, type="ms")

# ---

sephora_fpca(x_sephora = cluster1_ayarin, x=test_ayarin)

plot(cluster1_ayarin, type="derivatives")

cluster1_ayarin$phenoparams

cluster1_ayarin$usedTotal

cluster1_ayarin$series

# ---

sephora_fpca(x_sephora = cluster2_ayarin, x=test_ayarin)

plot(cluster2_ayarin, type="derivatives")

cluster2_ayarin$phenoparams

cluster2_ayarin$usedTotal

cluster2_ayarin$series

# --- PLOT base

image(subset(mohinora_DATA_interpol, 453), col=rev(terrain.colors(255)))
legend("topright", legend = usv_NAMES, col = usv_COLORS, lwd=4)
lines(mohinora_USV7, col=COLOR_USV, lwd=6)

# --- Exploración sobre Pino-Encino

# --- SELECCIONA un píxel al interior del polígono "Agro"
# --- YA TE LA SABES:
# --- 1. Ejecuta la línea de abajo
XY <- locator()
# --- 2. Haz click (SOLO UNA VEZ) en algún píxel en la imagen
# --- 3. Presiona la tecla ESC de tu teclado
# --- 4. Continúa con el script a partir de la línea 76

xy <- get_timeSeries_byClicking(c(XY$x, XY$y), 
                                df=mohinora_DATA_interpol_rTp$coords) #xy$coord = 210

pixel_pinoEncino <- mohinora_DATA_interpol_rTp$values[xy$coord, ]

# ---

BASIS <- drbasis(n=50, q=2)

output_pinoEncino <- phenopar(x=pixel_pinoEncino * 1e-4,
                          startYear=2000, endYear=2023, 
                          frequency=23, basis=BASIS,
                          distance="dtw_basic", 
                          clusterSize=15, numFreq=2, k=2)

# ---

plot(pixel_pinoEncino,
     startYear=2000, endYear=2023, frequency=23)

plot(c(t(output_pinoEncino$x_smooth)), 
     startYear=2000, endYear=2023, frequency=23)

plot(c(t(output_pinoEncino$x_smooth)),
     startYear=2000, endYear=2023, frequency=23,
     type="profiles", position_legend="bottom")

plot(output_pinoEncino, type="ms")

sephora_fpca(x_sephora = output_pinoEncino, x=pixel_pinoEncino)

plot(output_pinoEncino, type="derivatives")

output_pinoEncino$phenoparams

# =============================================================================

# Fecha fenológica | Notation |  Definición
#       Green up   |      GU  |  máximo (global)  2da derivada
#   Start of season|      SoS |  máximo (global) 1era derivada 
#         Maturity |      Mat |  mínimo (global) 2da derivada

#       Senescence |      Sen |  mínimo (local) 2da derivada 
#   End of season  |      EoS |  mínimo (global) 1era. derivada
#        Dormancy  |      Dor |  máximo (local) 2da. derivada

# =============================================================================


# ----------------------------------------
# --- sephora: cómputo a gran escala --- #
# ----------------------------------------

mohinora_pastizal <- crop(mohinora_DATA_interpol,
                          dplyr::filter(mohinora_USV7,
                                        DESCRIPCIO == "PASTIZAL INDUCIDO"),
                          mask=TRUE)

mohinora_pastizal_rTp <- spRast_valueCoords(mohinora_pastizal, na_rm = TRUE)

# --- CÓMPUTO EN PARALELO

source( paste0( getwd(), "/Rscripts/basic_auxFUN.R" ) )
source( paste0( getwd(), "/Rscripts/adv_auxFUN.R" ) )

BASIS <- drbasis(n=50, q=2)

# --- Antes de ejecutar la siguiente línea
# --- asegurarse de crear "/RData/mohinora_sephora/output_numFreq_2_k_" 

phenopar_polygon_test(product = "independent",
                      data = mohinora_pastizal_rTp$values,
                      numFreq = 2, distance = "dtw_basic",
                      clusterSize = 15, basis = BASIS,
                      k=2, numCores = 7, 
                      save_some = c("clustering", "fpca", "phenoparams"),
                      outputFileBaseName = "pastizal",
                      dirToSave = paste0( getwd(), "/RData/mohinora_sephora/output_numFreq_2_k_2" ))


# --- EXTRACTION

paramsList_round1 <- extractPhenoParams(data=mohinora_params, 
                                        L=nrow(mohinora_pastizal_rTp$values))

GU_round1 <- getParam(x=paramsList_round1$phenoParams, phenoparam = "GU")
SoS_round1 <- getParam(x=paramsList_round1$phenoParams, phenoparam = "SoS")
Mat_round1 <- getParam(x=paramsList_round1$phenoParams, phenoparam = "Mat")
Sen_round1 <- getParam(x=paramsList_round1$phenoParams, phenoparam = "Sen")
EoS_round1 <- getParam(x=paramsList_round1$phenoParams, phenoparam = "EoS")
Dor_round1 <- getParam(x=paramsList_round1$phenoParams, phenoparam = "Dor")

# GU_round1 <- GU_round1[!sapply(GU_round1, is.null)]
# GU_round1 <- unlist(GU_round1)
# length(GU_round1)
# 
# SoS_round1 <- SoS_round1[!sapply(SoS_round1, is.null)]
# length(unlist(SoS_round1))
# SoS_round1 <- unlist(SoS_round1)
# 
# Mat_round1 <- Mat_round1[!sapply(Mat_round1, is.null)]
# length(unlist(Mat_round1))
# Mat_round1 <- unlist(Mat_round1)
# 
# EoS_round1 <- EoS_round1[!sapply(EoS_round1, is.null)]
# length(unlist(EoS_round1))
# EoS_round1 <- unlist(EoS_round1)


phenoParams <- list(nrow( mohinora_pastizal_rTp$coords ))
phenoParams$coords <- matrix(nrow=nrow(mohinora_pastizal_rTp$coords),ncol=2)
phenoParams$coords <- mohinora_pastizal_rTp$coords

phenoParams$GU <- numeric(nrow(mohinora_pastizal_rTp$coords))
phenoParams$GU <- GU_round1

phenoParams$SoS <- numeric(nrow(mohinora_pastizal_rTp$coords))
phenoParams$SoS <- SoS_round1

phenoParams$Mat <- numeric(nrow(mohinora_pastizal_rTp$coords))
phenoParams$Mat <- Mat_round1

phenoParams$EoS <- numeric(nrow(mohinora_pastizal_rTp$coords))
phenoParams$EoS <- EoS_round1

phenoParams$Sen <- numeric(nrow(mohinora_pastizal_rTp$coords))
phenoParams$Sen <- Sen_round1

phenoParams$Dor <- numeric(nrow(mohinora_pastizal_rTp$coords))
phenoParams$Dor <- Dor_round1

# =====================
# --- RASTERIZATION ---
# =====================

PROYECTION <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs"

GU_mat <- cbind(phenoParams$coords, phenoParams$GU)
SoS_mat <- cbind(phenoParams$coords, phenoParams$SoS)
Mat_mat <- cbind(phenoParams$coords, phenoParams$Mat)
EoS_mat <- cbind(phenoParams$coords, phenoParams$EoS)
# Sen_mat <- cbind(phenoParams$coords, phenoParams$Sen)
# Dor_mat <- cbind(phenoParams$coords, phenoParams$Dor)

GU <- matrixToRaster(matrix = GU_mat,
                     projection = PROYECTION)

SoS <- matrixToRaster(matrix = SoS_mat,
                     projection = PROYECTION)

Mat <- matrixToRaster(matrix = Mat_mat,
                      projection = PROYECTION)

EoS <- matrixToRaster(matrix = EoS_mat,
                      projection = PROYECTION)

# Sen <- matrixToRaster(matrix = Sen_mat,
#                       projection = PROYECTION)
# 
# Dor <- matrixToRaster(matrix = Dor_mat,
#                       projection = PROYECTION)

# =====================
# --- SPIRALIZATION ---
# =====================

# --- colors used in spiralPlot below
cgu <- rgb(173/255,221/255,142/255)
csos <- rgb(120/255,198/255,121/255)
cmat <- rgb(49/255, 163/255,84/255)
csen <- rgb(217/255, 95/255, 14/255)
ceos <- rgb(254/255, 153/255, 41/255)
cdor <- rgb(208/255, 209/255, 230/255)

colores <- c(cgu,csos,cmat,csen,ceos,cdor)

# ---

Sen_intervenido <- sapply(Sen_round1, function(s) ifelse(length(s) == 0, NA, s) )

Dor_intervenido <- sapply(Dor_round1, function(s) ifelse(length(s) == 0, NA, s ) )

getSpiralPlot(MAT=cbind(GU_round1, SoS_round1, Mat_round1,
                        Sen_intervenido, EoS_round1, Dor_intervenido), 
              LABELS=month.name,
              vp_param=list(width=0.5, height=0.7))
vcd::grid_legend(x=1.215, y=0.125, pch=18, col=colores,
                 frame=FALSE,
                 labels=c("GU","SoS","Mat","Sen","EoS","Dor"),
                 title="Params")
















































