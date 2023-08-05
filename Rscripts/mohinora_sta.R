# --- Diplomado Geomática, IG, UNAM, 2022
# --- Módulo VI Análisis de series de tiempo de imágenes satelitales con R

# --- Elaborado: Junio 2, 2022
# --- Impartido: Julio 1, 2022
# --- Actualizado: Agosto 4, 2023

# --- En este script presentamos un análisis de tendencias estacionales basado en 
# --- rutinas de los paquetes geoTS y sta.

# --- DATASET: NDVI MOD13Q1 en Cerro Mohinora, Chihuahua

# --- Preámbulo
library(raster)
library(mapview)
library(RColorBrewer)
library(foreach)
library(doParallel)
library(Kendall)
library(trend)

library(geoTS)
library(sta)

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

SHP_anp <- list.files( path = paste0( getwd(), "/data/anp_2021" ),
                       pattern = ".shp", 
                       full.names = TRUE)

shp_anp <- shapefile( SHP_anp[1] )

shp_anp_sinu <- spTransform(shp_anp, crs(stack_primeras3Imagenes))

# ---

# ----------------------------------------------------------
# --- Análisis de tendencias estacionales: Exploración --- #
# ----------------------------------------------------------

XY <- list(x=-10698785, y=2893289)

pixel_full <- as.numeric(get_timeSeries_byClicking(c(XY$x, XY$y),
                                        df=mohinora_NDVI_rTp_full)$ts) * 1e-4

pixel_harm <- haRmonics(y=pixel_full[1:23], numFreq = 3, delta = 0)

pixel_harm

yRan <- range(pixel_full[1:23], pixel_harm$fitted)

plot(pixel_full[1:23], type="l", xlab="", ylab="",
     ylim=yRan)
par(new=TRUE)
plot(pixel_harm$fitted, type="l", col="red", ylim=yRan,
     xlab="DoY", ylab="NDVI")

# ---

i <- 4

pixel_harm <- haRmonics(y=pixel_full[ (23*(i-1) + 1):(23 * i) ], numFreq = 3, 
                        delta = 0)

yRan <- range(pixel_full[(23*(i-1) + 1):(23 * i)], pixel_harm$fitted)

plot(pixel_full[(23*(i-1) + 1):(23 * i)], type="l", 
     xlab="", ylab="", ylim=yRan)
par(new=TRUE)
plot(pixel_harm$fitted, type="l", col="red", ylim=yRan,
     xlab="DoY", ylab="NDVI")

# ---

pixel_harm_full <- sta(data = pixel_full, freq=23, endYear = 2009)

yRan <- range(pixel_full, pixel_harm_full$fit)

plot(pixel_full, type="l", 
     xlab="", ylab="", ylim=yRan)
par(new=TRUE)
plot(pixel_harm_full$fit, type="l", col="red", ylim=yRan,
     xlab="Tiempo (años)", ylab="NDVI")

# --- como 'ts' (serie de tiempo)

pixel_ts <- ts(pixel_full, start = c(2000,1), end = c(2009,23), frequency = 23)

fit_ts <- ts(pixel_harm_full$fit, start = c(2000,1), end = c(2009,23), 
             frequency = 23)

plot(pixel_ts, type="l", col="lightgray", 
     xlab="", ylab="", ylim=yRan, lwd=3)
par(new=TRUE)
plot(fit_ts, type="l", col="red", ylim=yRan,
     xlab="Tiempo (años)", ylab="NDVI")

# ---------------------------------
# --- Seasonal Trend Analysis --- #
# ---------------------------------

str(pixel_harm_full)

# ---

MEAN <- pixel_harm_full$sta$mean$harmCoeffs

MEAN_MK <- mk.test(MEAN)

MEAN_MK$p.value

MEAN_ST <- sens.slope(MEAN)

MEAN_ST$estimates

b_hat <- as.numeric(MEAN_ST$estimates)
a_hat <- median( MEAN - b_hat * 1:length(MEAN) )

lineaTS <- a_hat + b_hat * 1:length(MEAN)

yRan <- range(MEAN, lineaTS)

plot(ts(MEAN, start=c(2000,1), end=c(2009,1), frequency = 1), 
     type="p", ylab = "NDVI", ylim=yRan)
par(new=TRUE)
plot(ts(lineaTS, start=c(2000,1), end=c(2009,1), frequency = 1), 
     type="l", col = "lightcoral", ylab = "NDVI", lwd=5, ylim=yRan)

# --- 

ANNUAL <- pixel_harm_full$sta$annual$harmCoeffs

ANNUAL_MK <- mk.test(ANNUAL)

ANNUAL_MK$p.value

ANNUAL_ST <- sens.slope(ANNUAL)

ANNUAL_ST$estimates

b_hat <- as.numeric(ANNUAL_ST$estimates)
a_hat <- median( ANNUAL - b_hat * 1:length(ANNUAL) )

lineaTS <- a_hat + b_hat * 1:length(ANNUAL)

yRan <- range(ANNUAL, lineaTS)

plot(ts(ANNUAL, start=c(2000,1), end=c(2009,1), frequency = 1), 
     type="p", ylab = "NDVI", pch=16, ylim=yRan)
par(new=TRUE)
plot(ts(lineaTS, start=c(2000,1), end=c(2009,1), frequency = 1), 
     type="l", col = "lightcoral", ylab = "NDVI", lwd=5, ylim=yRan)

# ---

sANNUAL <- pixel_harm_full$sta$semiannual$harmCoeffs

sANNUAL_MK <- mk.test(sANNUAL)

sANNUAL_MK$p.value

sANNUAL_ST <- sens.slope(sANNUAL)

sANNUAL_ST$estimates

b_hat <- as.numeric(sANNUAL_ST$estimates)
a_hat <- median( sANNUAL - b_hat * 1:length(sANNUAL) )

lineaTS <- a_hat + b_hat * 1:length(sANNUAL)

yRan <- range(sANNUAL, lineaTS)

plot(ts(sANNUAL, start=c(2000,1), end=c(2009,1), frequency = 1), 
     type="p", ylab = "NDVI", pch=16, ylim=yRan)
par(new=TRUE)
plot(ts(lineaTS, start=c(2000,1), end=c(2009,1), frequency = 1), 
     type="l", col = "lightcoral", ylab = "NDVI", lwd=5, ylim=yRan)

# ---

plot(pixel_harm_full)

# ------------------------------------
# --- STA: cómputo a gran escala --- #
# ------------------------------------

output <- sta(data=mohinora_NDVI_rTp_full * 1e-4,
    freq = 23,
    numFreq = 3,
    delta=0,
    endYear = 2009, save = TRUE,
    dirToSaveSTA = paste0(getwd(), "/RData/mohinora_sta"),
    numCores = 6)

str(output)

df_slope_mean <- matrix(nrow=nrow(mohinora_NDVI_rTp_full), ncol=3)
df_slope_mean[,1:2] <- mohinora_NDVI_rTp_full[,1:2]
df_slope_mean[,3] <- output$sta$mean[,3]

df_pVal_mean <- matrix(nrow=nrow(mohinora_NDVI_rTp_full), ncol=3)
df_pVal_mean[,1:2] <- mohinora_NDVI_rTp_full[,1:2]
df_pVal_mean[,3] <- output$sta$mean[,4]

# ---

df_slope_annual <- matrix(nrow=nrow(mohinora_NDVI_rTp_full), ncol=3)
df_slope_annual[,1:2] <- mohinora_NDVI_rTp_full[,1:2]
df_slope_annual[,3] <- output$sta$annual[,3]

df_pVal_annual <- matrix(nrow=nrow(mohinora_NDVI_rTp_full), ncol=3)
df_pVal_annual[,1:2] <- mohinora_NDVI_rTp_full[,1:2]
df_pVal_annual[,3] <- output$sta$annual[,4]

# ---

df_slope_sannual <- matrix(nrow=nrow(mohinora_NDVI_rTp_full), ncol=3)
df_slope_sannual[,1:2] <- mohinora_NDVI_rTp_full[,1:2]
df_slope_sannual[,3] <- output$sta$semiannual[,3]

df_pVal_sannual <- matrix(nrow=nrow(mohinora_NDVI_rTp_full), ncol=3)
df_pVal_sannual[,1:2] <- mohinora_NDVI_rTp_full[,1:2]
df_pVal_sannual[,3] <- output$sta$semiannual[,4]


# -----------------------
# --- RASTERIZACION --- #
# -----------------------

PROJECTION <- raster::projection(stack_primeras3Imagenes)

map_pVal_mean <- matrixToRaster_test(matrix=df_pVal_mean, projection=PROJECTION)
map_slope_mean <- matrixToRaster_test(matrix=df_slope_mean, projection=PROJECTION)

map_pVal_annual <- matrixToRaster_test(matrix=df_pVal_annual, projection=PROJECTION)
map_slope_annual <- matrixToRaster_test(matrix=df_slope_annual, projection=PROJECTION)

map_pVal_sannual <- matrixToRaster_test(matrix=df_pVal_sannual, projection=PROJECTION)
map_slope_sannual <- matrixToRaster_test(matrix=df_slope_sannual, projection=PROJECTION)

# -----------------------
# --- VISUALIZACION --- #
# -----------------------

pVal_mean_SHP <- map_pVal_mean
pVal_mean_SHP <- crop(pVal_mean_SHP, shp_anp_sinu[165,])
pVal_mean_SHP <- mask(pVal_mean_SHP, shp_anp_sinu[165,])

slope_mean_SHP <- map_slope_mean
slope_mean_SHP <- crop(slope_mean_SHP, shp_anp_sinu[165,])
slope_mean_SHP <- mask(slope_mean_SHP, shp_anp_sinu[165,])

meanMap <- slope_mean_SHP
meanMap[ pVal_mean_SHP > 0.05 ] <- NA

mpMean <- getMapview(meanMap, colPal = brewer.pal(256, "Reds"), 
                     nameLayer = "media", typeQuery = "click")

shp_mohinora_mp <- mapview(shp_anp_sinu[165,], layer.name="shp", 
                           legend=FALSE,
                           color="darkblue", lwd=4,
                           alpha.regions=0, homebutton=FALSE)

mpMean + shp_mohinora_mp

# --- 

pVal_annual_SHP <- map_pVal_annual
pVal_annual_SHP <- crop(pVal_annual_SHP, shp_anp_sinu[165,])
pVal_annual_SHP <- mask(pVal_annual_SHP, shp_anp_sinu[165,])

slope_annual_SHP <- map_slope_annual
slope_annual_SHP <- crop(slope_annual_SHP, shp_anp_sinu[165,])
slope_annual_SHP <- mask(slope_annual_SHP, shp_anp_sinu[165,])

annualMap <- slope_annual_SHP
annualMap[ pVal_annual_SHP > 0.05 ] <- NA

mpAnnual <- getMapview(annualMap, colPal = brewer.pal(256, "Greens"),
                     nameLayer = "anual", typeQuery = "click")

mpAnnual + shp_mohinora_mp

# ---

pVal_sannual_SHP <- map_pVal_sannual
pVal_sannual_SHP <- crop(pVal_sannual_SHP, shp_anp_sinu[165,])
pVal_sannual_SHP <- mask(pVal_sannual_SHP, shp_anp_sinu[165,])

slope_sannual_SHP <- map_slope_sannual
slope_sannual_SHP <- crop(slope_sannual_SHP, shp_anp_sinu[165,])
slope_sannual_SHP <- mask(slope_sannual_SHP, shp_anp_sinu[165,])

sannualMap <- slope_sannual_SHP
sannualMap[ pVal_sannual_SHP > 0.05 ] <- NA

mpSemiAnnual <- getMapview(sannualMap, colPal = brewer.pal(256, "Blues"),
                       nameLayer = "anual", typeQuery = "click")

mpSemiAnnual + shp_mohinora_mp




