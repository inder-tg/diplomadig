# --- Diplomado Geomática, IG, UNAM, 2023
# --- Módulo X: Análisis de series de tiempo de imágenes satelitales con R

# --- Elaborado: Agosto 8, 2023
# --- Impartido: Agosto 11, 2023

# --- En este script hacemos:
# --- Análisis de tendencias (Mann-Kendall test y Theil-Sen estimador de pendiente)
# --- Análisis de cambio abrupto (bfast01 y bfast)

# --- La teoría detrás de las pruebas Mann-Kendall y Theil-Sen se puede
# --- encontrar en el archivo trendAlaysis.pdf en el folder PDF

# --- DATASET: LP_NDMI_S2_2019_2022.tif pegado en la nube del diplomado

library(terra)
library(rgdal)
library(raster)
library(mapview)
library(bfast)
library(Kendall)
library(trend)

source( paste0( getwd(), "/Rscripts/auxFUN.R" ) )

# ---

tifFILES <- list.files(path = paste0( getwd(), "/data/LaPiedad" ),
                       pattern = ".tif",
                       full.names = TRUE)

LPsentinel <- stack(tifFILES)

LPsentinel_rTp <- rasterToPoints(LPsentinel)

sentinel_LP <- rast(tifFILES)

# LPsentinel_rTp <- spRast_valueCoords(spRaster=LPsentinel, 
#                                      na_rm=TRUE)

# ---

plot(subset(LPsentinel,25))

# Agricola
# 20.365365, -102.024281

xy <- SpatialPoints(cbind(-102.024281, 20.365365))
proj4string(xy) <- CRS('+init=epsg:4326')
xy_UTM <- spTransform(xy, CRS('+init=epsg:32613'))

plot(xy_UTM, add=TRUE, pch="M")

# Agave
# 20.363361, -102.034041

xyT <- SpatialPoints(cbind(-102.034041, 20.363361))
proj4string(xyT) <- CRS('+init=epsg:4326')
xyT_UTM <- spTransform(xyT, CRS('+init=epsg:32613'))

plot(xyT_UTM, add=TRUE, pch="A")

# ---

mp <- mapview::mapview(subset(LPsentinel, 25))
mp_maiz <- mapview(xy_UTM)
mp_agave <- mapview(xyT_UTM)

mp + mp_maiz + mp_agave

# ---
sentinel_LP <- rast(tifFILES)
xy_UTM_vect <- vect(xy_UTM, 
                    crs="+proj=utm +zone=13 +datum=WGS84")

pixel_agro <- raster::extract(x=LPsentinel, y=xy_UTM)
agro_pixel <- terra::extract(x=sentinel_LP, y=xy_UTM_vect)

# --- comparación pixel_agro vs agro_pixel
cbind(as.numeric(pixel_agro), as.numeric(agro_pixel)[-1])

pixel_agro_ts <- ts(as.numeric(agro_pixel)[-1],
                    start = c(2019, 1),
                    end = c(2022, 12),
                    frequency = 12)

plot(pixel_agro_ts, ylab="NDMI", main="Producto sin escalar")

pixel_agro_clima <- climatology(x=as.numeric(agro_pixel)[-1],
                                lenPeriod = 12)

str(pixel_agro_clima)

boxplot(pixel_agro_clima$matrix)

agro_pixel_filled <- gapfill_climatology(y=as.numeric(agro_pixel)[-1],
                                         box="median", gapType = NaN,
                                         lenPeriod = 12)

filled_pixel <- climatology(x=agro_pixel_filled$filled,
                            lenPeriod = 12)

boxplot(filled_pixel$matrix)

agro_pixel_filled_ts <- ts(as.numeric(agro_pixel_filled$filled),
                           start = c(2019, 1),
                           end = c(2022, 12),
                           frequency = 12)

plot(agro_pixel_filled_ts, ylab="NDMI", 
     main=paste0("lon.: ", round(extent(xy)[1],2), 
                 "  lat.: ", round(extent(xy)[3],2)))

# --- anomalias
mu <- apply(filled_pixel$matrix, 2, mean)

sigma <- apply(filled_pixel$matrix , 2, sd)

anomalia <- (filled_pixel$matrix - mu)/sigma

anomalia_ts <- ts(c(t(anomalia)),
                  start = c(2019, 1),
                  end = c(2022, 12), 
                  frequency = 12)

plot(anomalia_ts, ylab="", 
     main=paste0("lon.: ", round(extent(xy)[1],2), 
                 "  lat.: ", round(extent(xy)[3],2)))
abline(h=2, col="darkgreen")
abline(h=-2, col="darkgreen")

# --- tendencia
pixel_MK <- mk.test(as.numeric(agro_pixel_filled$filled) * 1e-4)

pixel_MK$p.value

pixel_SenTheil <- sens.slope(as.numeric(agro_pixel_filled$filled) * 1e-4)

pixel_SenTheil$estimates

b_hat <- as.numeric(pixel_SenTheil$estimates)
a_hat <- median( as.numeric(agro_pixel_filled$filled) * 1e-4 - b_hat * 1:length(agro_pixel_filled$filled) )

lineaTheilSen <- ts(a_hat + b_hat * 1:length(agro_pixel_filled$filled), 
                    start = c(2019, 1),
                    end = c(2022, 12), frequency = 12)

par(mfrow=c(1,1), mar = c(2,2,1,2), adj=0)
plot(pixel_agro_ts * 1e-4, type="l", col = "gray", ylab = "NDMI")
lines(lineaTheilSen, lwd = 5, col = "lightcoral")
legend("topright", legend = c("raw data", "linear trend"),
       col = c("gray", "lightcoral"), lty = rep(1,2), lwd = c(1,5), bty = "n")

# --- cambio abrupto

agro_bp <- bfast(agro_pixel_filled_ts,
                 season = "harmonic")

plot(agro_bp)

# --- AGAVE

xyT_UTM_vect <- vect(xyT_UTM, 
                     crs="+proj=utm +zone=13 +datum=WGS84")

agave_pixel <- terra::extract(x=sentinel_LP, y=xyT_UTM_vect)

agave_pixel_ts <- ts(as.numeric(agave_pixel)[-1],
                     start = c(2019, 1),
                     end = c(2022, 12),
                     frequency = 12)

par(mfrow=c(1,1), mar = c(5,4,4,2))
plot(agave_pixel_ts)

agave_pixel_clima <- climatology(x=as.numeric(agave_pixel),
                                 lenPeriod = 12)

boxplot(agave_pixel_clima$matrix)

# --- anomalias

mu <- apply(agave_pixel_clima$matrix, 2, mean)

sigma <- apply(agave_pixel_clima$matrix , 2, sd)

anomalia <- (agave_pixel_clima$matrix - mu)/sigma

anomalia_ts <- ts(c(t(anomalia)),
                  start = c(2019, 1),
                  end = c(2022, 12), 
                  frequency = 12)

plot(anomalia_ts)
abline(h=2, col="darkgreen")
abline(h=-2, col="darkgreen")

# --- tendencia
pixel_MK <- mk.test(as.numeric(agave_pixel)[-1] * 1e-4)

pixel_MK$p.value

pixel_SenTheil <- sens.slope(as.numeric(agave_pixel)[-1] * 1e-4)

pixel_SenTheil$estimates

b_hat <- as.numeric(pixel_SenTheil$estimates)
a_hat <- median( as.numeric(agave_pixel)[-1] * 1e-4 - b_hat * 1:length(agave_pixel[-1]) )

lineaTheilSen <- ts(a_hat + b_hat * 1:length(agave_pixel[-1]), 
                    start = c(2019, 1),
                    end = c(2022, 12), frequency = 12)

par(mfrow=c(1,1), mar = c(2,2,1,2), adj=0)
plot(agave_pixel_ts * 1e-4, type="l", col = "gray", ylab = "NDVI")
lines(lineaTheilSen, lwd = 5, col = "lightcoral")
legend("topright", legend = c("raw data", "linear trend"),
       col = c("gray", "lightcoral"), lty = rep(1,2), lwd = c(1,5), bty = "n")

# --- cambio abrupto

agave_bp <- bfast(agave_pixel_ts,
                  season = "harmonic")

plot(agave_bp)

agave_bp01 <- bfast01(agave_pixel_ts, 
                      bandwidth = .15, 
                      order = 1)

plot(agave_bp01)

# --- pixel_tsAnalysis

x <- -102.0215001 #-102.0366292 #
y <- 20.3654926 # 20.363994 #

XY <- SpatialPoints(cbind(x, y))
proj4string(XY) <- CRS('+init=epsg:4326')
XY_UTM <- spTransform(XY, CRS('+init=epsg:32613'))

XY_UTM_vect <- vect(XY_UTM, 
                    crs="+proj=utm +zone=13 +datum=WGS84")

pixel <- terra::extract(x=sentinel_LP, y=XY_UTM_vect)

XY_list <- list(x=810948.1, y=2254742)
# XY_list <- list(x=extent(XY_UTM)[1], y=extent(XY_UTM)[3])

pixel_homemade <- get_timeSeries_byClicking(toPlot = XY_list, 
                                            df = LPsentinel_rTp)

# cbind(as.numeric(pixel_homemade$ts),
#       as.numeric(pixel)[-1])

pixel_ts <- ts(as.numeric(pixel_homemade$ts),
               start = c(2019, 1),
               end = c(2022, 12),
               frequency = 12)

pixel_bp01 <- bfast01(pixel_ts)

pixel_bp01 <- bfast01(pixel_ts, order = 1)

pixel_bps <- bfast(pixel_ts, season = "harmonic")

# --- visualización

par(mfrow=c(2,2))

plot(subset(LPsentinel,25))
plot(XY_UTM, add=TRUE, pch="N")

plot(pixel_ts, ylab="NDMI", 
     main=paste0("lon.: ", round(extent(XY)[1],2), 
                 "  lat.: ", round(extent(XY)[3],2)))

trendAnalysis(x=as.numeric(pixel_homemade$ts) * 1e-4,
              startYear = 2019,
              endYear = 2022,
              frequency = 12,
              productName = "NDWI")

plot(pixel_bp01)

plot(pixel_bps)
