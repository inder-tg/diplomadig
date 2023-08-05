
# --- Diplomado Geomática, IG, UNAM, 2022
# --- Módulo VI Análisis de series de tiempo de imágenes satelitales con R

# --- Elaborado: Junio 2, 2022
# --- Impartido: Junio 25, 2022
# --- Actualizado: Agosto 3, 2023

# --- En este script presentamos un análisis de tendencias
# --- basado en los estadísticos no-paramétricos Mann-Kendall para tendencia
# --- y Theil-Sen para estimar pendiente e intercepto de tendencia lineal.
# --- Usamos código en paralelo para eficientar el cómputo a alta escala.

# --- La teoría detrás de las pruebas Mann-Kendall y Theil-Sen se puede
# --- encontrar en el archivo trendAlaysis.pdf en el folder PDF

# --- DATASET: NDVI MOD13Q1 en Cerro Mohinora, Chihuahua

# --- Preámbulo
library(raster)
library(mapview)
library(RColorBrewer)
library(gtools)
library(foreach)
library(doParallel)

library(Kendall)
library(trend)

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


# -----------------------------------------
# --- Análisis de tendencias: Exploración #
# -----------------------------------------

XY <- list(x=-10698785, y=2893289)

# xy <- get_timeSeries_byClicking(c(XY$x, XY$y),
#                                 df=mohinora_NDVI_rTp_full)

# pixel_full <- mohinora_NDVI_rTp_full[xy$coord, 3:ncol(mohinora_NDVI_rTp_full)]

pixel_full <- get_timeSeries_byClicking(c(XY$x, XY$y),
                                        df=mohinora_NDVI_rTp_full)

pixel_full_ts <- ts(as.numeric(pixel_full$ts) * 1e-4, 
                    start = c(2000,1), 
                    end = c(2009,23),
                    frequency = 23)

plot(pixel_full_ts)

pixel_MannKendall <- mk.test(as.numeric(pixel_full$ts) * 1e-4)

pixel_MannKendall$p.value

pixel_SenTheil <- sens.slope(as.numeric(pixel_full$ts) * 1e-4)

pixel_SenTheil$estimates

b_hat <- as.numeric(pixel_SenTheil$estimates)
a_hat <- median( as.numeric(pixel_full$ts) * 1e-4 - b_hat * 1:length(pixel_full) )

lineaTheilSen <- ts(a_hat + b_hat * 1:length(pixel_full), start = c(2000, 1),
                    end = c(2009, 23), frequency = 23)

par(mfrow=c(1,1), mar = c(2,2,1,2), adj=0)
plot(pixel_full_ts, type="l", col = "gray", ylab = "NDVI")
lines(lineaTheilSen, lwd = 5, col = "lightcoral")
legend("topright", legend = c("raw data", "linear trend"),
       col = c("gray", "lightcoral"), lty = rep(1,2), lwd = c(1,5), bty = "n")

# ---

# --------------------------------------------
# --- Tendencias: análisis a gran escala --- #
# --------------------------------------------

# --- CODIGO EN PARALELO

# df_pvalue guardará los pvalues de la prueba Mann-Kendall
# df_slope guardará los estimadores de pendiente Theil-Sen

df_pvalue <- matrix(nrow=nrow(mohinora_NDVI_rTp_full), ncol=3)
df_pvalue[,1:2] <- mohinora_NDVI_rTp_full[,1:2]

df_slope <- matrix(nrow=nrow(mohinora_NDVI_rTp_full), ncol=3)
df_slope[,1:2] <- mohinora_NDVI_rTp_full[,1:2]

numCores <- detectCores()

# --- progress report file (to check out on the process)
progressReportFile <- paste0(getwd(), "/RData/progressReports/mohinora_trendAnalysis.txt" )
file.create(path=progressReportFile, showWarnings=FALSE)

write("===TREND analysis began at===",
      file=progressReportFile, append=TRUE)
write(as.character(Sys.time()[1]), file=progressReportFile,
      append=TRUE)

kluster <- parallel::makeCluster(numCores-1, outfile="")
registerDoParallel(kluster)

output <- foreach(i=1:nrow(mohinora_NDVI_rTp_full), .combine="rbind",
                  .packages=c("trend") ) %dopar% { # nrow(sp_ndvi_rTp)
                    
                    pixel <- mohinora_NDVI_rTp_full[i, 3:ncol(mohinora_NDVI_rTp_full)] * 1e-4
                    
                    # pixel_interpol <- na_interpolation(pixel * 1e-4)
                    
                    pixel_MannKendall <- mk.test(pixel)$p.value
                    
                    pixel_SenTheil <- sens.slope(pixel)$estimates
                    
                    s <- c(as.numeric(pixel_MannKendall), as.numeric(pixel_SenTheil))
                    
                    if(i %% 100 ==0){
                      texto <- paste0("Working on ROW: ", i)
                      write(texto, file=progressReportFile, append=TRUE)
                    }
                    
                    return(s)
                  }
stopCluster(kluster)

write( as.character(Sys.time()[1]), file=progressReportFile, append=TRUE)
write( "===TREND analysis ended here===", file=progressReportFile, append=TRUE)
# ---

# --- Guardando los estadísticos como objetos matrix
df_pvalue[,3] <- output[,1]
df_slope[,3] <- output[,2]

# --- Asegurarse de crear /mohinora_trendAnalysis

save(df_pvalue, file=paste0(getwd(),"/RData/mohinora_trendAnalysis/pvalue.RData"))
save(df_slope, file=paste0(getwd(),"/RData/mohinora_trendAnalysis/slope.RData"))
# ---

# -----------------------
# --- RASTERIZACION --- #
# -----------------------

PROJECTION <- raster::projection(stack_primeras3Imagenes)

df_pvalue <- LoadToEnvironment(paste0(getwd(),"/RData/mohinora_trendAnalysis/pvalue.RData"))$df_pvalue
df_slope <- LoadToEnvironment(paste0(getwd(),"/RData/mohinora_trendAnalysis/slope.RData"))$df_slope

map_pvalue <- matrixToRaster_test(matrix=df_pvalue, projection=PROJECTION)
map_slope <- matrixToRaster_test(matrix=df_slope, projection=PROJECTION)

writeRaster(map_pvalue,
            filename = paste0( getwd(), "/data/mohinora_trendAnalysis/pvalueMap"),
            format="GTiff", datatype="FLT4S", overwrite=TRUE)

writeRaster(map_slope,
            filename = paste0( getwd(), "/data/mohinora_trendAnalysis/slopeMap"),
            format="GTiff", datatype="FLT4S", overwrite=TRUE)
# ---

# ------------------
# --- ANÁLISIS --- #
# ------------------

pvalueMap <- map_pvalue
pvalueMap[pvalueMap > 0.05] <- NA

slopeMap <- map_slope
slopeMap[is.na(pvalueMap)] <- NA

plot(map_slope, main="Pendiente estimada")

plot(slopeMap, main="Pendiente significativa estimada al 95% de significancia")

hist(slopeMap)

slopeValues <- getValues(slopeMap) # extrae los numeritos en slopeMap
slopeValues_sinNA <- slopeValues[!is.na(slopeValues)] # limpia 'slopeValues' de NA 

slope_densidad <- KernSmooth::bkde(x=slopeValues_sinNA) 

plot(slope_densidad, type = "l", main="Densidad de Pendiente", xlab="Pendiente",
     ylab="Frecuencia")
abline(v=0, col="blue", lty=3)

QUANT <- quantile(slopeMap, probs= seq(0,1,by=0.1) )

plot(slopeMap, col=c("red", "green"),
     breaks=c(QUANT[1], 0, QUANT[11]))
plot(shp_anp_sinu[165,], add=TRUE, lwd=4)


slopeValues_parte_negativa <- slopeValues_sinNA[slopeValues_sinNA<=0]
slopeValues_parte_positiva <- slopeValues_sinNA[slopeValues_sinNA>0]

plot(slope_densidad, type = "l", main="Densidad de Pendiente", xlab="Pendiente",
     ylab="Frecuencia")
abline(v=0, col="blue", lty=3)
abline(v=median(slopeValues_parte_negativa), col="red", lty=3)
abline(v=median(slopeValues_parte_positiva), col="red", lty=3)

QUANT_negativa <- quantile( slopeValues_parte_negativa, probs= seq(0, 1, by=0.1) )
QUANT_positiva <- quantile( slopeValues_parte_positiva, probs= seq(0, 1, by=0.05) )

plot(slope_densidad, type = "l", main="Densidad de Pendiente", xlab="Pendiente",
     ylab="Frecuencia")
abline(v=0, col="blue", lty=3)
abline(v=QUANT_negativa[11], col="orange", lwd=2, lty=2)
abline(v=QUANT_positiva[1], col="orange", lwd=2, lty=3)
legend("topright", legend = c(round(QUANT_negativa[10],6), 
                              round(QUANT_positiva[2],6)),
       lty=c(2,3), col=rep("orange", 2), lwd=2)

plot(slopeMap, col=c("red", "yellow", "green"),
     breaks=c(QUANT[1], QUANT_negativa[10], QUANT_positiva[2], QUANT[11]))
plot(shp_anp_sinu[165,], add=TRUE, lwd=4)

writeRaster(slopeMap,
            filename=paste0( getwd(), "/data/mohinora_trendAnalysis/slopeMap_sign"),
            format="GTiff", datatype="FLT4S", overwrite=TRUE)
# --- 

# -----------------------
# --- VISUALIZACION --- #
# -----------------------

slope_mp <- mapview(slopeMap, na.color = "transparent",
                    col.regions = c("red", "yellow", "green"),
                    at = c(QUANT[1], QUANT_negativa[10], QUANT_positiva[2], 
                           QUANT[11]))

shp_mohinora_mp <- mapview(shp_anp_sinu[165,], layer.name="shp", 
                           legend=FALSE,
                           color="darkblue", lwd=4,
                           alpha.regions=0, homebutton=FALSE)

slope_mp + shp_mohinora_mp


