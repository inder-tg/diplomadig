# --- Elaborado: Mar 8, 2024
# --- Actualizado: Mar 6, 2025; Abril 5, 2025
# --- DATASET: NDVI MOD13Q1 v061 en Cerro Mohinora, Chihuahua, 2000-2024

# --- Este script interpola linealmente todos los datos faltantes
# --- en DATASET después de aplicar capa de calidad, es decir, después de aplicar
# --- el código del script mohinora_QA.R

library(terra)
library(sf)

library(geoTS)
library(foreach)
library(doParallel)

install.packages("imputeTS", dependencies = TRUE)
library(imputeTS)
# library(raster)

source("Rscripts/auxFUN.R")

# ---

SHPfiles <- list.files(path = paste0( getwd(), "/data/outputs" ),
                       pattern = ".shp$",
                       full.names = TRUE)

mohinora_shp <- read_sf(SHPfiles[1])


FILES_NDVI_imputation <- list.files(path = paste0( getwd(), "/data/outputs/mohinora_imputation" ), #mohinora_NDVI_imputation[1],
                                    pattern = ".tif$",
                                    full.names = TRUE)

FILES_NDVI_QA <- list.files(path = paste0( getwd(), "/data/mohinora/250m_16_days_NDVI_QA" ),
                            pattern = ".tif$",
                            full.names = TRUE)

FILES_NDVI <- c(FILES_NDVI_imputation, FILES_NDVI_QA)

mohinora_DATA <- rast(FILES_NDVI) # raster::stack(FILES_NDVI)

mohinora_DATA_rTp <- spRast_valuesCoords(mohinora_DATA) # rasterToPoints(mohinora_DATA)

# mohinora_DATA_rTp_coords <- mohinora_DATA_rTp[,1:2]
# mohinora_DATA_rTp_values <- mohinora_DATA_rTp[,3:551]

# --- Ejemplo en un pixel

plot( subset(mohinora_DATA, 10) )
lines( mohinora_shp, lwd=4)

# -----------------------------------------------------------------------------
# --- Para analizar la serie de tiempo de cualquier píxel en la imagen sigue estos
# --- pasos:
# --- 1. Ejecuta la línea de abajo
XY <- locator()
# --- 2. Haz click (SOLO UNA VEZ) en algún píxel en la imagen
# --- 3. Presiona la tecla ESC de tu teclado
# --- 4. Continúa con el script a partir de la línea 107

xy <- get_timeSeries_byClicking(c(XY$x, XY$y),
                                df=mohinora_DATA_rTp$coords)

pixel <- mohinora_DATA_rTp$values[xy$coord, ]

# pixel <- mohinora_DATA_rTp$values[295, ]

pixel_ts <- ts(pixel, start = c(2000,1), end = c(2024,23),
               frequency = 23)

plot(pixel_ts, xlab="Años", ylab="NDVI", col="darkgreen", 
     main="pixel como objeto 'ts'")

pixel_interpol <- na_interpolation(pixel) 

pixel_interpol_ts <- ts(pixel_interpol, start = c(2000,1), end = c(2024,23), 
                        frequency = 23 )

plot(pixel_interpol_ts)


# CODIGO en paralelo

mohinora_interpol_linear <- matrix(nrow=nrow(mohinora_DATA_rTp$values), 
                                   ncol=ncol(mohinora_DATA_rTp$values))

# --- progress report file (to check out the process)

progressReportFile <- paste0(getwd(), "/RData/progressReports/mohinora/progress_temporal_gapfilling.txt" )
file.create(path=progressReportFile, showWarnings=FALSE)

write("===TEMPORAL GAPFILLING began at===",
      file=progressReportFile, append=TRUE)
write(as.character(Sys.time()[1]), file=progressReportFile,
      append=TRUE)

numCores <- detectCores()

kluster <- parallel::makeCluster(numCores-1, outfile="")
registerDoParallel(kluster)

output <- foreach(i=1:nrow(mohinora_DATA_rTp$values), .combine="rbind",
                  .packages="imputeTS") %dopar% { # nrow(sp_ndvi_rTp)
                    
                    pixel <- mohinora_DATA_rTp$values[i,]
                    
                    out_linear <- pixel
                    
                    if(length( is.na(pixel) ) > 0){
                      out_linear <- na_interpolation(pixel) 
                    }
                    
                    s <- c(as.numeric(out_linear))
                    
                    if(i %% 100 ==0){
                      texto <- paste0("Working on ROW: ", i)
                      write(texto, file=progressReportFile, append=TRUE)
                    }
                    
                    return(s)
                  }
stopCluster(kluster)

write("===TEMPORAL GAPFILLING ended at===",
      file=progressReportFile, append=TRUE)
write(as.character(Sys.time()[1]), file=progressReportFile,
      append=TRUE)

str(output)

mohinora_interpol_linear <- output

# --- rasterization

dirTIFS_toGet_names <- paste0(getwd(), "/data/mohinora/250m_16_days_NDVI")
listTIFnames <- list.files(path = dirTIFS_toGet_names,
                           pattern = ".tif$",
                           full.names = TRUE)

vectorNAMES <- character(length( listTIFnames ))
for(i in 1:length( listTIFnames )){
  temp <- listTIFnames[i]
  aux <- strsplit(temp, "/")
  basename <- aux[[1]][ length(aux[[1]]) ]
  nameBASE <- strsplit(basename, ".tif", fixed=TRUE)[[1]]
  vectorNAMES[i] <- paste0(nameBASE, 
                           "_interpol")
}


# --- asegurarse de crear /outputs/mohinora_interpolation
DIR_interpol <- paste0( getwd(), "/data/outputs/mohinora_interpolation" )
dir.create(DIR_interpol, recursive = TRUE)

# --- las primeras 3 columnas contienen valores de NDVI imputados
# --- a través del procedimiento de climatología, por tanto, no es necesario
# --- guardar esas capas nuevamente

PROJECTION <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs"

for(i in 4:ncol(mohinora_interpol_linear)){
  
  if( i %% 100 == 0){
    cat("Working on layer ", i, "\n")
  }
  
  mat <- cbind(mohinora_DATA_rTp$coords[,1:2], 
               mohinora_interpol_linear[,i])
  
  layer <- matrixToRaster(matrix=mat, 
                          projection=PROJECTION) 
  
  raster::writeRaster(x=layer,
                      filename = paste0(DIR_interpol,
                                        "/",
                                        vectorNAMES[i-3]),
                      format="GTiff",
                      datatype="INT4S",
                      overwrite=TRUE)
  
}

# --- verificando q todo está OK

maskTIF <- list.files(path=paste0(getwd(), "/data/mohinora/250m_16_days_NDVI_QA"),
                      pattern = ".tif$",
                      full.names = TRUE)

interpolTIFS <- list.files(path=paste0(getwd(), "/data/outputs/mohinora_interpolation"),
                           pattern = ".tif$",
                           full.names = TRUE)


a <- rast(maskTIF) # mask

b <- rast(interpolTIFS) # interpol

x <- 450

par(mfrow=c(1,2))
plot(subset(a,x), main="Sin interpolación")
# lines( mohinora_SHP_st, lwd=4)
plot(subset(b,x), main="Con interpolación")
# lines( mohinora_SHP_st, lwd=4)

