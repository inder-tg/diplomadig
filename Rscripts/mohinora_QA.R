
# --- Elaborado Mar 9, 2023
# --- Código para calcular % de dato faltante -a nivel pixel-
# --- y maxGapLength en rasterStack
# --- DATASET: NDVI MOD13Q1 en Cerro Mohinora, Chihuahua, 2000-2023
# --- Actualizado: Feb 24, 2024, Abril 5. 2025
# --- DATASET: NDVI MOD13Q1 v061 en Cerro Mohinora, Chihuahua, 2000-2023 

# --- NOTA: Crear subdirectorios necesarios (/data/mohinora)

library(terra)
library(sf)

library(geoTS)
library(foreach)
library(doParallel)

source("Rscripts/auxFUN.R")

# ---

DIRS <-  list.dirs( path = paste0( getwd(), "/data"  ) ) #"C:/Users/inder/OneDrive/Desktop/proyectoCONABIO2025/mohinora"

NDVIfiles <- list.files(path = DIRS[4], # checar numero
                        pattern = ".tif", 
                        full.names = TRUE)

mohinora_DATA <- rast(NDVIfiles)

RELIABILITYfiles <- list.files(path = DIRS[5], # checar numero
                               pattern = ".tif", 
                               full.names = TRUE)

mohinora_DATA_reliability <- rast(RELIABILITYfiles)

# mohinoraRDataDIR <- paste0( mestiperDIR, "/RData" )

SHPfiles <- list.files(path = DIRS[7],
                       pattern = ".shp$",
                       full.names = TRUE)

mohinora_shp <- read_sf(SHPfiles[1])

# --- Ejemplos

TEMP <- subset(mohinora_DATA, 159)
AUX <- subset(mohinora_DATA_reliability, 159)

par(mfrow=c(1,2))
plot(TEMP)
lines(mohinora_shp)

plot(AUX)
lines(mohinora_shp, col = "cyan")

TEMP[ AUX >= 2 ] <- NA

plot(TEMP)
lines(mohinora_shp)

plot(AUX)
lines(mohinora_shp, col = "cyan")

# ---

whereToSave <- paste0(DIRS[3], "/250m_16_days_NDVI_QA") #paste0( getwd(), "/data/mohinora/250m_16_days_NDVI_QA" ) 
dir.create(whereToSave, recursive = TRUE)

TEMP <- subset(mohinora_DATA, 1)
AUX <- subset(mohinora_DATA_reliability, 1)
TEMP[ AUX >= 2 ] <- NA

nameFILE <- basename( NDVIfiles[1]  )
nameFILE <- paste0(strsplit( nameFILE, ".tif" )[[1]][1], "_QA.tif")
writeRaster(TEMP, 
            filename = paste0( whereToSave, "/", nameFILE ),
            datatype = datatype(mohinora_DATA)[1],
            overwrite = TRUE)


# mohinora_DATA_mask <- TEMP
for(i in 2:nlyr(mohinora_DATA)){
  TEMP <- subset(mohinora_DATA, i)
  AUX <- subset(mohinora_DATA_reliability, i)
  TEMP[ AUX >= 2 ] <- NA 
  # add(mohinora_DATA_mask) <- TEMP
  
  nameFILE <- basename( NDVIfiles[i]  )
  nameFILE <- paste0(strsplit( nameFILE, ".tif" )[[1]][1], "_QA.tif")
  writeRaster(TEMP, 
              filename = paste0( whereToSave, "/", nameFILE ),
              datatype = datatype(mohinora_DATA)[1],
              overwrite = TRUE)
}

# --- 

ndviQAFILES <- list.files( path = paste0( getwd(), "/data/mohinora/250m_16_days_NDVI_QA" ),
                           pattern = ".tif$",
                           full.names = TRUE )

mohinora_DATA_QA <- rast(ndviQAFILES)

mohinora_DATA_QA_shp <- crop(mohinora_DATA_QA, mohinora_shp,
                             mask=TRUE)

mohinora_DATA_QA_rTp <- spRast_valuesCoords(mohinora_DATA_QA_shp)

# --- EJEMPL sobre un pixel

pixel <- mohinora_DATA_QA_rTp$values[1700,]

(pixel_percentMiss <- sum(is.na(pixel)) / length(pixel)) * 100 # length de cualquier pixel es 483

(pixel_maxgap <- maxLagMissVal(x=pixel)$maxLag)

# --- COMPUTO en PARALELO

df_perc_miss <- matrix(nrow=nrow(mohinora_DATA_QA_rTp$values), ncol=3)
df_perc_miss[,1:2] <- mohinora_DATA_QA_rTp$coords[,1:2]

maxgap_df <- matrix(nrow=nrow(mohinora_DATA_QA_rTp$values), ncol=3)
maxgap_df[,1:2] <- mohinora_DATA_QA_rTp$coords[,1:2]

# --- progress report file (to check out on the process)
# --- antes de ejecutar, crear /RData/progressReports/mohinora (sólo en caso de que los directorios no existan)

dir.create( paste0(getwd(), "/RData/progressReports/mohinora"), recursive = TRUE )

numCores <- detectCores()

progressReportFile <- paste0(getwd(), "/RData/progressReports/mohinora/progress_QA.txt" )
file.create(path=progressReportFile, showWarnings=FALSE)

write("===QA analysis began at===",
      file=progressReportFile, append=TRUE)
write(as.character(Sys.time()[1]), file=progressReportFile,
      append=TRUE)

kluster <- parallel::makeCluster(numCores-1, outfile="")
registerDoParallel(kluster)

output <- foreach(i=1:nrow(mohinora_DATA_QA_rTp$values), .combine="rbind",
                  .packages="geoTS") %dopar% { 
                    
                    pixel <- mohinora_DATA_QA_rTp$values[i,]
                    
                    pixel_percentMiss <- sum(is.na(pixel)) / length(pixel) # length de cualquier pixel es 483
                    
                    pixel_maxgap <- maxLagMissVal(x=pixel)$maxLag
                    
                    s <- c(as.numeric(pixel_percentMiss), as.numeric(pixel_maxgap))
                    
                    if(i %% 100 ==0){
                      texto <- paste0("Working on ROW: ", i)
                      write(texto, file=progressReportFile, append=TRUE)
                    }
                    
                    return(s)
                  }
stopCluster(kluster)

write( as.character(Sys.time()[1]), file=progressReportFile, append=TRUE )
write( "===QA analysis ended here===", file=progressReportFile, append=TRUE )

# ---

# --- Saving output

df_perc_miss[,3] <- output[,1]
maxgap_df[,3] <- output[,2]

# --- antes de ejecutar, crear /RData/mohinora_QA

dir.create( paste0(getwd(),"/RData/mohinora_QA"), recursive = TRUE )

save(df_perc_miss, file=paste0(getwd(),"/RData/mohinora_QA/percent_missingValue.RData"))
save(maxgap_df, file=paste0(getwd(),"/RData/mohinora_QA/maxgap.RData"))

# --- Rasterization

PROJECTION <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R=6371007.181 +units=m +no_defs" #projection(mohinora_DATA) # crs(mohinora_mask) # raster::projection(TEMP) # crs(STACK_sp_ndvi_subset)

map_percentMissing <- matrixToRaster(matrix=df_perc_miss, projection=PROJECTION)
map_maxgap <- matrixToRaster(matrix=maxgap_df, projection=PROJECTION)

map_percentMissing
map_maxgap

# --- antes de ejecutar, crear /outputs/mohinora_QA

dir.create( paste0( getwd(), "/data/outputs/mohinora_QA" ),
            recursive = TRUE )

writeRaster(map_percentMissing,
            filename = paste0( getwd(), "/data/outputs/mohinora_QA/missingValue" ),
            format="GTiff", datatype="FLT4S", overwrite=TRUE)

writeRaster(map_maxgap,
            filename = paste0( getwd(), "/data/outputs/mohinora_QA/maxGap" ),
            format="GTiff", datatype="INT2U", overwrite=TRUE)

# ---

QAfiles <- list.files(path=paste0(getwd(), "/data/outputs/mohinora_QA"),
                      pattern=".tif",
                      full.names=TRUE)

maxGap <- rast(QAfiles[1])
percent <- rast(QAfiles[2])

maxGap <- crop(maxGap, mohinora_shp, mask=TRUE)
percent <- crop(percent, mohinora_shp, mask=TRUE)

par(mfrow=c(1,2))
plot(maxGap, main="max-gap length")
plot(percent * 100, main="% missing values")

# ---
