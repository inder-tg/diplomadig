# --- Diplomado Geomática, IG, UNAM, 2022
# --- Módulo VI Análisis de series de tiempo de imágenes satelitales con R

# --- Elaborado: Junio 2, 2022
# --- Impartido: Julio 1, 2022
# --- ACtualizado: Agosto 4, 2023

# --- En este script presentamos un análisis de puntos de cambio basado en 
# --- rutinas del paquete bfast.
# --- Usamos código en paralelo para eficientar el cómputo a alta escala.

# --- DATASET: NDVI MOD13Q1 en Cerro Mohinora, Chihuahua

# --- Preámbulo
library(raster)
library(rasterVis)
library(mapview)
library(RColorBrewer)
library(gtools)
library(foreach)
library(doParallel)

library(bfast)
library(tmap)

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


# -------------------------------------------------
# --- Análisis de cambio abrupto: Exploración --- #
# -------------------------------------------------

# --- Clase de tendencia: 5; Año: 2002

XY <- list(x=-10698785, y=2893289)


pixel_full_type5 <- get_timeSeries_byClicking(c(XY$x, XY$y),
                                              df=mohinora_NDVI_rTp_full)

pixel_ts <- ts(as.numeric(pixel_full_type5$ts), 
               start = c(2000,1), 
               end = c(2009, 23),
               frequency = 23)

pixel_bfast <- bfast(pixel_ts, season = "harmonic", h=0.2)

plot(pixel_bfast)

pixel_bfast01 <- bfast01(pixel_ts, bandwidth=0.2)

plot(pixel_bfast01)

bfast01classify(pixel_bfast01)

(temp <- getBreak(data = pixel_ts, start=2000, end=2009, bw=0.2))

getYear(start = 2000, end = 2009, bp=temp$bPs)

# --- Clase de tendencia: 1

XY <- list(x=--10689992, y=2891202)

pixel_full_type1 <- get_timeSeries_byClicking(c(XY$x, XY$y),
                                              df=mohinora_NDVI_rTp_full)

pixel_ts <- ts(as.numeric(pixel_full_type1$ts), 
               start = c(2000,1), end = c(2009, 23),
               frequency = 23)

pixel_bfast01 <- bfast01(pixel_ts, bandwidth=0.2)

plot(pixel_bfast01)

bfast01classify(pixel_bfast01)

(temp <- getBreak(data = pixel_full_type1, start=2000, end=2009, bw=0.2))

getYear(start = 2000, end = 2009, bp=temp$bPs)


# --- Clase de tendencia: 5; Año: 2006

XY <- list(x=--10698418, y=2887497)

pixel_full_type5_2006 <- get_timeSeries_byClicking(c(XY$x, XY$y),
                                                   df=mohinora_NDVI_rTp_full)

pixel_ts <- ts(as.numeric(pixel_full_type5_2006), 
               start = c(2000,1), end = c(2009, 23),
               frequency = 23)

pixel_bfast01 <- bfast01(pixel_ts, bandwidth = 0.2)

plot(pixel_bfast01)

bfast01classify(pixel_bfast01)

(temp <- getBreak(data = pixel_full_type5_2006, start=2000, end=2009, bw=0.2))

getYear(start = 2000, end = 2009, bp=temp$bPs)


# --- Clase de tendencia: 7

XY <- list(x=-10687435, y=2894415)

pixel_full <- get_timeSeries_byClicking(c(XY$x, XY$y),
                                        df=mohinora_NDVI_rTp_full)

pixel_ts <- ts(as.numeric(pixel_full$ts), start = c(2000,1), end = c(2009, 23),
               frequency = 23)

pixel_bfast01 <- bfast01(pixel_ts, bandwidth = 0.2)

plot(pixel_bfast01)

bfast01classify(pixel_bfast01)

(temp <- getBreak(data = pixel_full, start=2000, end=2009, bw=0.2))

getYear(start = 2000, end = 2009, bp=temp$bPs)

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

# -----------------------------------------------
# --- Cambio abrupto: cómputo a gran escala --- #
# -----------------------------------------------

# --- CODIGO EN PARALELO

df_YEARS <- matrix(nrow=nrow(mohinora_NDVI_rTp_full), ncol=3)
df_YEARS[,1:2] <- mohinora_NDVI_rTp_full[,1:2]

df_CP <- matrix(nrow=nrow(mohinora_NDVI_rTp_full), ncol=3)
df_CP[,1:2] <- mohinora_NDVI_rTp_full[,1:2]

df_TYPE <- matrix(nrow=nrow(mohinora_NDVI_rTp_full), ncol=3)
df_TYPE[,1:2] <- mohinora_NDVI_rTp_full[,1:2]

df_SIGN <- matrix(nrow=nrow(mohinora_NDVI_rTp_full), ncol=3)
df_SIGN[,1:2] <- mohinora_NDVI_rTp_full[,1:2]

df_STABLE <- matrix(nrow=nrow(mohinora_NDVI_rTp_full), ncol=3)
df_STABLE[,1:2] <- mohinora_NDVI_rTp_full[,1:2]
# ---
numCores <- detectCores()

# --- progress report file (to check out on the process)
progressReportFile <- paste0(getwd(), "/RData/progressReports/mohinora_cps.txt" )
file.create(path=progressReportFile, showWarnings=FALSE)

write("===CHANGE-POINT analysis began at===",
      file=progressReportFile, append=TRUE)
write(as.character(Sys.time()[1]), file=progressReportFile,
      append=TRUE)

kluster <- parallel::makeCluster(numCores-1, outfile="")
registerDoParallel(kluster)

# nrow(mohinora_NDVI_rTp_full)

output <- foreach(j=1:nrow(mohinora_NDVI_rTp_full), .combine="rbind",
                  .packages=c("bfast", "imputeTS")) %dopar% {
                    # ---
                    pixel <- mohinora_NDVI_rTp_full[j,3:ncol(mohinora_NDVI_rTp_full)]
                    
                    # ---
                    BREAKS <- getBreak(data=pixel, start=2000, end=2009, bw=0.2)
                    
                    YEARS <- getYear(bp=BREAKS$bPs, start=2000, end=2009)
                    # ---
                    s <- c(YEARS, BREAKS$bPs, BREAKS$type, 
                           BREAKS$significance, BREAKS$stability)
                    
                    if(j %% 250==0){
                      text <- paste0("Working on ROW: ", j)
                      write(text, file=progressReportFile, append=TRUE)
                    }
                    
                    # if(i %% 100==0){
                    # }
                    
                    return(s)
                  }
stopCluster(kluster)
# --- END parallel processing
write( as.character(Sys.time()[1]), file=progressReportFile, append=TRUE)
write( "===CHANGE-POINT analysis ended here===", file=progressReportFile, 
       append=TRUE)
# ---

# --- Guardando el output de bfast01 como objetos matrix
df_YEARS[,3] <- output[,1]
df_CP[,3] <- output[,2]
df_TYPE[,3] <- output[,3]
df_SIGN[,3] <- output[,4]
df_STABLE[,3] <- output[,5]


# --- Asegurarse de haber creado el folder /mohinora_cps

save(df_YEARS, file=paste0(getwd(),"/RData/mohinora_cps/YEARS.RData"))
save(df_CP, file=paste0(getwd(),"/RData/mohinora_cps/CP.RData"))
save(df_TYPE, file=paste0(getwd(),"/RData/mohinora_cps/TYPE.RData"))
save(df_SIGN, file=paste0(getwd(),"/RData/mohinora_cps/SIGN.RData"))
save(df_STABLE, file=paste0(getwd(),"/RData/mohinora_cps/STABLE.RData"))
# --- 

# -----------------------
# --- RASTERIZACION --- #
# -----------------------

PROJECTION <- raster::projection(stack_primeras3Imagenes)

# df_YEARS <- LoadToEnvironment(paste0(getwd(),"/RData/mohinora_cps/YEARS.RData"))$df_YEARS
# df_CP <- LoadToEnvironment(paste0(getwd(),"/RData/mohinora_cps/CP.RData"))$df_CP
# df_TYPE <- LoadToEnvironment(paste0(getwd(),"/RData/mohinora_cps/TYPE.RData"))$df_TYPE
# df_SIGN <- LoadToEnvironment(paste0(getwd(),"/RData/mohinora_cps/TYPE.RData"))$df_SIGN
# df_STABLE <- LoadToEnvironment(paste0(getwd(),"/RData/mohinora_cps/TYPE.RData"))$df_STABLE

map_YEARS <- matrixToRaster_test(matrix=df_YEARS, projection=PROJECTION)
map_CP <- matrixToRaster_test(matrix=df_CP, projection=PROJECTION)
map_TYPE <- matrixToRaster_test(matrix=df_TYPE, projection=PROJECTION)
map_SIGN <- matrixToRaster_test(matrix=df_SIGN, projection=PROJECTION)
map_STABLE <- matrixToRaster_test(matrix=df_STABLE, projection=PROJECTION)


# --- Asegurarse de haber creado /mohinora_cps

writeRaster(map_YEARS,
            filename = paste0( getwd(), "/data/mohinora_cps/map_YEARS"),
            format="GTiff", datatype="INT2U", overwrite=TRUE)

writeRaster(map_CP,
            filename = paste0( getwd(), "/data/mohinora_cps/map_CP"),
            format="GTiff", datatype="INT1U", overwrite=TRUE)

writeRaster(map_TYPE,
            filename = paste0( getwd(), "/data/mohinora_cps/map_TYPE"),
            format="GTiff", datatype="INT1U", overwrite=TRUE)

writeRaster(map_SIGN,
            filename = paste0( getwd(), "/data/mohinora_cps/map_SIGN"),
            format="GTiff", datatype="FLT4S", overwrite=TRUE)

writeRaster(map_STABLE,
            filename = paste0( getwd(), "/data/mohinora_cps/map_STABLE"),
            format="GTiff", datatype="FLT4S", overwrite=TRUE)

# ----------------
# --- ANALISIS ---
# ----------------

# --- TYPES
# --- 1 - Enverdecimiento
# --- 2 - Pardeamiento
# --- 3 - Enver. sostenido
# --- 4 - Pardea. sostenido
# --- 5 - Enver. demorado
# --- 6 - Pardea. demorado
# --- 7 - Enver. a Pardea.
# --- 8 - Pardea. a Enver.

hist(map_TYPE) 

# TYPE 4 no está presente
# TYPE1 > TYPE5 > TYPE7 > TYPE2 > TYPE6 > TYPE3 > TYPE8

TYPE_SHP <- map_TYPE
TYPE_SHP <- crop(TYPE_SHP, shp_anp_sinu[165,])
TYPE_SHP <- mask(TYPE_SHP, shp_anp_sinu[165,])

plot(map_TYPE)

plot(TYPE_SHP)

hist(TYPE_SHP) 

# --- SIGN
# --- 0 - Ambos segmentos significativos o sin break y tendencia significativa
# --- 1 - Sólo el primer segmento significativo
# --- 2 - Sólo el segundo segmento significativo
# --- 3 - Ambos segmentos no significativos o sin break y no significativo

SIGN_SHP <- map_SIGN
SIGN_SHP <- crop(SIGN_SHP, shp_anp_sinu[165,])
SIGN_SHP <- mask(SIGN_SHP, shp_anp_sinu[165,])

YEARS_SHP <- map_YEARS
YEARS_SHP <- crop(YEARS_SHP, shp_anp_sinu[165,])
YEARS_SHP <- mask(YEARS_SHP, shp_anp_sinu[165,])

hist(SIGN_SHP)

hist(YEARS_SHP)

# --- TYPE1_SIGN_YEARS

TYPE1_SIGN <- SIGN_SHP
TYPE1_SIGN[TYPE_SHP != 1] <- NA

hist(TYPE1_SIGN)

TYPE1_SIGN[ SIGN_SHP == 3 ] <- NA

plot(TYPE1_SIGN)

# ---

YEARS_TYPE1_SIGN <- YEARS_SHP
YEARS_TYPE1_SIGN[is.na(TYPE1_SIGN)] <- NA

plot(YEARS_TYPE1_SIGN)

hist(YEARS_TYPE1_SIGN)

# --- TYPE5_SIGN_YEARS

TYPE5_SIGN <- SIGN_SHP
TYPE5_SIGN[TYPE_SHP != 5] <- NA

hist(TYPE5_SIGN)

TYPE5_SIGN[ SIGN_SHP == 1] <- NA

plot(TYPE5_SIGN)

# ---

YEARS_TYPE5_SIGN <- YEARS_SHP
YEARS_TYPE5_SIGN[is.na(TYPE5_SIGN)] <- NA

plot(YEARS_TYPE5_SIGN)

hist(YEARS_TYPE5_SIGN)

# -----------------------
# --- VISUALIZACION --- #
# -----------------------

COLORES <- c("#50C878", "#A52A2A",
             "#006B3C", "#E32636",
             "#87A96B", "#9F8170",
             # "#FBEC5D", 
             "orange",
             "#66FF00")

LABELS <- c("Enverdecimiento", # 1
            "Pardeamiento",    # 2
            "E. sostenido",    # 3
            "P. sostenido",    # 4
            "E. demorado",     # 5
            "P. demorado",     # 6
            "E. a P.",         # 7
            "P. a E.")         # 8

# ---
  
shp_mohinora <- tm_shape(shp_anp_sinu[165,]) +
  tm_borders(lwd = 2) + 
  tm_layout(main.title = "", 
            # outer.margins = c(B,L,T,R)
            outer.margins = c(0.015,0.015,0,0.015),
            inner.margins = c(0,0.015,0,0.45),
            main.title.size = 1)

TYPE1 <- TYPE_SHP
TYPE1[ is.na(TYPE1_SIGN) ] <- NA

map_TYPE1_SIGN <- tm_shape(TYPE1) + 
  tm_raster(style = "fixed", 
            title = "",
            breaks=c(1:2),
            palette = COLORES[1],
            labels = LABELS[1])

shp_mohinora + map_TYPE1_SIGN 

# ---

TYPE5 <- TYPE_SHP
TYPE5[ is.na(TYPE5_SIGN) ] <- NA

map_TYPE5_SIGN <- tm_shape(TYPE5) + 
  tm_raster(style = "fixed", 
            title = "",
            breaks=c(5:6),
            palette = COLORES[5],
            labels = LABELS[5])

shp_mohinora + map_TYPE5_SIGN 


shp_mohinora + map_TYPE1_SIGN + map_TYPE5_SIGN 
