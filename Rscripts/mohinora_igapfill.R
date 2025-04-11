# --- Diplomado Geomática, IG, UNAM, 2025
# --- Bloque 2: Sistemas de Información Geográfica
# --- Módulo V: R como herramienta de SIG
# --- ELABORADO: Abril 11, 2025
# --- DATASET: NDVI MOD13Q1 en Cerro Mohinora, Chihuahua, 2000-2024
# --- EN este script mostramos el potencial del paquete igapfill
# --- para rellenar huecos en STIS utilizando una estrategia espacio-temporal
# --- EL paquete igapfill se encuentra en etapa de pruebas

# --- Preámbulo END

library(igapfill) # este paquete hay que instalarlo localmente
library(heatmaply)

# --- Creando materiales in /250m_16_days_NDVI_QA_byYears

NDVIQAdir <- paste0( getwd(), "/data/mohinora/250m_16_days_NDVI_QA" )

FILES_NDVI_imputation <- list.files(path = paste0( getwd(), "/data/outputs/mohinora_imputation" ), #mohinora_NDVI_imputation[1],
                                    pattern = ".tif$",
                                    full.names = TRUE)

FILES_NDVI_QA <- list.files(path = NDVIQAdir,
                            pattern = ".tif$",
                            full.names = TRUE)

FILES_NDVI <- c(FILES_NDVI_imputation, FILES_NDVI_QA)

YEARS <- 2000:2024

for (i in 1:length(YEARS)) {
  
  if(!dir.exists( paste0( NDVIQAdir, "_byYears/", YEARS[i] ) )){
    dir.create( paste0( NDVIQAdir, "_byYears/", YEARS[i] ), 
                recursive = TRUE )
  }
  
  file.copy(from = FILES_NDVI[((i-1)*23+1):(i*23)], 
            to = paste0( NDVIQAdir, "_byYears/", YEARS[i] ) )
  
}

# --- First iteration

NDVIbyYearsDIR <- paste0( NDVIQAdir, "_byYears")

DIRS <- list.dirs(path = NDVIbyYearsDIR)[-1]

COLNAMES <- paste0( "DoY-", seq(1,365,by=16) )

missingValueSieve <- mvSieve(dirs=DIRS,
                             filesPerDir = 23,
                             startPeriod = 2000,
                             endPeriod = 2024,
                             colNames = COLNAMES)

heatmaply(missingValueSieve/(54*88) * 100, 
          limits = c(0,100),
          colors = cool_warm,
          dendrogram = "none",
          xlab = "", ylab = "",
          main = "% of missing values in NDVI TS of Cerro Mohinora (2000-2024)",
          scale = "none",
          draw_cellnote=TRUE,
          cellnote_textposition = "middle center",
          cellnote_size = 10,
          margins = c(60, 100, 40, 20),
          grid_color = "white",
          # grid_width = 1e-26,
          # grid_size = 1e-26,
          titleX = FALSE,
          hide_colorbar = TRUE,
          labRow = rownames(missingValueSieve),
          labCol = colnames(missingValueSieve),
          heatmap_layers = theme(axis.line = element_blank()))

# --- BLOCK0

SIEVE <- (missingValueSieve/(54*88)) * 100

sieveMinBlock(sieve=SIEVE, rank=10)

# --- gapfilliing block0

DIRS <- list.dirs(path = NDVIbyYearsDIR)[-1]

files2019 <- list.files(path = DIRS[20],
                        pattern = ".tif$",
                        full.names = TRUE)

files2020 <- list.files(path = DIRS[21],
                        pattern = ".tif$",
                        full.names = TRUE)

block0DIR <- paste0( getwd(), "/data/outputs/mohinora_igapfill/block0" )

if( !dir.exists(block0DIR) ){
  dir.create( block0DIR, recursive = TRUE )
}

file.copy(from=files2019[10:11], to=block0DIR)

file.copy(from=files2020[10:11], to=block0DIR)

# --- C:/Users/inder/OneDrive/Desktop/teaching/diplomado/diplomado2025/diplomadig/data/outputs/mohinora_igapfill/block0

igapfill()

# --- backup code si se presenta algún error en igapfill()

allDIRS <- list.dirs(path=block0DIR, full.names = TRUE)

block0FILES <- list.files(path = block0DIR,
                          pattern = ".tif$",
                          full.names = TRUE)

LAT <- get_LAT(stack=raster::stack(block0FILES))

LON <- get_LON(stack=raster::stack(block0FILES))

applyGapfill(inputDir = allDIRS[7],
             outputDir = allDIRS[5],
             progressDir = allDIRS[6],
             lat = LAT,
             lon = LON,
             days = c(145, 161),
             years = 2019:2020,
             numCores = 4,
             scale = 1e0,
             clipRange = c(-1e4, 1e4))

parallel_mosaic(inputDirImages = block0DIR,
                inputDirRData = allDIRS[5],
                inputDirMaster = allDIRS[4],
                outputDir = allDIRS[3],
                progressReportDir = allDIRS[6],
                scaleFactor = 1e0,
                dataType = "INT4S",
                numCores = 2)

# --- checkin before and after gapfill

originalImages <- list.files(path = block0DIR,
                             pattern = ".tif$",
                             full.names = TRUE)

gapfillImages <- list.files(path = paste0(block0DIR, "/gapfill/filled"),
                            pattern = ".tif$",
                            full.names = TRUE)

originalRasters <- rast(originalImages)
gapfillRasters <- rast(gapfillImages)

par(mfrow=c(2,2))

plot(originalRasters, cex.main=0.75,
     main=sapply(c("Day145, 2019", "Day161, 2019", "Day145, 2020", "Day161, 2020"), 
                 function(s) paste("Original:", s)))
plot(gapfillRasters, cex.main=0.75,
     main=sapply(c("Day145, 2019", "Day161, 2019", "Day145, 2020", "Day161, 2020"), 
                 function(s) paste("Original:", s)))


