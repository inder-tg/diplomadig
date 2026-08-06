# -----------------------------------------------------------------------------
#
# Elaborado por Inder Tecuapetla, May 31, 2023
#
# Modificado Julio 28, 2023
#
# Instalación de paquetes/bibliotecas a utilizar en este módulo
# 
# Hecho para SELPER/CEOS Working Group Chapter D Training Group  
#
# Modificado Agosto 6, 2026
# 
# Instalación de paquetes a utilizar en Diplomado en Geomática, UNAM, 2026

# -----------------------------------------------------------------------------

# neededPackages <- c("raster", 
#                     # "mapview", 
#                     "sp", "RColorBrewer", "gtools",
#                     "foreach", "doParallel", "Kendall", "trend",
#                     "geoTS", "sta", "rasterVis", "bfast", "tmap", 
#                     "ggplot2", "dplyr", "terra", "imputeTS",
#                     "numbers", "itertools", "gapfill", "heatmaply",
#                     "vcd", "sephora")
# 
# packagesToInstall <- neededPackages[!(neededPackages %in% installed.packages()[,"Package"])]
# 
# if( length(packagesToInstall) ){
#   for( i in 1:length(packagesToInstall) ){
#     message("Installing package", packagesToInstall[i], "\n")
#     install.packages(packagesToInstall[i], dependencies = TRUE)
#   }
# } 


# =============================================================================
# --- Agosto 2026
# =============================================================================

neededPackages <- c("here", "readr", "sf", "tidyverse",
                    "leaflet", "dplyr", "ggplot2",
                    "ggcorrplot", "kableExtra",
                    "tmap", "terra", "mapview", "knitr", "tinytex",
                    "RColorBrewer", "gtools", "foreach", "doParallel",
                    "bfast", "geoTS"
                    )

packagesToInstall <- setdiff(neededPackages, rownames(installed.packages()))

if( length(packagesToInstall) > 0 ){
  for( i in 1:length(packagesToInstall) ){
    message("Installing package: ", packagesToInstall[i], "\n")
    install.packages(packagesToInstall[i], dependencies = TRUE)
  }
} 

# --- Requerido para compilar RMD -> PDF
library(tinytex)
tinytex::install_tinytex()

# --- Requerido para algunos plots de sephora
(!require("BiocManager", quietly = TRUE))
install.packages("BiocManager", source = TRUE)

library(BiocManager)
BiocManager::install("ComplexHeatmap")


