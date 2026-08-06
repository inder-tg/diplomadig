############################################################
# Unidad 1 - Sesión 1
# R como herramienta para SIG
# Instalación, conceptos básicos y manejo de archivos
############################################################

# 0. Instalar (una sola vez) y cargar paquetes (en cada sesión nueva)

library(here)
library(tidyverse)
library(sf)
library(terra)
library(mapview)
library(knitr)
library(kableExtra)
library(ggcorrplot)

# 1. Manejo de archivos: .csv 
# Proyecto CalCOFI: https://calcofi.org/
calcofi <- read_csv( here("data", "datos_calcofi.csv") )

# 2. Conceptos básicos: usando CalCOFI como base
# Tipos de objetos
vector_num <- calcofi$mean_T_degC
vector_chr <- calcofi$Dist_cat

df_ej <- data.frame(Temperatura = vector_num, Distancia = vector_chr)
lista_ej <- list(vec = vector_num, df = df_ej)

# Funciones y paquetes
glimpse(calcofi)

mean(vector_num, na.rm = TRUE)   # función base
calcofi %>% summarise(mean_temp = mean(mean_T_degC, na.rm = TRUE)) # tidyverse

# 3. Data wrangling: Parte 1

calcofi_subset <- calcofi %>%
  select(
    - Bottom_D, - Water_density_class, - Cst_Cnt, -IntChl,
    - mean_Phaeop
  )

glimpse(calcofi_subset)

calcofi_subset <- calcofi_subset %>%
  rename(
    Profundidad = Depth_zone,
    ChlorA = mean_ChlorA,
    NO2 = mean_NO2uM,
    NO3 = mean_NO3uM,
    O2 = mean_O2ml_L,
    O2Sat = mean_O2Sat,
    PO4 = mean_PO4uM,
    Salinidad = mean_Salnty,
    SiO3 = mean_SiO3uM,
    Temperatura = mean_T_degC,
    Dist_costa = Distance,
    Lat = Lat_Dec,
    Lon = Lon_Dec,
    Trimestre = Quarter,
    Dist_cost_cat = Dist_cat
  )

glimpse(calcofi_subset)

# Estadísticas globales
calcofi_subset %>%
  summarise(mean_temp = mean(Temperatura, na.rm = TRUE),
            sd_temp   = sd(Temperatura, na.rm = TRUE),
            mean_oxy  = mean(O2, na.rm = TRUE),
            sd_oxy    = sd(O2))

num_vars <- calcofi_subset %>%
  select( where( is.numeric ) )

glimpse(num_vars)

stats <- data.frame(
  mean = sapply(num_vars, function(x) mean(x, na.rm = TRUE) ),
  sd = sapply(num_vars, function(x) sd(x, na.rm = TRUE) )
)

print(stats)

calcofi_stats <- calcofi_subset %>%
  select( -c(Dist_costa, Lat, Lon, Year, Trimestre)) %>%
  summarise( across( where(is.numeric),
                     list(
                       mean = ~mean(.x, na.rm=TRUE),
                       sd = ~sd(.x, na.rm=TRUE)
                     )) )

kable(calcofi_stats,
      caption = "Estadíscticas básicas de las variables numéricas del (CalCOFI)")

kable(calcofi_stats, caption = "Estadísticas descriptivas de variables numéricas (CalCOFI)") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"),
                full_width = FALSE, position = "center")

names(calcofi_stats) <- gsub(".*_", "", names(calcofi_stats))


kable(calcofi_stats, digits = 2, caption = "Estadísticas descriptivas de variables físico-químicas (CalCOFI)") %>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>%
  add_header_above(c(
                     "ChlorA" = 2,
                     "NO2" = 2,
                     "NO3" = 2,
                     "O2" = 2,
                     "O2Sat" = 2,
                     "PO4" = 2,
                     "Salinidad"   = 2,
                     "SiO3" = 2,
                     "Temperatura" = 2))   # cada variable agrupa sus columnas
# --- Hacer ejercicio de esto de arriba pero que usen median y mad

# 4. Gráficos Parte 1
# Boxplot
boxplot( calcofi_subset$Temperatura ~ calcofi_subset$Trimestre,
         main = "Temperatura por trimestre",
         xlab = "Trimestre",
         ylab = "Temperatura (°C)")

boxplot( calcofi_subset$Temperatura ~ calcofi_subset$Profundidad,
         main = "Temperatura por profundidad",
         xlab = "Trimestre",
         ylab = "Temperatura (°C)")

boxplot( calcofi_subset$Temperatura ~ calcofi_subset$Dist_cost_cat,
         main = "Temperatura por distancia a la costa",
         xlab = "Trimestre",
         ylab = "Temperatura (°C)")

# 5. Data wrangling Parte 2
# Seleccionar variable y filtrar por Q2, very shallow & Costera

calcofi_Q2_veryShallow_costera <- calcofi_subset %>%
  filter(Trimestre == 2, Profundidad == "Very shallow", 
         Dist_cost_cat == "Costera") %>%
  select(Year, Temperatura)

glimpse(calcofi_Q2_veryShallow_costera)

# Serie de tiempo anual de temperatura
temp_year <- calcofi_Q2_veryShallow_costera %>%
  group_by(Year) %>%
  summarise(mean_temp = mean(Temperatura, na.rm = TRUE))

glimpse(temp_year)

# 6. Gráficos Parte 2
# Serie de tiempo con ggplot
ggplot(temp_year, aes(x = Year, y = mean_temp)) +
  geom_line(color = "red", linewidth = 1) +
  geom_point(color = "blue") +
  labs(title = "Temperatura promedio anual (CalCOFI)",
       x = "Año", y = "Temperatura (°C)") +
  theme_minimal()

# Correlation matrix con GGally
calcofi_Q2_veryShallow_costera <- calcofi_subset %>%
  filter(Trimestre == 2, Profundidad == "Very shallow", 
         Dist_cost_cat == "Costera") %>%
  select(ChlorA, NO2, NO3, O2, O2Sat, PO4, Salinidad, SiO3,
         Temperatura)

glimpse(calcofi_Q2_veryShallow_costera)

corMatrix <- cor(calcofi_Q2_veryShallow_costera)

ggcorrplot(corMatrix, hc.order = TRUE, type = "lower",
           lab = TRUE, lab_size = 3,
           colors = c("red", "white", "blue"),
           title = "Matriz de correlación entre variables ambientales")

# Ejercicio: filtrar dataset por primer trimestre, Deep, Oceánica y graficar
# matriz de correlación de todas las variables numéricas
# y graficar serie de tiempo de temperatura promedio anual

# 7. Manejo de archivos: .shp y visualizaciones espaciales

# Shapefile de línea de costa (ejemplo: Natural Earth)
# Natural Earth: https://www.naturalearthdata.com/downloads/10m-physical-vectors/10m-coastline/
  
coast <- st_read( here( "data", "ne_10m_coastline", "ne_10m_coastline.shp") )

# Recortar al área de CalCOFI
bbox_calcofi <- st_bbox(calcofi_subset %>% 
                          select(Lon, Lat) %>% 
                          drop_na() %>% 
                          st_as_sf(coords = c("Lon","Lat"), crs = st_crs(4326) ))

coast_crop <- st_crop(coast, bbox_calcofi)


calcofi_sf <- st_as_sf(calcofi_subset,
                       coords = c("Lon", "Lat"),
                       crs = st_crs(4326))

mapview(coast_crop, color = "black") +
  mapview(calcofi_sf, col.regions = "red")

# ----

prueba <- calcofi_subset %>% 
  filter( Trimestre==2, 
          Profundidad == "Very shallow",
          Dist_cost_cat == "Costera")

calcofi_sf <- st_as_sf(prueba,
                       coords = c("Lon", "Lat"),
                       crs = 4326)


mapview(coast_crop, color = "black") +
  mapview(calcofi_sf, col.regions = "red")

