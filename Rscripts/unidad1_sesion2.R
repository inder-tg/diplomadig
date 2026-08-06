############################################################
# Unidad 1 - Sesión 2
# R como herramienta SIG
# Agosto 1, 2026
############################################################

############################################################
# 1. Carga de paquetes 
############################################################

# --- Paquetes necesarios
library(leaflet)
library(tmap)
library(mapview)
library(sf)
library(readr)
library(ggplot2)
library(dplyr)
library(here)

############################################################
# 2. Manejo de archivos (.csv y .shp)
############################################################

# Cargar CSV de cines
csvFILES <- list.files(path = here("data", "cinesCDMX"),
                       pattern = ".csv",
                       full.names = TRUE)

cines <- read_csv(file = csvFILES,
                  locale = locale(encoding = "latin1"))

# Cargar shapefile de colonias
shpFILES <- list.files(path = here("data", "cinesCDMX", "coloniascdmx"),
                       pattern = ".shp",
                       full.names = TRUE)

SHP <- st_read( shpFILES )

############################################################
# 3. Data wrangling (base vs tidyverse)
############################################################

# Filtrar solo CDMX
cinesCDMX <- cines %>%
  filter(`Entidad federativa` == "CIUDAD DE MÉXICO")

# Conteo en R base
table(cinesCDMX$Municipio)

# Conteo de cines por municipio (tidyverse)
conteo <- cinesCDMX %>%
  group_by(Municipio) %>%
  summarise(num_cines = n()) %>%
  arrange(desc(num_cines))
# arrange(Municipio)

conteo
############################################################
# 4. Visualización
############################################################

### 4.1 R base
par(mar = c(5, 12, 4, 2))
barplot(conteo$num_cines,
        names.arg = conteo$Municipio,
        col = "steelblue",
        main = "Distribución de cines por alcaldía en CDMX",
        horiz = TRUE, las = 1, cex.names = 0.7)

### 4.2 ggplot
conteo %>%
  ggplot(aes(x = reorder(Municipio, num_cines), y = num_cines)) +
  geom_col(fill = "steelblue") +
  labs(title = "Distribución de cines por alcaldía en CDMX",
       x = "Municipio",
       y = "Número de cines") +
  coord_flip() +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

ggplot(conteo, aes(x = Municipio, y = num_cines)) +
  geom_col(fill = "steelblue") +
  labs(title = "Distribución de cines por alcaldía en CDMX",
       x = "Municipio", y = "Número de cines") +
  coord_flip()

conteo_ord_alfabetico <- conteo %>%
  arrange(Municipio) %>%   # ordena alfabéticamente
  mutate(Municipio = factor(Municipio, levels = rev(Municipio)))

ggplot(conteo_ord_alfabetico, aes(x = Municipio, y = num_cines)) +
  geom_col(fill = "steelblue") +
  geom_text(aes(label = num_cines), 
            hjust = -0.1, size = 3) +   # etiquetas al lado de cada barra
  labs(title = "Distribución de cines por alcaldía en CDMX",
       x = "Municipio", y = "Número de cines") +
  coord_flip() +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

### 4.1 Leaflet
leaflet(data = cinesCDMX) %>%
  addTiles() %>%
  addCircleMarkers(~Longitud, ~Latitud,
                   popup = ~`Nombre de la Unidad Económica`,
                   radius = 4, color = "red")

### 4.4 mapview

# Convertir cinesCDMX a sf
cinesCDMX_sf <- st_as_sf(cinesCDMX,
                         coords = c("Longitud", "Latitud"),
                         crs = 4326)

mapview(SHP) + mapview(cinesCDMX_sf)

### 4.5 tmap
tmap_mode("view")
tm_shape(SHP) +
  tm_polygons() +
  tm_shape(cinesCDMX %>% st_as_sf(coords = c("Longitud","Latitud"), crs = 4326)) +
  tm_dots(col = "red")


# 5. Data wrangling

# --- agregamos info de num cines x alcaldia a cinesCDMX

cinesCDMX <- cinesCDMX %>%
  add_count( Municipio, name = "n_cines" )

leaflet( data = cinesCDMX ) %>%
  addTiles() %>%
  addCircleMarkers( ~Longitud, ~Latitud,
                    radius = 4, 
                    color = "red", 
                    popup = ~paste0(
                      `Nombre de la Unidad Económica`, "<br>",
                      "<b>Alcaldía: </b>", Municipio, "<br>",
                      "<b>No. cines por Alcaldía: </b>", n_cines
                    ))

# --- Cines en la Cuauhtemoc

SHP_alcaldia_cuauhtemoc <- SHP %>%
  filter( alcaldi == "CUAUHTEMOC" )

plot( st_geometry(SHP_alcaldia_cuauhtemoc) )

cinesCuauhtemoc <- cinesCDMX %>%
  filter( Municipio == "Cuauhtémoc" )

leaflet( data = SHP_alcaldia_cuauhtemoc ) %>%
  addTiles() %>%
  addPolygons( fillColor = "lightblue",
               popup = ~paste("Alcaldía: ", alcaldi)) %>%
  addCircleMarkers( data = cinesCuauhtemoc,
                    ~Longitud, ~Latitud,
                    radius = 4, 
                    color = "red", 
                    
                    popup = ~paste0(
                      `Nombre de la Unidad Económica`, "<br>",
                      "<b>Alcaldía: </b>", Municipio, "<br>",
                      "<b>No. cines por Alcaldía: </b>", n_cines
                    ))

# --- cines por colonia en la Cuauhtémoc

cinesCuauhtemoc_sf <- st_as_sf(cinesCuauhtemoc,
                               coords = c("Longitud", "Latitud"),
                               crs = 4326)

plot(st_geometry(cinesCuauhtemoc_sf))

cinesCuauhtemoc_colonia <- st_join(cinesCuauhtemoc_sf, SHP_alcaldia_cuauhtemoc)

conteo_colonias <- cinesCuauhtemoc_colonia %>%
  group_by( nombre ) %>%
  summarise( num_cines_colonia = n() ) %>%
  st_drop_geometry()


SHP_alcaldia_cuauhtemoc <- SHP_alcaldia_cuauhtemoc %>%
  left_join(conteo_colonias, by = "nombre" )


leaflet( data = SHP_alcaldia_cuauhtemoc ) %>%
  addTiles() %>%
  addPolygons( fillColor = "lightblue",
               popup = ~paste("Alcaldía: ", alcaldi, "<br>",
                              "Colonia: ", nombre, "<br>",
                              "Número cines en la colonia: ", num_cines_colonia)) %>%
  addCircleMarkers( data = cinesCuauhtemoc_colonia,
                    # ~Longitud, ~Latitud,
                    radius = 4, 
                    color = "red", 
                    popup = ~paste0(
                      `Nombre de la Unidad Económica`, "<br>",
                      "<b>Alcaldía: </b>", Municipio, "<br>",
                      "<b>No. cines por Alcaldía: </b>", n_cines
                    ))

# --- tmap

tm_shape(SHP_alcaldia_cuauhtemoc) +
  tm_polygons("num_cines_colonia",
              palette = "Reds",
              title = "Número de cines",
              border.col = "gray") +
  tm_shape(cinesCuauhtemoc_colonia) +
  tm_dots(col = "blue", size = 0.8,
          popup.vars = c("Nombre de la Unidad Económica" = "Nombre de la Unidad Económica",
                         "Colonia" = "nombre"))



tm_shape(SHP_alcaldia_cuauhtemoc) +
  tm_polygons("num_cines_colonia",
              palette = "brewer.reds",
              fill.legend = "Número de cines",
              border.col = "gray",
              colorNA = "transparent",
              textNA = "Sin cines",
              popup.vars = c( "Colonia" = "nombre", "Alcaldía" = "alcaldi" )) +
  tm_shape(cinesCuauhtemoc_colonia) +
  tm_dots(fill = "blue", size = 0.8,
          popup.vars = c("Nombre de la Unidad Económica" = "Nombre de la Unidad Económica",
                         "Colonia" = "nombre"))
