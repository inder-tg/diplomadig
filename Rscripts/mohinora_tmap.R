
# --- Diplomado Geomática, IG, UNAM, 2024
# --- Módulo IX: Percepcion Remota: Análisis de series de tiempo de imágenes satelitales con R
# --- Diplomado Geomática, IG, UNAM, 2025
# --- Bloque 2: Sistemas de Información Geográfica
# --- Módulo V: R como herramienta de SIG

# --- Elaborado: Agosto 8, 2024
# --- Actualizado: Abril 4, 2025

# --- En este script usamos el paquete tmap para generar mapas estáticos de Cerro
# --- Mohinora, Chihuahua y de los tipos de uso de suelo y vegetación de esta
# --- región de interés.

# --- DATASET: 
# --- El directorio ANP debe descargarse de la nube del Diplomado

# --- Preámbulo END

library(sf)
library(tmap)

# ---

DIR <- paste0(getwd(), "/data")

listDIRS <- list.dirs(path=DIR)

anpFILES <- list.files(path = listDIRS[2],
                   full.names = TRUE,
                   pattern = ".shp$")

usvFILES <- list.files(path = listDIRS[6],
                  full.names = TRUE,
                  pattern = ".shp$")

ndviFILES <- list.files(path = listDIRS[4],
                        full.names = TRUE,
                        pattern = ".tif$")

shpANP <- read_sf(anpFILES)

mohinora_USV <- read_sf(usvFILES)

# ---

which(shpANP$NOMBRE == "Cerro Mohinora")

plot(shpANP[144,])

plot( st_geometry(shpANP[144,]), lwd = 4 )

str(mohinora_USV)

mohinora_USV$geometry

shpANP$geometry

mohinora_shp_GCS <- shpANP[144,]

mohinora_shp_sinu <- st_transform(x=mohinora_shp_GCS, 
                                  crs=st_crs(mohinora_USV))

# ---

usv_COLORS <- c("#A1E5A5", "#E9D66B", "#00A877", 
                "#66B032", "#83A4F0", "#FC8FAB", "#F500A1")

usv_NAMES <- c("Pino", "Pastizal", "Pino-Encino",
               "Ayarin", "Agro", "Arbustiva", "Arborea")

COLOR_USV <- c(usv_COLORS[1], 
               rep(usv_COLORS[3],2), 
               usv_COLORS[4], 
               rep(usv_COLORS[2],2),
               rep(usv_COLORS[5], length(7:12)), 
               rep(usv_COLORS[6], 9),
               usv_COLORS[7])

# --- Adding variable COLOR to sf-object mohinora_USV 
mohinora_USV$COLOR <- COLOR_USV

shp_mohinora_rect <- tm_shape(mohinora_USV) +
  tm_borders(lwd = 3) +
  tm_graticules(n.x=4,
                labels.size=1.5) +
  tm_compass(type = "rose", size=4,
             position = c("right", "top")) +
  tm_scale_bar(text.size = 0.75,
    position = c("right", "bottom")) +
  tm_fill(col= "COLOR") +
  tm_add_legend("symbol", 
                labels=usv_NAMES, 
                col=usv_COLORS,
                border.col = "grey40",
                size=4,
                shape=15,
                is.portrait = FALSE) +
  tm_layout(legend.outside = TRUE,
            legend.outside.position = "bottom",
            legend.outside.size = 0.115,
            legend.text.size = 4,
            legend.text.fontface = 2)

shp_mohinora_rect

# ---

bbox_new <- st_bbox(mohinora_USV) # current bounding box

xrange <- bbox_new$xmax - bbox_new$xmin # range of x values
yrange <- bbox_new$ymax - bbox_new$ymin # range of y values

bbox_new[1] <- bbox_new[1] - (0.25 * xrange) # xmin - left
bbox_new[3] <- bbox_new[3] + (0.25 * xrange) # xmax - right
bbox_new[2] <- bbox_new[2] - (0.3 * yrange) # ymin - bottom
bbox_new[4] <- bbox_new[4] + (0.25 * yrange) # ymax - top

bbox_new <- bbox_new %>%  # take the bounding box ...
  st_as_sfc() # ... and make it a sf polygon

# ---

shp_mohinora_fill <- tm_shape(mohinora_USV, bbox=bbox_new) +
  tm_borders(lwd = 3) +
  tm_graticules(n.x=4,
                labels.size=1.5) +
  tm_compass(type = "rose", size=4,
             position = c("right", "top")) +
  tm_scale_bar(text.size = 0.75,
    position = c("right", "bottom")) +
  tm_fill(col= "COLOR") +
  tm_add_legend("symbol", 
                labels=usv_NAMES, 
                col=usv_COLORS,
                border.col = "grey40",
                size=4,
                shape=15,
                is.portrait = FALSE) +
  tm_layout(legend.outside = TRUE,
            legend.outside.position = "bottom",
            legend.outside.size = 0.115,
            legend.text.size = 4,
            legend.text.fontface = 2)#,

shp_mohinora_fill  

# ---

mohinora_USV_lines <- mohinora_USV %>%
  st_cast("MULTILINESTRING")


shp_mohinora_lines <- tm_shape(mohinora_USV_lines, 
                      bbox = bbox_new) +
  tm_lines(col = "COLOR", lwd=3) +
  tm_compass(type = "8star", position = c("right", "bottom")) +
  tm_scale_bar(text.size = 0.65, position = c("right", "bottom")) +
  tm_add_legend("line",
                labels=usv_NAMES,
                col=usv_COLORS,
                border.col = "grey40",
                lwd=3) +
  tm_layout(frame = FALSE, bg.color = NA,
            legend.outside = FALSE,
            legend.stack = "horizontal",
            legend.position = c("left", "bottom"))

shp_mohinora_lines

# ---

st_write(mohinora_shp_sinu,
         dsn = paste0( getwd(), "/data/outputs/mohinora_shp_sinu.shp" ))

