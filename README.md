# diplomadig

Este repositorio contiene materiales presentados en diversos módulos del Diplomado
en Geomática organizado por el Instituto de Geografía de la UNAM

<details>

<summary>XVI Edición, 2023</summary>

Módulo X "Análisis de series de tiempo de imágenes satelitales con R"

## 2023
- Creación de este repositorio
- Directorio /rspatial usado en ```intro_RSIG.R``` se descarga desde el folder asignado en la nube
- Directorio /LaPiedad usado en ```prelim_trendAnalysis_maiz_agave.R``` se descarga desde el folder asignado en la nube
- Directorio /muni_2018 usado en ```trendAnalysis_maiz_agave.R``` se descarga desde el folder asignado en la nube
</details>

<details>
<summary>XVII Edición, 2024</summary>

Módulo IX "Percepción Remota: Análisis de series de tiempo de imágenes satelitales con R"

## Primera Sesión (Agosto 10): 

- Directorio /data/rspatial usado en ```intro_RSIG.R```

- Directorio /data/ANP distribuido a través de la nube del Diplomado y usado en ```mohinora_tmap.R```

- Directorio /data/USV7 distribuido a través de la nube del Diplomado y usado en ```mohinora_tmap.R```

## Segunda Sesión (Agosto 16)

- Directorio /data/mohinora usado en ```mohinora_imputation.R```

- Uso de ```mohinora_imputation.R```, ```mohinora_anomalies.R``` y ```mohinora_trendAnalysis.R``` 

## Tercera sesión (Agosto 17)

- Uso de ```mohinora_trendAnalysis.R``` y ```mohinora_cps.R```

## Cuarta sesión (Agosto 23)

- Uso de ```mohinora_cps.R``` y ```mohinora_sephora.R```.

**NOTA:** El portal de residencia de ```ComplexHeatmap``` es [Bioconductor](https://www.bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html). Por esta razón,
las instrucciones para instalar el paquete ```ComplexHeatmap``` son ligeramente distintas a las discutidas
hasta ahora. De acuerdo al portal mencionado, las instrucciones para installar ```ComplexHeatmap``` en R (en versiones superiores a la 4.4) son:

```{r, eval=FALSE}
(!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")
```

</details>


<details>
<summary>XVIII Edición, 2025</summary>

# Bloque 2, Módulo V "R como herramienta de SIG"

## Primera Sesión (Abril 4): 

- Directorio /data/rspatial usado en ```intro_RSIG.R```

- Directorio /data/ANP distribuido a través de la nube del Diplomado y usado en ```mohinora_tmap.R```

- Se empleó el script ```intro_RSIG.R``` y materiales auxiliares

## Segunda Sesión (Abril 5):

- Creación del directorio /data/mohinora/250m_16_days_NDVI_QA

- Creación del directorio /data/outputs

- Creación del directorio /data/outputs/mohinora_QA

- Creación del directorio /data/outputs/mohinora_imputation

- Se emplearon los scripts ```mohinora_QA.R``` y ```mohinora_imputation.R``` junto con algunos materiales auxiliares

## Tercera Sesión (Abril 11):

- Creación del directorio /data/mohinora/250m_16_days_NDVI_QA_byYears y los sub-directorios auxiliares

- Creación del directorio /data/outputs/mohinora_interpolation

- Creación del directorio /data/outputs/mohinora_igapfill

- Creación del directorio /data/outputs/mohinora_anomalies

- Se emplearon los scripts ```mohinora_interpolation.R```, ```mohinora_igapfill.R``` y ```mohinora_anomalies.R``` junto con algunos materiales auxiliares

## Cuarta Sesión (Abril 12):

- Creación del directorio /data/outputs/mohinora_trends

- Creación del directorio /data/outputs/mohinora_cps

- Se emplearon los scripts ```mohinora_trends.R```, ```mohinora_cps.R``` y ```mohinora_sephora.R``` junto con algunos materiales auxiliares
</details>



