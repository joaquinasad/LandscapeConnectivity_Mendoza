library(dplyr)
library(sf)

rm(list = ls())

setwd("C:/Users/ASUS/OneDrive - Facultad de Agronomía - Universidad de Buenos Aires/Escritorio/Objs_3y4")
datos <- read.csv("Nodos_IPC_AreaDeEstudio_1985_2024_v2.csv")

#Solo año 2024
nodos_2024 <- datos %>% filter(year == 2024)

#Preparación espacial y cálculo de matriz de distancias
#Convertir el dataframe a un objeto espacial (puntos) usando las coordenadas
nodos_sf <- st_as_sf(nodos_2024, coords = c("x_lon", "y_lat"), crs = 4326)

#Proyectar a un CRS métrico local para Mendoza (EPSG:5344) para calcular distancias en metros
nodos_sf_proj <- st_transform(nodos_sf, crs = 5344)

#Calcular matriz de distancias euclidianas entre todos los centroides
dist_matrix <- st_distance(nodos_sf_proj)
dist_matrix <- as.numeric(dist_matrix)
dist_matrix <- matrix(dist_matrix, nrow = nrow(nodos_2024), ncol = nrow(nodos_2024))

#Parametrización de la dispersión ecológica (Probabilidad p_ij)
#La distancia a la cual la probabilidad de dispersión es del 50% (0.5)
distancia_media_dispersion <- 5000

#Constante de decaimiento exponencial
k <- -log(0.5) / distancia_media_dispersion

#Matriz de probabilidades de dispersión (de 0 a 1)
p_ij <- exp(-k * dist_matrix)
diag(p_ij) <- 1 # La probabilidad de un parche consigo mismo es 1

#Cálculo de la Importancia del Nodo (Aproximación de Flujo dPC)
areas <- nodos_2024$area_ha
n_nodos <- length(areas)
dPC_flujo <- numeric(n_nodos)

#Iterar sobre cada nodo para sumar su contribución ponderada a la red
for (i in 1:n_nodos) {
  flujo_i <- 0
  for (j in 1:n_nodos) {
    # area del nodo de origen * area del destino * probabilidad de conexión
    flujo_i <- flujo_i + (areas[i] * areas[j] * p_ij[i, j])
  }
  dPC_flujo[i] <- flujo_i
}

#Transformación logarítmica y normalización (0 a 1)
#Aplicamos log10 para atenuar el sesgo del parche masivo (Piedemonte)
#Se suma 1 para evitar problemas matemáticos con log(0) si existieran valores muy bajos
dPC_log <- log10(dPC_flujo + 1)

#Normalización estricta de 0 a 1
min_val <- min(dPC_log)
max_val <- max(dPC_log)

nodos_2024$dPC_Bruto <- dPC_flujo
nodos_2024$Importancia_Relativa <- (dPC_log - min_val) / (max_val - min_val)

#Limpiar columnas y guardar el CSV final para subir a GEE
nodos_exportar <- nodos_2024 %>% 
  select(id_nodo, x_lon, y_lat, area_ha, Importancia_Relativa)

View(nodos_exportar)

#Histograma
hist(nodos_exportar$Importancia_Relativa)

#Guardar
write.csv(nodos_exportar, "Nodos_Importancia_2024_v2.csv", row.names = FALSE)
