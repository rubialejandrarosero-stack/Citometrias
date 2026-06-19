# 01_leer_fcs.R — ejemplo mínimo de lectura y visualización de un archivo FCS
# Entorno: conda activate citometrias   (paquetes: flowCore, ggcyto, tidyverse)

suppressPackageStartupMessages({
  library(flowCore)
  library(ggcyto)
})

# Busca el primer .fcs en data/raw/
archivos <- list.files("data/raw", pattern = "\\.fcs$", full.names = TRUE, ignore.case = TRUE)

if (length(archivos) == 0) {
  message("No hay archivos .fcs en data/raw/. Copia alguno ahí y vuelve a ejecutar.")
  message("Mientras tanto, verificamos que flowCore y ggcyto cargan correctamente:")
  message("  flowCore  ", as.character(packageVersion("flowCore")))
  message("  ggcyto    ", as.character(packageVersion("ggcyto")))
  quit(save = "no", status = 0)
}

fcs <- archivos[1]
message("Leyendo: ", fcs)
ff <- read.FCS(fcs, transformation = FALSE)

# Resumen rápido
print(ff)
cat("\nMarcadores/canales:\n")
print(pData(parameters(ff))[, c("name", "desc")])

# Gráfico de ejemplo: primeros dos canales (ajusta los nombres a tus marcadores)
canales <- colnames(ff)
if (length(canales) >= 2) {
  dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
  p <- ggcyto(ff, aes(x = !!canales[1], y = !!canales[2])) +
    geom_hex(bins = 128) +
    theme_bw()
  ggsave("output/figures/ejemplo_scatter.png", plot = as.ggplot(p),
         width = 6, height = 5, dpi = 150)
  message("Gráfico guardado en output/figures/ejemplo_scatter.png")
}
