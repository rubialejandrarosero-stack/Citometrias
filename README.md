# Citometrias

Análisis de **citometría de flujo** en R (ecosistema Bioconductor: `flowCore`, `ggcyto`,
`flowWorkspace`, `openCyto`).

## Estructura
```
Citometrias/
├── data/raw/        # archivos .fcs de entrada (NO se versionan)
├── R/               # scripts de análisis
├── output/figures/  # gráficos y resultados generados
├── environment.yml  # entorno conda reproducible del proyecto
└── README.md
```

## Puesta en marcha
```bash
# 1) Crear el entorno del proyecto (una vez)
conda env create -f environment.yml

# 2) Activarlo
conda activate citometrias

# 3) Ejecutar el script de ejemplo
Rscript R/01_leer_fcs.R
```

En VS Code, abre la carpeta y acepta las extensiones recomendadas (R, etc.). El proyecto ya apunta
al intérprete de R del entorno `citometrias` (ver `.vscode/settings.json`).

## Datos
Coloca tus `.fcs` en `data/raw/`. Están excluidos de git por defecto (suelen ser pesados); ajusta
`.gitignore` si quieres versionar alguno pequeño de ejemplo.
