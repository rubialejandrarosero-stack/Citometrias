# R/10_infiltracion_gating.R
# Gating jerarquico del panel de infiltracion (manuscrito) sobre data/raw/Infiltracion.
# Estrategia (definida por el usuario):
#   CD45+ -> CD3+  -> CD4, CD8, CD4/HLA-DR+, CD8/HLA-DR+
#            CD3- CD14+ -> CD64+/HLA-DR+, CD64+/CD11b+
#            CD3- CD14- -> CD19+ (B), CD16+ (NK)
# Umbrales objetivos = p99.9 del control Unstained (validados con single-stains).
# NOTA: gating automatico por umbrales; puede diferir de compuertas manuales (FlowJo).

suppressPackageStartupMessages({library(flowCore); library(dplyr); library(readr)})

DIR <- "data/raw/Infiltracion"
REF <- file.path(DIR, "Reference Group")
out_dir <- "output/Infiltracion"; dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ch <- c(CD45="APC-Cy7-A", CD3="PE-Cy7-A", CD4="BV421-A", CD8="BV510-A", CD14="PE-A",
        CD16="FITC-A", CD19="StarBrightViolet 475-A", CD64="PerCP-Cy5.5-A",
        CD11b="APC-A", HLADR="Pacific Blue-A")
FSC_MIN <- 5e5
rd <- function(f) read.FCS(f, transformation = FALSE, truncate_max_range = FALSE)

# --- Calibracion: Logicle + umbrales p99.9 del unstained ---
lg <- estimateLogicle(rd(file.path(DIR, "48h/INF_ACT/donante_3.fcs")), channels = unname(ch))
tdf <- function(f) { e <- as.data.frame(exprs(transform(rd(f), lg))); e[e[["FSC-A"]] > FSC_MIN, , drop = FALSE] }
un <- tdf(file.path(REF, "Unstained (Cells).fcs"))
th <- setNames(sapply(ch, function(c) as.numeric(quantile(un[[c]], 0.999, na.rm = TRUE))), names(ch))
message("Umbrales calibrados: ", paste(sprintf("%s=%.2f", names(th), th), collapse = "  "))

pc <- function(n, d) if (d > 0) round(n / d * 100, 3) else NA_real_

gate <- function(fpath) {
  e <- tdf(fpath)
  P <- function(m) e[[ch[[m]]]]
  n_cells <- nrow(e)
  cd45 <- e[P("CD45") > th["CD45"], , drop = FALSE]
  g <- function(d, m, op = `>`) d[op(d[[ch[[m]]]], th[m]), , drop = FALSE]
  # T cells
  cd3  <- g(cd45, "CD3")
  cd3n <- cd45[cd45[[ch[["CD3"]]]] <= th["CD3"], , drop = FALSE]
  cd4  <- g(cd3, "CD4");  cd8 <- g(cd3, "CD8")
  cd4h <- g(cd4, "HLADR"); cd8h <- g(cd8, "HLADR")
  # Myeloid
  cd14 <- g(cd3n, "CD14")
  cd64h <- cd14[cd14[[ch[["CD64"]]]] > th["CD64"] & cd14[[ch[["HLADR"]]]] > th["HLADR"], , drop = FALSE]
  cd64c <- cd14[cd14[[ch[["CD64"]]]] > th["CD64"] & cd14[[ch[["CD11b"]]]] > th["CD11b"], , drop = FALSE]
  # CD3- CD14-
  dn <- cd3n[cd3n[[ch[["CD14"]]]] <= th["CD14"], , drop = FALSE]
  b  <- g(dn, "CD19"); nk <- g(dn, "CD16")

  n <- list(Cells=n_cells, CD45=nrow(cd45), CD3=nrow(cd3), CD4=nrow(cd4), CD8=nrow(cd8),
            CD4_HLADR=nrow(cd4h), CD8_HLADR=nrow(cd8h), Mono_CD14=nrow(cd14),
            CD64_HLADR=nrow(cd64h), CD64_CD11b=nrow(cd64c),
            CD3neg_CD14neg=nrow(dn), B_CD19=nrow(b), NK_CD16=nrow(nk))
  data.frame(
    n_Cells=n$Cells, n_CD45=n$CD45,
    n_CD3=n$CD3, n_CD4=n$CD4, n_CD8=n$CD8, n_CD4_HLADR=n$CD4_HLADR, n_CD8_HLADR=n$CD8_HLADR,
    n_Mono_CD14=n$Mono_CD14, n_CD64_HLADR=n$CD64_HLADR, n_CD64_CD11b=n$CD64_CD11b,
    n_CD3neg_CD14neg=n$CD3neg_CD14neg, n_B_CD19=n$B_CD19, n_NK_CD16=n$NK_CD16,
    # % de CD45
    pctCD45_CD3=pc(n$CD3,n$CD45), pctCD45_CD4=pc(n$CD4,n$CD45), pctCD45_CD8=pc(n$CD8,n$CD45),
    pctCD45_CD4HLADR=pc(n$CD4_HLADR,n$CD45), pctCD45_CD8HLADR=pc(n$CD8_HLADR,n$CD45),
    pctCD45_Mono=pc(n$Mono_CD14,n$CD45), pctCD45_CD64HLADR=pc(n$CD64_HLADR,n$CD45),
    pctCD45_CD64CD11b=pc(n$CD64_CD11b,n$CD45), pctCD45_B=pc(n$B_CD19,n$CD45), pctCD45_NK=pc(n$NK_CD16,n$CD45),
    # % del parental
    pctPar_CD3_de_CD45=pc(n$CD3,n$CD45), pctPar_CD4_de_CD3=pc(n$CD4,n$CD3), pctPar_CD8_de_CD3=pc(n$CD8,n$CD3),
    pctPar_CD4HLADR_de_CD4=pc(n$CD4_HLADR,n$CD4), pctPar_CD8HLADR_de_CD8=pc(n$CD8_HLADR,n$CD8),
    pctPar_Mono_de_CD45=pc(n$Mono_CD14,n$CD45), pctPar_CD64HLADR_de_Mono=pc(n$CD64_HLADR,n$Mono_CD14),
    pctPar_CD64CD11b_de_Mono=pc(n$CD64_CD11b,n$Mono_CD14),
    pctPar_B_de_DN=pc(n$B_CD19,n$CD3neg_CD14neg), pctPar_NK_de_DN=pc(n$NK_CD16,n$CD3neg_CD14neg)
  )
}

# --- Procesar todas las muestras del manifest ---
man <- read_csv(file.path(DIR, "manifest.csv"), show_col_types = FALSE) %>% filter(tipo == "muestra")
message("Procesando ", nrow(man), " muestras...")
res <- bind_rows(lapply(seq_len(nrow(man)), function(i) {
  r <- man[i, ]
  cbind(data.frame(Tiempo = r$tiempo, Condicion = r$condicion, Activacion = r$activacion,
                   Donante = r$donante, Archivo = r$archivo_original),
        gate(file.path(DIR, r$archivo_destino)))
})) %>% arrange(Tiempo, Condicion, Activacion, as.integer(Donante))

# --- QC: exclusión de CD4/CD8 por spreading violeta (BV421/BV510) ---
# La muestra NO_INF/NO_ACT/96h/donante 1 colapsa en diagonal CD4/CD8 (CD4%+CD8% de CD3 = 188%,
# casi todo doble-positivo): su reparto CD4/CD8 no es fiable. Se anulan SOLO las columnas
# específicas de CD4/CD8 (CD3, monocitos, B, NK quedan intactos). Ver auditoría de spreading.
cd48_cols <- c("n_CD4","n_CD8","n_CD4_HLADR","n_CD8_HLADR",
               "pctCD45_CD4","pctCD45_CD8","pctCD45_CD4HLADR","pctCD45_CD8HLADR",
               "pctPar_CD4_de_CD3","pctPar_CD8_de_CD3","pctPar_CD4HLADR_de_CD4","pctPar_CD8HLADR_de_CD8")
qc_bad <- with(res, Tiempo == 96 & Condicion == "NO_INF" & Activacion == "NO_ACT" & as.integer(Donante) == 1)
res[qc_bad, cd48_cols] <- NA
message("QC: CD4/CD8 anulado en ", sum(qc_bad), " muestra(s) por spreading violeta severo.")

write_csv(res, file.path(out_dir, "infiltracion_conteos_porcentajes.csv"))
message("Guardado: ", file.path(out_dir, "infiltracion_conteos_porcentajes.csv"))
print(res %>% select(Tiempo, Condicion, Activacion, Donante, n_CD45, n_CD3, pctCD45_CD3, pctCD45_Mono) %>% head(8))
