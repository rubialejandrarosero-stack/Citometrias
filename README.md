# Citometrias — Immune infiltration in 3D lung adenocarcinoma spheroids

Analysis code for the manuscript **"Immune Infiltration and Modulation in Three-Dimensional
Lung Adenocarcinoma Spheroids Co-cultured with Human PBMCs"** (A549/MRC-5 heterotypic spheroids
co-cultured with human PBMCs). All figures in the article and supplementary material are generated
by the scripts in this repository from the raw cytometry, microscopy, histology and Luminex data.

## Environments (conda / Miniforge, conda-forge)

Two environments are used:

```bash
# R — flow cytometry (Bioconductor: flowCore, ggcyto, flowWorkspace, openCyto) + tidyverse
conda env create -f environment.yml            # env name: citometrias
conda activate citometrias

# Python — image processing, figure composition, secretome/PCA/volcano, stats
conda env create -f environment-py.yml         # env name: citometrias-py
conda activate citometrias-py
```

## Repository layout

```
Citometrias/
├── data/raw/            # raw inputs (.fcs, Luminex .xlsm, microscopy) — NOT versioned (heavy)
├── R/                   # R analysis + per-panel plotting (flowCore, tidyverse, lmerTest)
├── python/              # image processing + figure composition (Pillow, numpy, sklearn, scipy)
├── output/
│   ├── tidy/            # processed/tidy data (CSV) — reproducible without raw data
│   ├── figuras/         # generated panels
│   │   └── submission/  # final TIFFs (300 dpi, RGB) used in the manuscript
│   └── manuscrito/      # working copy of the manuscript + reference-update notes
├── environment.yml      # R environment
├── environment-py.yml   # Python environment
└── README.md
```

## Reproducing the figures

Panels are produced in R (`R/*.R`), then composed into the final multi-panel figures in Python
(`python/compose_*.py`), and finally exported to Frontiers-ready TIFFs (300 dpi, RGB, 85/180 mm)
with `python/make_submission.py <figure_name> <width_mm>`.

### Figure → script map (final figures in `output/figuras/submission/`)

| Manuscript figure | Final file | Composition script | Upstream panel / data scripts |
|---|---|---|---|
| **Fig 1** Morphology + H&E + morphometrics | `Fig1_morphology_full.tiff` | `python/compose_fig1.py` | `R/31_fig_morfologia_dotplot.R`, `R/16_estadistica_morfologia.R`, `python/compose_morphology.py`, `python/img_norm.py` |
| **Fig 2** Ki-67 IHC | `Fig_IHC_Ki67.tiff` | `python/compose_ihc_full.py` | `R/40_fig_ihc_quant.R`, `python/compose_pptx_ihc.py`, `python/img_norm.py` |
| **Fig 3** Confocal infiltration | `Fig4_confocal_infiltration.tiff` | `python/compose_fig4.py` | `R/12_fig_confocal.R`, `R/17_fig_confocal_composicion.R`, `R/19_fig_confocal_radial_overlay.R`, `python/confocal_radial_analysis.py` |
| **Fig 4** CD3/CD14 IHC | `Fig_IHC_CD14CD3.tiff` | `python/compose_ihc_full.py` | `R/40_fig_ihc_quant.R`, `python/compose_pptx_ihc.py` |
| **Fig 5** Flow immune populations | `Fig6_populations.tiff` | `python/compose_fig6.py` | `R/27_fig_citometria_dotplot.R`, `R/13_fig_citometria.R` |
| **Fig 6** Immune/macrophage activation | `Fig7_activation.tiff` | `python/compose_fig7.py` | `R/34_fig_activacion_dotplot.R`, `R/30_dotplots_representativos.R` |
| **Fig 7** Secretome (Luminex) | `Fig_secretome_full.tiff` | `python/compose_secretome_full.py` | `python/compose_luminex_heatmap.py`, `python/compose_luminex_pca.py`, `python/compose_luminex_volcano.py`, `R/43_fig_cytokines.R`, `R/14_fig_luminex.R` |
| **Fig 8** PD-L1 IHC | `Fig_IHC_PDL1.tiff` | `python/compose_ihc_full.py` | `R/40_fig_ihc_quant.R`, `python/compose_pptx_ihc.py` |
| **Fig 9** Viability & checkpoints (flow) | `Fig_viability_full.tiff` | `python/compose_fig_viability_full.py` | `R/44_fig_flow_muerte.R`, `python/prep_flow_muerte.py`, `python/compose_gating_viability.py` |

### Supplementary figures

| Supplementary | Final file | Script |
|---|---|---|
| **Fig S1** Gating strategy (immunophenotyping) | `FigS_gating_composition.tiff` | `python/compose_gating_composition.py`, `python/compose_gating.py` |
| **Fig S2** Gating strategy (viability/blockade) | `FigS_gating_viability.tiff` | `python/compose_gating_viability.py` |
| **Fig S3** Brightfield spheroid grid | `FigS_spheroids_brightfield.tiff` | `python/compose_spheroid_supp.py` |
| **Fig S4** H&E spheroid grid | `FigS_spheroids_HE.tiff` | `python/compose_spheroid_supp.py` |
| **Fig S5–S7** Confocal matrices (24/48/96 h) | `FigS_confocal_{24,48,96}h.tiff` | `python/compose_confocal_montage.py` |

### Data ingestion / tidying
- `R/01_leer_fcs.R`, `R/10_infiltracion_gating.R` — read `.fcs`, gating, infiltration counts.
- `python/prep_flow_muerte.py` — tidy the death/immune-blockade cytometry table.
- `R/14_fig_luminex.R` — Luminex secretome ingestion.
- Tidy outputs are written to `output/tidy/*.csv`.

## Notes
- Statistical analyses: R (mixed-effects models via `lmerTest`, non-parametric tests) and Python
  (PCA via scikit-learn, Mann–Whitney volcano via SciPy). See the manuscript Methods.
- Raw `.fcs`/Luminex files are not versioned (size/privacy); the tidy CSVs in `output/tidy/`
  allow regenerating the figures from processed data.
