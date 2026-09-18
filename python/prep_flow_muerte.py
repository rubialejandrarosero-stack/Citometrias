#!/usr/bin/env python3
"""Tidy the death / immune-blockade cytometry table (rectified xlsm) into long CSVs.
Gating: Total -> Zombie- (live) -> {CD45- spheroid -> EpCAM+ -> PD-L1+ ; CD45+ PBMC -> CD3+ -> PD-1+}.
Outputs (output/tidy/): viabilidad.csv, pdl1_flow.csv, epcam_flow.csv, pd1_flow.csv, cd45_counts.csv"""
import pandas as pd, numpy as np, os

SRC = "/mnt/c/Users/57319/Downloads/IDCBIS/FUCS/GRAFICAS TESIS/resultados_citometria experimento muerte y bloqueo inmune.xlsm"
OUT = "output/tidy"; os.makedirs(OUT, exist_ok=True)

raw = pd.read_excel(SRC, sheet_name="RESULTADOS VIABILIDAD", header=0, engine="openpyxl")
# select by position (headers have irregular spacing)
c = {i: raw.columns[i] for i in range(raw.shape[1])}
df = pd.DataFrame({
    "tipo": raw[c[2]], "amb": raw[c[3]], "act": raw[c[4]], "time": raw[c[1]], "donor": raw[c[5]],
    "viability": pd.to_numeric(raw[c[8]], errors="coerce"),   # Celulas_Vivas %
    "cd45neg": pd.to_numeric(raw[c[10]], errors="coerce"),    # Esferoide CD45- (live count)
    "pdl1":    pd.to_numeric(raw[c[12]], errors="coerce"),    # PDL1_Pos % (rectified in place)
    "epcam":   pd.to_numeric(raw[c[13]], errors="coerce"),    # exp_EpCAM_Pos %
    "cd45pos": pd.to_numeric(raw[c[16]], errors="coerce"),    # PBMC CD45+ (live count)
    "pd1":     pd.to_numeric(raw[c[18]], errors="coerce"),    # exp_PD1_T %
})
df["environment"] = df.amb.map({"INF": "Inflammatory", "NO_INF": "Basal"})
df["condition"] = np.where(df.tipo == "Esferoide_Solo", "Control",
                    np.where(df.act == "NO_ACT", "Resting", "Activated"))
df["donor"] = df.donor.where(df.donor.isin(["D1", "D2"]), pd.NA)
df["time"] = pd.to_numeric(df.time, errors="coerce").astype(int)

base = ["time", "environment", "condition", "donor"]
df[base + ["viability"]].to_csv(f"{OUT}/viabilidad.csv", index=False)
df[base + ["pdl1"]].to_csv(f"{OUT}/pdl1_flow.csv", index=False)
df[base + ["epcam"]].to_csv(f"{OUT}/epcam_flow.csv", index=False)
df[base + ["pd1"]].to_csv(f"{OUT}/pd1_flow.csv", index=False)

# CD45 viable counts in long format (two lines: spheroid vs PBMC)
cd45 = df.melt(id_vars=base, value_vars=["cd45neg", "cd45pos"],
               var_name="population", value_name="count")
cd45["population"] = cd45.population.map({"cd45neg": "CD45- (spheroid)", "cd45pos": "CD45+ (PBMC)"})
cd45.to_csv(f"{OUT}/cd45_counts.csv", index=False)
print("tidy CSVs written to", OUT, "| rows:", len(df))
