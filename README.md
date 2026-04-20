# rocs-repro

Reproducible pipeline for the ROCS marine sedaDNA analyses.
Focuses on WGCNA co-expression modules and HMM ecological states,
with EMP/SAP metrics and a comparison against TEA indices.

## Pipeline

```
01_data_prep.R    Filter damaged reads → DESeq2/VST normalisation
02_wgcna.R        Consensus WGCNA across ST8/ST13/GeoB-R1; preservation test on GeoB-R2
03_hmm_states.R   HMM (K=5) on module eigengene PCA; glacial/interglacial states
04_emp.R          EMP = Σ(VST × module_completeness × |ΔG|); SAP ratio
05_tea_vs_emp.R   TEA indices (OAP, DCI, MII, SRPI, MGI, MFI) + comparison to EMP
```

## Quick start

```bash
bash run_all.sh           # all steps
bash run_all.sh --start 3 # resume from HMM
```

## EMP vs TEA

| Metric | Formula | What it captures |
|--------|---------|-----------------|
| EMP | Σ(VST × comp × \|ΔG\|) | Aggregate thermodynamic work capacity |
| OAP | Σ(rel_ab × E0') | Mean electron acceptor redox quality |
| DCI | log(nosZ / narG+norB) | Denitrification completeness |
| MII | log(bd-oxidase / aa3-oxidase) | Microoxic niche |
| SRPI | (log(M00596) + log(aprAB+dsrAB)) / 2 | Sulfate reduction potential |
| MGI | log(mcrABG completeness) | Methanogenesis |
| MFI | log(pmoA / mcrABG) | Methane filter |
