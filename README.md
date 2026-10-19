# hcc-signature-audit

Code and source data for "Random signatures reproduce the prognostic findings of published hepatocellular carcinoma cell death signatures": an audit of 149 published HCC cell-death prognostic signatures against size-matched and input-matched random gene signatures, followed by a two-gate audit framework.

Citation: [to be added]

Several downstream modules re-implement analysis steps of the template workflow commonly used in the audited literature (see Methods).

## Repository layout

```
analysis/              bulk, single-cell and simulation analyses (R, Python, shell)
corpus/                literature search, full-text retrieval and coefficient extraction (Python, R)
figures/               figure, Source Data and supplementary-table scripts (R)
  data/                small intermediate tables read by the figure scripts (see "Figures" below)
data/genesets/         immune gene sets used by modules M2 and M13
source_data/           Source Data for every main and supplementary figure panel (TSV)
  TableS3_numbers.tsv  every number quoted in the text, with its source file and extraction rule
  per_signature/       per-signature module calls and gate results for the 148 scorable signatures
supplementary_tables/  Supplementary Tables S1-S3 (xlsx)
env/                   R package versions, R sessionInfo(), Python package versions
```

## Paths

All scripts resolve paths from environment variables, so nothing needs to be edited:

| Variable | Meaning | Default |
|---|---|---|
| `HCC_AUDIT_ROOT` | working root; scripts are read from `analysis/` and `figures/` under it, and outputs are written to `N11_out/`, `N12_signatures/` and `N28_fig/` under it | `.` (run from the repository root) |
| `HCC_AUDIT_DATA` | raw public data (layout below) | `data/raw` |
| `HCC_AUDIT_RLIB` | optional R library path | unset |
| `GENEFORMER_MODEL_DIR` | Geneformer V1-10M pretrained model | `models/Geneformer-V1-10M` |
| `RSCRIPT`, `PYTHON` | interpreters used by the two shell run chains | `Rscript`, `python` |

Use absolute paths for these variables if scripts are launched from another working directory (the shell chains in `analysis/*.sh` change directory before running R).

## Data sources

All data are public. Place them under `$HCC_AUDIT_DATA` as follows.

| Dataset | Source | Expected files |
|---|---|---|
| TCGA-LIHC RNA-seq (RSEM, log2(x+1)), survival, phenotype, MC3 somatic mutations | UCSC Xena, TCGA hub (https://xenabrowser.net) | `TCGA/TCGA-LIHC_gene_expression.tsv.gz`, `TCGA/TCGA-LIHC_survival.tsv`, `TCGA/TCGA-LIHC_phenotype.tsv`, `TCGA/LIHC_mc3.txt.gz` |
| GSE14520 (GPL3921) | GEO (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE14520) | `GEO_bulk/GSE14520/GSE14520-GPL3921_series_matrix.txt.gz`, `GEO_bulk/GSE14520/GPL3921.annot.gz`, `GEO_bulk/GSE14520/GSE14520_Extra_Supplement.txt.gz` |
| LIRI-JP | HCCDB, dataset HCCDB18 (http://lifeome.net/database/hccdb) | `GEO_bulk/HCCDB18_LIRI-JP/HCCDB18_mRNA_level3.txt`, `HCCDB18.sample.txt`, `HCCDB18.patient.txt` |
| GSE149614 scRNA-seq | GEO (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149614) | `GSE149614/GSE149614_HCC.metadata.updated.txt.gz`; `GSE149614/GSE149614_counts_sparse.rds` (the GEO count matrix stored as a sparse genes x cells `dgCMatrix`) |
| GDSC2 training data for oncoPredict | oncoPredict data files (OSF repository linked from the oncoPredict documentation; `DataFiles.zip`) | `external/oncoPredict_osf/DataFiles.zip` |
| DSigDB | https://dsigdb.tanlab.org | `external/DSigDB/DSigDB_All_detailed.txt` |

Files read from `$HCC_AUDIT_ROOT`:

| File | Source |
|---|---|
| `N11_geneset/h.all.v2024.1.Hs.symbols.gmt`, `c2.cp.reactome.v2024.1.Hs.symbols.gmt`, `c2.cp.wikipathways.v2024.1.Hs.symbols.gmt`, `c5.go.bp.v2024.1.Hs.symbols.gmt` | MSigDB v2024.1 (https://www.gsea-msigdb.org) |
| `N12_signatures/ref/hgnc.txt` | HGNC complete set (https://www.genenames.org/download/archive/) |
| `N12_signatures/signatures_long.tsv`, `signatures_index.tsv` | produced by `corpus/` (contents are in Supplementary Table S1) |

Gene sets included in this repository (`data/genesets/`):

| File | Content |
|---|---|
| `immune29_ssgsea.gmt` | the 29 immune signatures (immune cell types and immune functions) commonly used for ssGSEA immune-infiltration scoring, from He et al. 2018, *J Exp Clin Cancer Res* 37:327, "Classification of triple-negative breast cancers based on immunogenomic profiling" (PMID 30594216; doi:10.1186/s13046-018-1002-1); used by module M2 |
| `immune_checkpoint_genes.txt` | 47 immune checkpoint genes, one symbol per line, used by module M13 (C10orf54 recorded under its current HGNC symbol VSIR) |

## Software

R 4.6.1; package versions are listed in `env/software_versions_R.tsv` and `env/sessionInfo_R.txt`. The Geneformer step (module M10) used Python 3.11 with Geneformer 0.1.0, PyTorch 2.11.0 and transformers 4.49.0 (`env/software_versions_python.txt`).

## Running order

Scripts are listed in the order of the Methods sections. Each script header states its inputs, outputs and decision rules. Several scripts accept environment variables for smoke runs (for example `R4_SMOKE=1`, `M1315_SMOKE=1`); the defaults reproduce the reported analyses.

### 1. Literature corpus and signature extraction (`corpus/`, run inside `$HCC_AUDIT_ROOT/N12_signatures`)

1. `search.py` - PubMed searches (one per cell-death programme)
2. `meta.py` - PubMed metadata and title screening
3. `fetch.py` - PMC full text (BioC XML, JATS fallback)
4. `totxt.py` - full text to plain text
5. `scan.py`, `parse2.py` - coefficient parsing (`hgnc_util.py` maps symbols to HGNC)
6. `authors.py` - author lists for author-group collapsing
7. `build.py`, `write_out.py` - signature tables
8. `module_frequency_fulltext.R` - use of downstream analyses in the corpus (Fig. 1b)

`parse.py` and `extract.py` are earlier parsing passes whose outputs are not read by `build.py` (see Methods).

### 2. Cohorts, gene harmonization, scores and proxies

- `analysis/load_cohorts.R` - loads the three bulk cohorts, harmonizes symbols, computes proxy scores and defines the template pipeline (univariable Cox, then LASSO-Cox). It is sourced by the scripts below and is not run on its own.

### 3. Cell-death gene sets (Fig. 2)

- `analysis/robust_proxy.R` - 50 internal splits, external validation and proxy-adjusted signal for the eight gene sets
- `analysis/adj_null.R` - proxy-adjusted random null in the external cohorts

### 4. Proxy analysis (Fig. 3)

- `analysis/ring2_proxy.R` - per-signature proxy correlations, R² and residual likelihood-ratio tests
- `analysis/ring2_analyze.R` - published versus random comparisons
- `analysis/n4_loo_proxy.R` - proxies recomputed without the signature genes

### 5. Template downstream modules (Fig. 4)

- `analysis/ring3_lib.R` - module functions (sourced)
- `analysis/ring3_full.R` - modules M0-M2 and M4-M7 for published signatures
- `analysis/ring3_armA_realrescale.R` - published signatures scored with coefficient x TCGA s.d.
- `analysis/ring3_armA_null.R` - size-matched and input-matched random signatures
- `analysis/ring3_analyze.R`, `analysis/ring3_analyze_armA_v3.R` - module positive rates and decisions (`ring3_analyze_armA.R` and `ring3_analyze_armA_v2.R` are earlier versions of v3)
- `analysis/hiend_m3_m12.R` - drug sensitivity (M3) and DSigDB enrichment (M12)
- `analysis/hiend_m13_m15.R` - immune checkpoints (M13), consensus clustering (M14), calibration and decision curves (M15)
- `analysis/cmh_modules.R` - Cochran-Mantel-Haenszel tests across the 16 modules
- `analysis/cmh_modules_fullcov95.R` - the same tests restricted to the 95 signatures fully covered in all three cohorts (sensitivity analysis; read by `figures/figS_supp.R`)
- `analysis/n3_cohort_use_strata.R` - external validation stratified by external-cohort use in the source article (descriptive)

### 6. Single-cell and in silico perturbation modules (Fig. 6)

- `analysis/hiend_m8_sc_pilot.R` - per-gene single-cell marker and tumour-versus-normal tests (M8a, M8b); despite its name, its per-gene output is the input of the formal call
- `analysis/hiend_m8_formal.R` - per-signature M8 calls
- `analysis/hiend_m10_gf_export.R` -> `analysis/hiend_m10_gf.py` -> `analysis/hiend_m10_lookup.py` - Geneformer in silico deletion (M10)
- `analysis/hiend_m9_knk_official_check.R` - agreement of the scTenifoldKnk implementation with the reference function
- `analysis/hiend_m9_knk_hvg.R`, `analysis/hiend_m9_discrim_check.R`, `analysis/hiend_m9_make_kolist.py`, `analysis/hiend_m9_lookup.py` - scTenifoldKnk virtual knockout; run chains `analysis/hiend_m9_chain_formal.sh` and `analysis/hiend_m9_chain_hep.sh`

### 7. Matched null references

- `analysis/hiend_m12_annot_control.R` - DSigDB annotation-matched null (M12)
- `analysis/hiend_m10_abundance.py` - detection-abundance-matched null (M10)
- `analysis/hiend_m14_prolif_check.R` - proliferation separation of consensus-clustering subtypes (M14)

### 8. Simulation and two-gate audit (Fig. 7)

- `analysis/ring4v2_selected_sim.R` - plasmode simulation with selected candidate signatures and the misspecified-null scenario
- `analysis/n2_published_gates.R` - the 148 published signatures through both gates
- `analysis/ring4_positive_control.R` - earlier simulation design superseded by `ring4v2_selected_sim.R` (see Methods)

### 9. Figures, Source Data and supplementary tables (`figures/`)

`00_config.R` holds shared settings and is sourced by every figure script.

1. Preparation: `prep_sc_embedding.R`, `prep_fig4c_random_showcase.R`, `check_fig5_showcase_typicality.R`
2. Panels: `fig1_design.R`, `fig2_celldeath_sets.R`, `fig3_proxy.R`, `fig4_modules.R`, `figA_funnel.R`, `figX_landscape.R`, `fig5_singlecell.R`, `figCD_sc_extra.R`, `fig7_framework_v2.R`, `figS_supp.R`, `figS5_showcase_typicality.R`
3. Assembly: `assemble_main.R` (main Figs. 4-6)
4. Numbers and tables: `make_numbers.R` (`source_data/TableS3_numbers.tsv`), `make_supp_tables.R` (Supplementary Tables S1-S3 and `source_data/per_signature/`)
5. `session_versions.R` - writes `env/` version files

Intermediate tables in `figures/data/` (written by the preparation scripts, included so that the panel scripts can run without the raw data):

| Files | Written by | Read by |
|---|---|---|
| `fig4c_km_curves.tsv`, `fig4c_km_logrank.tsv`, `fig4c_multicox.tsv`, `fig4c_timeroc_auc.tsv`, `fig4c_timeroc_curves.tsv` | `prep_fig4c_random_showcase.R` | `fig4_modules.R` |
| `fig5_showcase_all_Ak.tsv`, `fig5_showcase_typicality.tsv` | `check_fig5_showcase_typicality.R` | `figS5_showcase_typicality.R` |
| `sc_umap.tsv.gz`, `sc_sigscore_coverage.tsv` | `prep_sc_embedding.R` | `fig1_design.R`, `fig5_singlecell.R`, `figCD_sc_extra.R`, `make_numbers.R` |

R objects (`*.rds`) are not included. `fig5_singlecell.R` and `figCD_sc_extra.R` also need `figures/data/sc_sigscore.rds`, so run `prep_sc_embedding.R` (requires GSE149614) first. `assemble_main.R` needs the panel objects saved by `fig4_modules.R`, `fig5_singlecell.R`, `figX_landscape.R`, `figA_funnel.R` and `figCD_sc_extra.R`, and `figS_supp.R` needs the object saved by `figX_landscape.R`; run those scripts first.

Script names follow the order in which panels were drafted, not the final figure numbers; the header of each script names the panels it produces. `fig6_framework.R` is an earlier version of `fig7_framework_v2.R`.

### Supplementary figure numbering

Supplementary figure files and Source Data keep their internal names; the numbering in the manuscript follows the order of first citation:

| Manuscript | Output name in this repository (`figS_supp.R`, `figS5_showcase_typicality.R`) | Source Data prefix |
|---|---|---|
| Supplementary Fig. S1 | `FigS4_gene_reuse` | `FigS4*` |
| Supplementary Fig. S2 | `FigS2_consensus_proliferation` | `FigS2*` |
| Supplementary Fig. S3 | `FigS3_singlecell_k_dependence` | `FigS3*` |
| Supplementary Fig. S4 | `FigS1_fullcoverage_sensitivity` | `FigS1*` |
| Supplementary Fig. S5 | `FigS5_showcase_typicality` | `FigS5*` |

The `fig` column of `source_data/TableS3_numbers.tsv` uses the manuscript numbering.

### Exploratory scripts (not used for reported results)

`analysis/pilot_proxy.R`, `analysis/ring3_pilot.R`, `analysis/ring3_pilot2.R`, `analysis/ring3_pilot3_published.R`, `analysis/ring3_m4time.R`, `analysis/n16_gate_pilot.R`, `analysis/n16_znf180_check.R`, `analysis/hiend_m9_knk_pilot.R`, `analysis/hiend_m9_knk_shared.R`, `analysis/hiend_m9_knk_validate.R`, `analysis/eqprobe.R`, `analysis/eqprobe2.R`, `analysis/probe_oncopredict.R`.

These were feasibility checks, timing runs or superseded implementations. They are kept for transparency.

## Random seeds

All random draws use fixed seeds set inside the scripts; the seed scheme for each analysis is written in the script header or the corresponding analysis plan.

## License

MIT (see `LICENSE`).
