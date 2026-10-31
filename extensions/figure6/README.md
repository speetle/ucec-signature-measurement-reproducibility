# Figure 6: exploratory perturbation comparisons

Run from the repository root with Python, NumPy, pandas, SciPy and Matplotlib:

    python extensions/figure6/run.py

The entry point checks or downloads the GSE316734 count matrix from NCBI GEO with a fixed SHA-256, runs the two sample-score calculations and the 5,000 reference draws, ranks the 92 Tahoe conditions, verifies the displayed summaries and writes PDF/PNG images in `figures/`. The external input and outputs are intentionally ignored by git. It uses the frozen CC0 HGNC table already available in the root `input/hgnc.txt`.

## Frozen inputs and provenance

`assets/GSE73024_AN3CA基因表达与响应.csv`: gene-level data derived from the published GSE73024 processed GC-RMA matrix, four sample columns and the full combinations, using GPL570 annotation and the frozen HGNC table. Nakanishi et al., Mol Cancer Ther 2015, PMID 26438159; original accessions GSM1880917–GSM1880920. This is not a CEL-level rebuild. Data source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73024.

`assets/Tahoe板2_AN3官方响应.npz` and `assets/Tahoe板2官方响应完整性.csv`: project extraction of cell-line CVCL_0028, plate 2, from the Tahoe-100M public pseudobulk differential-expression files at dataset revision `349e2a55f5538a3e73b32a06c6bdd8e9e87011e8`; source: https://huggingface.co/datasets/tahoebio/Tahoe-100M . The data source is described as CC0-1.0 in the 2026 Cell article data-availability statement. Cite Zhang et al., Cell 2026, PMID 42753697. Full source count generation and the model behind this published differential table have not been reproduced here. Per-condition cell counts and selected values are retained as downloaded, with no claim of significance or mechanism from these tables.

`assets/Tahoe源48及共有29基因映射.csv` is the fixed signed project panel definition. `GSE316734_Raw_count_file_rsemgenes.txt.gz` is retrieved from GEO, not redistributed. It supplies the six samples, and GEO labels them biological replicate cultures; exact permutation assumes exchangeable labels. Study: Huang et al., PMID 41873535.

The GEO processed numerical matrix is provided for reproducibility with source attribution; the original publications and third-party materials are not relicensed by this repository. NCBI data use guidance and original study conditions apply. Frozen input and script hashes appear in `MANIFEST.json`. These are computational tests on the same cell line and distinct compounds, not validation of protein measurement, cell death or clinical benefit.
