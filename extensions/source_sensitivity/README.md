# Post hoc RSEM and statistical-model sensitivity

From the repository root:

    python extensions/source_sensitivity/prepare_input.py
    Rscript extensions/source_sensitivity/reproduce.R extensions/source_sensitivity

The R environment must provide edgeR, limma and statmod. The checked run used R 4.6.1, edgeR 4.10.5 and limma 3.68.5; see the emitted `results/运行环境.txt`. The input downloader obtains GSE235046_Count_table.txt.gz from the GEO accession and verifies its fixed SHA-256. The assets are the project-candidate definition, GEO sample map, mouse/human traceability, and the previously retained 101-gene edgeR result. Expected: each of the fractional edgeR, rounded edgeR and fractional limma-voom branches retains the same 53 of 101 candidates, including all 48 frozen historical members.

These are post hoc comparisons based on the same 15 libraries and common 12,705-gene filter. No external or functional validation is claimed. The available count table does not provide transcript-abundance/effective-length inputs for a gene-length offset, and donor or culture independence remains unresolved. Rounding is a numeric sensitivity test, not a recommendation. Individual-contrast FDRs do not imply 5% joint FDR for their intersection.
