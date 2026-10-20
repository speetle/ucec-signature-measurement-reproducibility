# UCEC signature measurement correspondence

This release examines a bounded measurement question: when patients and a partial protein endpoint are fixed, how does restricting RNA membership change ordering correspondence? Positive correspondence changes do not establish biological construct accuracy, cell-death activity, assay interchangeability or clinical usefulness.

## Reproduce

Use Python 3.9+ with the pinned requirements, R and survival 3.8.6. The tested Python numerical environment is in environment.json. The code uses two computational threads per process. No fee-based service or wet experiment is required. Runtime depends on hardware; the included validation is same-host computational verification, not independent hardware replication.

    python -m pip install -r requirements.txt
    python reproduce.py --fetch-missing

The public archive obtains the versioned LinkedOmics RNA, protein and metadata, TCPA RPPA matrix and publisher antibody workbook directly from their original sources. Source SHA-256 checks fail on changed versions instead of silently substituting new data. Existing full local inputs require no network:

    python reproduce.py
    python reproduce.py --verify-only

Figure-only use requires prior complete analysis outputs. The complete entry requires R; it does not silently substitute frozen clinical tables.

    python reproduce.py --figures-only

## What is recalculated

The single entry recomputes primary 20,000-draw patient contrasts, weight controls and batch mismatch, three structural references, 29 member deletions, hypothetical LOXL3 scenarios, actual MS benchmarks, 278,400 observation-mask rows, independent mapped-array projections, RPPA mapping and patient pairing, nested held-out prediction, lineage-count scores and exact donor tests, detailed/simplified adjustment, source48 clinical 500 refits, and association/PH diagnostics. Every displayed panel is reconstructed from those results. Numeric comparison files and stage logs identify the actual scope and failures.

Input layers remain explicit. RNA/protein matrices are publisher-processed or frozen gene-mapped analysis inputs. Cell inputs are raw-count and corrected-count lineage aggregates plus recorded QC/annotation evidence, not raw single-cell reads. The original upstream QC and marker-annotation scripts are supplied for inspection in upstream_reference; the entry does not rerun them. Array input is the deposited-scale, gene-mapped matrix, not CEL reconstruction. Clinical analyses start from adjudicated individual analysis records. No raw FASTQ alignment, raw MS quantification, new donor-identity mapping, measured AOI geometry, direct death assay or prospective clinical utility test is claimed.

Frozen definitions (gene sets, signed membership, observation lists, model covariates and endpoint adjudication) are part of the analysis specification. Frozen expected tables are validation references. They are not used as substitutes for the re-run displayed statistics. Historical non-core analyses are outside this minimal display release.

## Outputs

- results/primary: same-patient fixed-endpoint contrasts and individual scores.
- results/null and results/controls: weight, mismatch, structural matching, deletion and assumptions.
- results/ms and results/rppa: pathway comparisons and complete eligibility records.
- results/masks and results/arrays: measurement/background stress and actual RNA-array projections.
- results/heldout and results/cells: score/composition descriptions and exact donor inference.
- results/clinical: clinical increment and association/PH checks, with negative findings retained.
- results/review: matching balance/error strata, all 50 pathway qualification records and all adjustment designs.
- results/figures: five main multipanel figures, PDF/SVG/TIFF/PNG, source hashes and alignment reports.

The 31 external eligible Hallmarks and 19 unavailable sets are both retained; 11 eligible contrasts are negative. RPPA covers other pathways and zero source48 members. No clinical increment was established. These are findings, not errors to tune away.

## Provenance and attribution

See input_manifest.json, external_sources.json, optional_sources.json, THIRD_PARTY_DATA.md and panel_registry.json. A processed-data reconstruction is not a reconstruction of every raw experiment. The manuscript and author-specific declarations are maintained separately; this code repository does not supply or invent author affiliations, ethics determinations, grants, conflicts, journal decisions or a data DOI.
