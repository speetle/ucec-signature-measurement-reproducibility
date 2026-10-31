# Additional analyses for the Chinese measurement manuscript

These scripts extend the public v1.0.1 processed-data reproducibility archive. They reproduce the later ordering and RPPA summary, Figure 6, and the source-signature sensitivity analysis. They do not re-run primary sequencing, mass-spectrometry quantification, establish biological ground truth, or establish clinical utility. Analyses were performed after the original study results were known.

## Ordering and RPPA aggregate

Run `python reproduce.py` first to generate the original v1.0.1 `results/null` and `results/rppa` inputs. Then run `python extensions/v17_review/reproduce_v17.py`. The latter needs the four prior outputs listed in `v17_review/README.md`; it is a downstream rerun and does not replace the v1.0.1 acquisition/processing workflow. The independent CPTAC portal returned HTTP 403 in one check; if retrieval is unavailable, provide the original permitted version-3.0 inputs locally and verify the SHA-256 listed in `external_sources.json`. The original RPPA download is also subject to its provider's availability and terms.

## Figure 6

Follow `figure6/README.md`. `python extensions/figure6/run.py` downloads and verifies the GEO GSE316734 source file, calculates the six sample scores and 5,000 nested reference panels, reconstructs the 92-condition retrieval from frozen response summaries, and writes the figure. The GSE73024 and Tahoe data enter at processed response matrices. Tahoe raw single-cell counts and the generation of the published differential expression table are not independently reconstructed or validated here.

## Source-signature sensitivity

Follow `source_sensitivity/README.md`. A GEO input download, three alternative statistical runs and exact comparison against the original 101-candidate audit are supplied. This is a post hoc sensitivity comparison on the same samples. It cannot prove independent cultures, remove transcript-usage bias, or validate cell-death mechanisms.

Code is covered by this repository's MIT license. Third-party data are subject to the individual source terms in `THIRD_PARTY_DATA.md` and `figure6/README.md`.
