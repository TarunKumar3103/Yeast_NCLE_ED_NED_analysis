# Yeast NCLE–poor-exponential-fit association analysis

Tests whether *S. cerevisiae* proteins with a native non-covalent lasso entanglement (NCLE) are more likely to have a poor single-exponential degradation fit, using turnover data from Christiano et al., *Cell Reports* 2014. The repo contains the full pipeline, frozen input snapshots, results, and two independent verifications.

**Outcome definition.** NED = R² of the single-exponential fit below the cutoff; ED = R² at or above it. The primary cutoff is 0.90, with sensitivity cutoffs at 0.75, 0.85 and 0.95. This is an operational poor-fit label. It is not McShane's NED classification used in the human and mouse analyses, and it does not establish bi-exponential kinetics or misfolding, so report it separately rather than pooling it with those species.

## Results

Final cohort: 2,072 proteins (1,565 NCLE-positive, 507 NCLE-negative).

| R² cutoff | NED proteins | NCLE OR (length-adjusted) | 95% CI | P |
|---|---|---|---|---|
| 0.90 (primary) | 43 | 0.60 | 0.31–1.17 | 0.135 |
| 0.75 | 9 | 0.83 | 0.19–3.64 | 0.801 |
| 0.85 | 19 | 0.54 | 0.20–1.43 | 0.212 |
| 0.95 | 80 | 0.72 | 0.43–1.19 | 0.203 |

No statistically significant NCLE association was detected at any cutoff, and every CI includes 1. The 0.75 cutoff has only nine events and is very imprecise. Crude Fisher results are in `results/association_results.csv`.

**Methods in brief.** Model 1 is logistic regression `NED ~ NCLE + z(protein length)`. The crude analysis uses a Fisher exact P with a Woolf CI. Cohort filters, in order:

- single-accession turnover record with an exact match to the native annotation table
- supplied mean pLDDT ≥ 70
- not knotted, and no covalent lasso (AlphaLasso)
- no UniProt Membrane (KW-0472) or Secreted (KW-0964) keyword

Full details are in `README.txt`.

## Reproduce

Python 3.12 and the pinned packages in `requirements.txt` (statsmodels is needed only for the independent check). Everything runs offline from the frozen files in `raw_sources/`.

~~~
python -m pip install -r requirements.txt statsmodels

python scripts/run_yeast_r2.py \
  --turnover raw_sources/Christiano_yeast_mmc2.xlsx \
  --native raw_sources/AF_Yeast_entanglement.pkl \
  --knots raw_sources/knot_all.txt \
  --alphalasso raw_sources/AlphaLasso.zip \
  --keywords raw_sources/uniprot_S_cerevisiae_keyword_snapshot_2026_10_10.tsv \
  --primary-cutoff 0.90 --cutoffs 0.90 0.75 0.85 0.95 \
  --output-dir results/local_reproduction
python scripts/verify_results.py \
  --results results/local_reproduction \
  --raw-sources raw_sources \
  --output verification/local_reproduction
python scripts/independent_checks.py --yeast .
~~~

Expected output:

- `verify_results.py` prints `PASS: 88 numerical fields agree`.
- `independent_checks.py` reports **15/15 independent checks passed**.

`independent_checks.py` rebuilds the cohort from the raw files without importing any package code and refits the models with statsmodels. It also confirms that your rerun matches the shipped tables. Last-digit float rounding differences between platforms are reported but accepted. To rerun, delete the `local_*` output folders first; the scripts refuse to overwrite them.

## Layout

| Folder | Contents |
|---|---|
| `scripts/` | Pipeline (`run_yeast_r2.py`, `stats_model1.py`), verifiers |
| `raw_sources/` | Frozen inputs: Christiano supplement, native NCLE table, knots, AlphaLasso, UniProt keyword snapshot |
| `results/` | Shipped result tables |
| `verification/` | Shipped verification outputs |
| `Yeast_NCLE_ED_NED_R2_results.xlsx` | Formatted workbook |

`PACKAGE_MANIFEST.json` lists SHA-256 hashes for every file.

## Limitations

- **R² measures fit noise, not just decay shape:** R² is strongly tied to the slope's coefficient of variation (Spearman ρ = −0.88). Poor-fit proteins also cluster among the fastest degraders, while NCLE proteins turn over more slowly. Adding log(rate) to the model as a sensitivity check leaves the NCLE OR at about 0.57 (P = 0.12).
- **Grouped accessions:** the 100 grouped multi-accession records are left unresolved and excluded.
- **Keyword snapshot:** the UniProt keywords are a 2026_03 snapshot.
- **Interpretation:** the ORs are associations, not causal effects.

## Sources

Christiano et al., *Cell Reports* 9, 1959–1965 (2014), https://doi.org/10.1016/j.celrep.2014.10.065 · UniProt, https://www.uniprot.org
