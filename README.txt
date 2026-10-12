YEAST NCLE–NED ASSOCIATION: QUYEN-APPROVED R² CLASSIFICATION
Analysis prepared 10 October 2026 (America/New_York).

QUESTION AND OUTCOME
Are native NCLE-positive proteins more likely to have poor single-exponential
degradation fits, after controlling for protein length?
Quyen specified ED = R² >= cutoff; NED = R² < cutoff.
The primary cutoff is 0.90. Exactly the cutoff is ED.
Sensitivity cutoffs are 0.75, 0.85 and 0.95, using the same eligible proteins.
These operational NED labels do not establish bi-exponential decay,
age-dependent degradation, or misfolding. No separate UN group is defined.

RESULTS
Species: Saccharomyces cerevisiae. Final cohort: 2,072 proteins.
NCLE-positive: 1,565. NCLE-negative: 507.

Cutoff  Analysis             Method          NED    NCLE OR       95% CI           P
0.90    Primary              Crude            43    0.664853    0.348520–1.268303  0.212379
0.90    Primary              Length-adjusted  43    0.601104    0.308232–1.172257  0.135283
0.75    Cutoff sensitivity   Crude             9    0.646568    0.161119–2.594675  0.463970
0.75    Cutoff sensitivity   Length-adjusted   9    0.825951    0.187212–3.643978  0.800655
0.85    Cutoff sensitivity   Crude            19    0.551927    0.216120–1.409509  0.279919
0.85    Cutoff sensitivity   Length-adjusted  19    0.536691    0.202003–1.425901  0.211925
0.95    Cutoff sensitivity   Crude            80    0.746852    0.457982–1.217926  0.235022
0.95    Cutoff sensitivity   Length-adjusted  80    0.719345    0.433211–1.194470  0.202962

CONTINGENCY COUNTS
          Analysis  R2_cutoff  NCLE  NED   ED  Total  NED_fraction
           Primary       0.90     1   29 1536   1565      0.018530
           Primary       0.90     0   14  493    507      0.027613
Cutoff sensitivity       0.75     1    6 1559   1565      0.003834
Cutoff sensitivity       0.75     0    3  504    507      0.005917
Cutoff sensitivity       0.85     1   12 1553   1565      0.007668
Cutoff sensitivity       0.85     0    7  500    507      0.013807
Cutoff sensitivity       0.95     1   56 1509   1565      0.035783
Cutoff sensitivity       0.95     0   24  483    507      0.047337

No statistically significant NCLE association was detected at any tested
cutoff. Every NCLE 95% CI includes 1. This is not proof of no association.
The 0.75 cutoff yields only nine NED proteins and very imprecise estimates.

METHODS
Crude: sample cross-product OR, Woolf log-OR Wald 95% CI, and two-sided
Fisher exact P-value. Model 1: maximum-likelihood logistic regression
NED ~ NCLE + z(Protein_length), with an intercept, logistic Wald 95% CI
and two-sided Wald P-values. Reference group: NCLE=0. Outcome=1: NED.
The normal critical value is scipy.stats.norm.ppf(0.975).
Length z-score: (length − 465.986969112) / 333.010409316.
The SD is the sample SD (ddof=1) of the 2,072 eligible proteins.
All models pass convergence, design-rank and separation checks.
No complex-membership covariate or disorder-fraction filter is used.
No multiple-testing adjustment is applied. All sensitivity results are
reported without choosing a threshold based on its P-value.

COHORT CONSTRUCTION
                                                Step  Remaining  Removed_at_step  NCLE_1  NCLE_0  NCLE_missing  NED_0_9  NED_0_75  NED_0_85  NED_0_95
                         All turnover source records       3772                0    2591    1068           113       96        17        45       181
Retain unambiguous single-accession turnover records       3672              100    2591    1068            13       96        17        45       179
     Require exact native-annotation accession match       3659               13    2591    1068             0       94        17        43       176
                     Retain supplied mean pLDDT >=70       2896              763    2125     771             0       62        13        29       116
                   Exclude supplied knotted proteins       2879               17    2111     768             0       62        13        29       116
            Exclude supplied covalent-lasso proteins       2852               27    2102     750             0       62        13        29       115
 Require coverage in frozen UniProt keyword snapshot       2852                0    2102     750             0       62        13        29       115
            Exclude UniProt Membrane keyword KW-0472       2089              763    1581     508             0       45         9        20        83
            Exclude UniProt Secreted keyword KW-0964       2072               17    1565     507             0       43         9        19        80

The supplied plddt field is interpreted as average AlphaFold pLDDT and
preserved without recomputation. The protein cutoff is pLDDT >=70.
NCLE comes from entanglement_AF (Yes=1, No=0), checked against whGLN.
Knot IDs and covalent-lasso IDs are excluded using the supplied sets.
Membrane and Secreted exclusions use frozen UniProt KW-0472 and KW-0964.
Source keyword release/query metadata and original inputs are included.
Single-accession exact matching is used. The 100 grouped records are
retained for review and are not expanded into independent observations.
Thirteen other accessions lack an exact native match. No isoform or
secondary-accession remapping is performed.
The source does not provide its AlphaFold release/date or producer code.
length_AF supplies the length covariate to retain the structural sequence
basis. P40516 has AF length 257 versus current UniProt 229; it is flagged
and the supplied 257 is preserved. Forty retained corrected rates are
nonnumeric (n.d.); their R² is available and defines the present outcome.

FILES AND WHICH CODE PRODUCES THEM
Yeast_NCLE_ED_NED_R2_results.xlsx:
  Results: primary and cutoff-sensitivity NCLE ORs, CIs, P-values,
    event counts, formula-linked contingency tables and interpretation.
  Model coefficients: intercept, NCLE and standardized-length effects,
    coefficients, standard errors, ORs, CIs, P-values and diagnostics.
  Cohort flow: sequential exclusions and NED counts at every cutoff.
  Protein records: all 3,772 turnover records with classifications,
    annotations, inclusion status, source row and first exclusion reason.
    Source definitions are alongside the data.
  Excel calculates contingency counts and crude OR/CIs. Fisher P-values
    and regression values are Python snapshots. Rerun Python after edits.

scripts/run_yeast_r2.py validates and joins the sources, applies
  eligibility rules, creates R²-based ED/NED labels, and writes:
  results/association_results.csv — crude and all regression terms.
  results/contingency_counts.csv — NED/ED counts and frequencies.
  results/eligible_yeast_cohort.csv — final 2,072 proteins.
  results/all_turnover_record_audit.csv — every experimental record.
  results/cohort_flow.csv — sequential counts and exclusions.
  results/grouped_accession_review.csv — 100 unresolved groups.
  results/native_annotation_export.csv — all 6,039 native entries.
  results/covalent_lasso_evidence.csv — matching archive evidence.
  results/model_diagnostics.csv — event counts, convergence diagnostics.
  results/analysis_manifest.json — methods, source hashes and runtime.
  results/workbook_payload.json — typed formatter inputs.
scripts/stats_model1.py implements crude and logistic estimates.
scripts/verify_results.py independently rebuilds the cohort from raw
  sources and refits logistic models with IRLS/Fisher scoring without
  importing the production statistical routines.
scripts/build_workbook.mjs formats the XLSX using the Codex artifact
  runtime. It is not needed to reproduce any statistical values locally.

VERIFICATION
Independent numerical fields compared: 88; all passed.
Every eligible ID, NCLE call, protein length, pLDDT and R² independently
matches the frozen sources. The cohort is identical to the previous
exploratory analysis. The previous 0.75/0.95 estimates reproduce.

HOW TO REPRODUCE IN VS CODE ON YOUR MAC
1. Extract the ZIP. Open its Yeast_NCLE_ED_NED_analysis folder in VS Code.
   Open Terminal > New Terminal. Run these terminal commands there:

   conda activate human_ncle_repro
   python --version
   python -m pip install -r requirements.txt

   Or create an environment with Python 3.12:
   conda create -n yeast_ncle_repro python=3.12 -y
   conda activate yeast_ncle_repro
   python -m pip install -r requirements.txt

2. Rebuild from frozen inputs into a new output directory:

   python scripts/run_yeast_r2.py \
     --turnover raw_sources/Christiano_yeast_mmc2.xlsx \
     --native raw_sources/AF_Yeast_entanglement.pkl \
     --knots raw_sources/knot_all.txt \
     --alphalasso raw_sources/AlphaLasso.zip \
     --keywords raw_sources/uniprot_S_cerevisiae_keyword_snapshot_2026_10_10.tsv \
     --primary-cutoff 0.90 --cutoffs 0.90 0.75 0.85 0.95 \
     --output-dir results/local_reproduction

   The output folder must be new or empty. No live API query is needed.

3. Independently check your rebuilt outputs:

   python scripts/verify_results.py \
     --results results/local_reproduction \
     --raw-sources raw_sources \
     --output verification/local_reproduction

   Expected: PASS: 88 numerical fields agree; cohort inclusion checks passed.

4. Open results/local_reproduction/association_results.csv in Excel.
   Compare the NCLE rows to the results above. Original shipped outputs
   remain in results/. A second rerun should use a different empty folder.
   To explore additional thresholds, change --cutoffs while retaining
   the chosen --primary-cutoff in that list, and report every result.

REFERENCE
Christiano et al., Cell Reports 9, 1959–1965 (2014).
https://doi.org/10.1016/j.celrep.2014.10.065

PACKAGING REVISION — 11 OCTOBER 2026
Removed an unused partial AlphaLasso download. The complete AlphaLasso.zip,
other scientific inputs, analysis scripts, result tables and workbook are
preserved. This cleanup does not change any eligibility rule or reported number.

INTERPRETATION CLARIFICATIONS
The R2-defined yeast endpoint is poor single-exponential fit. It differs from
McShane's NED classification in human and mouse and should be reported separately.
Within the eligible yeast cohort, R2 and slope CV have Spearman rho=-0.878226
(N=2,072); 21 of the 38 primary NED proteins with numeric positive rates are
in the fastest fifth of eligible numeric positive rates. These are descriptive
checks consistent with fit-quality concerns, not proof of a specific decay model.
The 100 grouped source records are correctly left unresolved. None has R2<0.90,
but adding ED-only observations can change the OR by changing ED denominators.
Their exclusion should not be described as mathematically unable to affect OR.
The primary upper confidence limit of 1.172 applies to this operational outcome
and filtered cohort, not to validated bi-exponential kinetics or misfolding.
The auditor's log(rate)-adjusted analysis is an additional sensitivity analysis,
separate from the requested length-only Model 1 reported in this package.
For cross-species methods, human membrane and secreted exports have filenames
dated 2025-03-19 and 2026-09-21, respectively. Mouse and yeast use the frozen
2026_03 UniProt snapshots. Filename/export dates are not release-version labels.

Revision verification: all 9 scientific CSVs and the workbook are unchanged.
Independent refit: 88 numerical comparisons passed; no failures.
