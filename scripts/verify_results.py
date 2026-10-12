#!/usr/bin/env python3
"""Independent numerical refit and direct source/cohort checks for the yeast report."""
from pathlib import Path
import argparse
import json
import zipfile
import numpy as np
import pandas as pd
from scipy.special import expit
from scipy.stats import fisher_exact, norm


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--results", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--raw-sources", type=Path,
                   default=Path(__file__).resolve().parents[1] / "raw_sources",
                   help="Frozen source directory; defaults to the package raw_sources folder")
    args = p.parse_args()
    data = pd.read_csv(args.results / "eligible_yeast_cohort.csv")
    reported = pd.read_csv(args.results / "association_results.csv")
    audit = pd.read_csv(args.results / "all_turnover_record_audit.csv")
    assert len(data) == len(audit.loc[audit.Included.eq(1)])
    assert set(data.UniProt_ID) == set(audit.loc[audit.Included.eq(1), "UniProt_ID"])
    assert data.UniProt_ID.is_unique
    assert data.Mean_pLDDT.ge(70).all()
    assert data.NCLE.isin([0, 1]).all()
    assert data[["Knotted", "Covalent_lasso", "Membrane_keyword", "Secreted_keyword"]].eq(0).all().all()
    raw = args.raw_sources
    source = pd.read_excel(raw / "Christiano_yeast_mmc2.xlsx", sheet_name=0)
    native = pd.read_pickle(raw / "AF_Yeast_entanglement.pkl")
    candidate = source.loc[~source.UNIPROT.str.contains(";", regex=False)].merge(
        native, left_on="UNIPROT", right_on="Uniprot", how="inner", validate="one_to_one")
    candidate = candidate.loc[candidate.plddt.ge(70)].copy()
    knot_ids = set(pd.read_csv(raw / "knot_all.txt", sep=";", usecols=[0], dtype=str).iloc[:, 0])
    lasso_ids = set()
    with zipfile.ZipFile(raw / "AlphaLasso.zip") as archive:
        for name in archive.namelist():
            if not name.endswith(".tsv") or "__MACOSX" in name:
                continue
            with archive.open(name) as handle:
                for chunk in pd.read_csv(handle, sep="\t", usecols=["Uniprot"], dtype=str, chunksize=100000):
                    lasso_ids.update(chunk.Uniprot.str.replace(r"-F\d+$", "", regex=True))
    candidate = candidate.loc[~candidate.Uniprot.isin(knot_ids | lasso_ids)]
    keyword = pd.read_csv(raw / "uniprot_S_cerevisiae_keyword_snapshot_2026_10_10.tsv", sep="\t", dtype=str)
    candidate = candidate.merge(keyword, left_on="Uniprot", right_on="Entry", how="inner", validate="one_to_one")
    kw = candidate["Keyword ID"].fillna("").str.split(";").map(lambda xs: {v.strip() for v in xs})
    candidate = candidate.loc[~kw.map(lambda v: "KW-0472" in v or "KW-0964" in v)]
    assert set(candidate.Uniprot) == set(data.UniProt_ID)
    provenance = candidate.merge(data, left_on="Uniprot", right_on="UniProt_ID", validate="one_to_one",
                                 suffixes=("_source", "_result"))
    assert np.allclose(provenance.plddt, provenance.Mean_pLDDT, rtol=1e-12)
    assert np.array_equal(provenance.length_AF, provenance.Protein_length)
    assert np.array_equal(provenance.entanglement_AF_source.map({"Yes": 1, "No": 0}), provenance.NCLE)
    assert np.allclose(provenance["R2 (quality of curve fitting)_source"], provenance.R2, rtol=1e-12)
    comparisons = []
    manifest = json.loads((args.results / "analysis_manifest.json").read_text())
    primary = manifest["Primary_R2_cutoff"]
    assert np.array_equal(data.NED, data.R2.lt(primary).astype(int))
    assert np.array_equal(data.Degradation_profile, np.where(data.R2.lt(primary), "NED", "ED"))
    for cutoff in manifest["R2_cutoffs"]:
        assert np.array_equal(data[f"NED_{str(cutoff).replace(chr(46), chr(95))}"], data.R2.lt(cutoff).astype(int))
        y = data.R2.lt(cutoff).to_numpy(dtype=float)
        ncle = data.NCLE.to_numpy(dtype=float)
        length = data.Protein_length.to_numpy(dtype=float)
        standardized = (length - np.mean(length)) / np.std(length, ddof=1)
        X = np.column_stack([np.ones(len(y)), ncle, standardized])
        beta = np.zeros(3)
        # Independent Fisher scoring / IRLS, without importing production routines.
        for iteration in range(100):
            probabilities = expit(X @ beta)
            information = X.T @ ((probabilities * (1 - probabilities))[:, None] * X)
            step = np.linalg.solve(information, X.T @ (y - probabilities))
            beta += step
            if np.max(np.abs(step)) < 1e-11:
                break
        else:
            raise ValueError("Independent IRLS did not converge")
        probabilities = expit(X @ beta)
        information = X.T @ ((probabilities * (1 - probabilities))[:, None] * X)
        se = np.sqrt(np.diag(np.linalg.inv(information)))
        for i, term in enumerate(["Intercept", "NCLE", "z_Protein_length"]):
            values = {"Coefficient": beta[i], "SE": se[i], "OR": np.exp(beta[i]),
                      "CI_lower": np.exp(beta[i] - norm.ppf(.975) * se[i]),
                      "CI_upper": np.exp(beta[i] + norm.ppf(.975) * se[i]),
                      "P_value": 2 * norm.sf(abs(beta[i] / se[i]))}
            row = reported.loc[reported.R2_cutoff.eq(cutoff) & reported.Method.eq("Model 1")
                               & reported.Term.eq(term)].iloc[0]
            for field, value in values.items():
                passed = bool(np.isclose(row[field], value, rtol=1e-6, atol=1e-8))
                comparisons.append({"R2_cutoff": cutoff, "Method": "Model 1", "Term": term,
                                    "Field": field, "Reported": float(row[field]),
                                    "Independent": float(value), "Pass": passed})
        table = np.array([[int(((ncle == n) & (y == status)).sum())
                           for status in [1, 0]] for n in [1, 0]])
        a, b, c, d = table.flatten()
        crude_or = (a * d) / (b * c)
        log_se = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
        values = {"OR": crude_or,
                  "CI_lower": np.exp(np.log(crude_or) - norm.ppf(.975) * log_se),
                  "CI_upper": np.exp(np.log(crude_or) + norm.ppf(.975) * log_se),
                  "P_value": fisher_exact(table, alternative="two-sided").pvalue}
        row = reported.loc[reported.R2_cutoff.eq(cutoff) & reported.Method.eq("Crude")].iloc[0]
        for field, value in values.items():
            passed = bool(np.isclose(row[field], value, rtol=1e-10, atol=1e-12))
            comparisons.append({"R2_cutoff": cutoff, "Method": "Crude", "Term": "NCLE",
                                "Field": field, "Reported": float(row[field]),
                                "Independent": float(value), "Pass": passed})
    checks = pd.DataFrame(comparisons)
    args.output.mkdir(parents=True, exist_ok=True)
    checks.to_csv(args.output / "independent_numerical_checks.csv", index=False)
    summary = {"Numerical_fields_compared": len(checks), "All_pass": bool(checks.Pass.all()),
               "Independent_refit": "Fisher scoring / IRLS implemented without production functions",
               "Regression_rtol": 1e-6, "Regression_atol": 1e-8,
               "Cohort_integrity_and_inclusion_checks": True,
               "Independent_raw_source_reconstruction_exact_ID_match": True,
               "Proteins_compared_to_native_and_turnover_sources": len(provenance),
               "Per_protein_NCLE_length_pLDDT_and_R2_match": True}
    (args.output / "verification_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    if not checks.Pass.all():
        print(checks.loc[~checks.Pass].to_string(index=False))
        raise ValueError("Independent numerical verification failed")
    print(f"PASS: {len(checks)} numerical fields agree; cohort inclusion checks passed.")


if __name__ == "__main__":
    main()
