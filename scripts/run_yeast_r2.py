#!/usr/bin/env python3
"""Native-NCLE association with Quyen's R2-defined yeast ED/NED groups.

ED: R2 >= cutoff; NED: R2 < cutoff. The primary cutoff is 0.90.
These operational labels do not establish bi-exponential kinetics or misfolding.
Only crude association and Model 1 with standardized protein length are fitted.
"""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
import pickle
import platform
import re
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import scipy

from stats_model1 import crude, fit_model1


ALLOWED_GLOBALS = {
    ("pandas.core.frame", "DataFrame"),
    ("pandas.core.internals.managers", "BlockManager"),
    ("pandas._libs.internals", "_unpickle_block"),
    ("numpy.core.numeric", "_frombuffer"),
    ("numpy._core.numeric", "_frombuffer"),
    ("numpy", "dtype"), ("builtins", "slice"),
    ("numpy.core.multiarray", "_reconstruct"),
    ("numpy._core.multiarray", "_reconstruct"),
    ("numpy", "ndarray"), ("pandas.core.indexes.base", "_new_Index"),
    ("pandas.core.indexes.base", "Index"),
    ("pandas.core.indexes.range", "RangeIndex"),
}


class DataFrameUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        if (module, name) not in ALLOWED_GLOBALS:
            raise pickle.UnpicklingError(f"Unsupported global: {module}.{name}")
        return super().find_class(module, name)


def fingerprint(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return {"filename": path.name, "bytes": path.stat().st_size, "sha256": h.hexdigest()}


def read_native(path):
    with path.open("rb") as f:
        d = DataFrameUnpickler(f).load()
    required = {"Uniprot", "fragment", "uniprot_fragment", "length_AF",
                "plddt", "whGLN", "entanglement_AF"}
    if not isinstance(d, pd.DataFrame) or not required.issubset(d.columns):
        raise ValueError("Unexpected native annotation schema")
    d = d.copy()
    d["Uniprot"] = d.Uniprot.astype("string").str.strip()
    if d.Uniprot.isna().any() or not d.Uniprot.is_unique:
        raise ValueError("Native table requires one nonmissing row per UniProt accession")
    if d.Uniprot.str.contains(r"[;|]", regex=True).any():
        raise ValueError("Grouped native identifiers need a documented mapping")
    if not d.fragment.eq("F1").all():
        raise ValueError("Multiple/non-F1 fragments require a protein-level length/quality policy")
    if not d.uniprot_fragment.eq(d.Uniprot + "-F1").all():
        raise ValueError("Fragment identifier disagrees with UniProt accession")
    if not d.entanglement_AF.isin(["Yes", "No"]).all():
        raise ValueError("Native entanglement calls must be explicit Yes/No")
    for row in d.itertuples():
        if row.whGLN == "No":
            if row.entanglement_AF != "No":
                raise ValueError(f"Inconsistent native call for {row.Uniprot}")
        else:
            regions = ast.literal_eval(row.whGLN)
            if not isinstance(regions, list) or not regions or row.entanglement_AF != "Yes":
                raise ValueError(f"Inconsistent native region record for {row.Uniprot}")
    for column in ["length_AF", "plddt"]:
        d[column] = pd.to_numeric(d[column], errors="raise")
        if not np.isfinite(d[column]).all():
            raise ValueError(f"Nonfinite native field {column}")
    if d.length_AF.le(0).any() or not np.equal(d.length_AF, np.floor(d.length_AF)).all():
        raise ValueError("Invalid structural sequence lengths")
    if d.plddt.lt(0).any() or d.plddt.gt(100).any():
        raise ValueError("Native pLDDT outside [0,100]")
    d["NCLE"] = d.entanglement_AF.map({"Yes": 1, "No": 0})
    d["Protein_length"] = d.length_AF
    d["Mean_pLDDT"] = d.plddt
    return d


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--turnover", required=True, type=Path)
    parser.add_argument("--native", required=True, type=Path)
    parser.add_argument("--knots", required=True, type=Path)
    parser.add_argument("--alphalasso", required=True, type=Path)
    parser.add_argument("--keywords", required=True, type=Path)
    parser.add_argument("--primary-cutoff", type=float, default=.90)
    parser.add_argument("--cutoffs", type=float, nargs="+", default=[.90, .75, .85, .95])
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    if any(not 0 < cutoff < 1 for cutoff in args.cutoffs):
        raise ValueError("Each R2 cutoff must be strictly between 0 and 1")
    if len(set(args.cutoffs)) != len(args.cutoffs) or args.primary_cutoff not in args.cutoffs:
        raise ValueError("Cutoffs must be unique and include the primary cutoff")
    args.cutoffs = [args.primary_cutoff] + sorted(c for c in args.cutoffs if c != args.primary_cutoff)
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise ValueError("Choose a new or empty output directory")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    source = pd.read_excel(args.turnover, sheet_name=0)
    if source.UNIPROT.isna().any() or not source.UNIPROT.is_unique:
        raise ValueError("Turnover records must have unique, nonmissing accession groups")
    source["UniProt_ID"] = source.UNIPROT.astype("string").str.strip()
    source["Source_excel_row"] = range(2, len(source) + 2)
    source["Accession_count"] = source.UniProt_ID.str.split(";").map(len)
    source["Grouped_accession"] = source.Accession_count.gt(1)
    source["R2"] = pd.to_numeric(source["R2 (quality of curve fitting)"], errors="raise")
    if source.R2.isna().any() or source.R2.lt(0).any() or source.R2.gt(1).any():
        raise ValueError("Missing or invalid turnover R2 values")
    native = read_native(args.native)
    fields = ["Uniprot", "fragment", "uniprot_fragment", "NCLE", "Protein_length",
              "Mean_pLDDT", "entanglement_AF"]
    d = source.merge(native[fields], left_on="UniProt_ID", right_on="Uniprot",
                     how="left", validate="one_to_one")
    if len(d) != len(source):
        raise ValueError("Unexpected merge row expansion")

    knot_table = pd.read_csv(args.knots, sep=";", usecols=[0], dtype=str)
    knot_ids = set(knot_table.iloc[:, 0].dropna().str.strip())
    source_tokens = {token.strip() for value in source.UniProt_ID for token in value.split(";")}
    evidence = []
    with zipfile.ZipFile(args.alphalasso) as archive:
        for name in archive.namelist():
            if not name.endswith(".tsv") or "__MACOSX" in name:
                continue
            with archive.open(name) as f:
                for chunk in pd.read_csv(f, sep="\t", dtype=str, chunksize=100000):
                    # Fragment suffix normalization preserves meaningful isoform suffixes.
                    normalized = chunk.Uniprot.str.strip().str.replace(r"-F\d+$", "", regex=True)
                    mask = normalized.isin(source_tokens)
                    if mask.any():
                        found = chunk.loc[mask].copy()
                        found["UniProt_ID"] = normalized[mask]
                        found["AlphaLasso_table"] = name
                        evidence.append(found)
    evidence = pd.concat(evidence, ignore_index=True) if evidence else pd.DataFrame(columns=["UniProt_ID"])
    lasso_ids = set(evidence.UniProt_ID)
    d["Knot_member_IDs"] = [";".join(sorted(set(v.split(";")) & knot_ids)) for v in d.UniProt_ID]
    d["Covalent_lasso_member_IDs"] = [";".join(sorted(set(v.split(";")) & lasso_ids)) for v in d.UniProt_ID]
    d["Knotted"] = d.Knot_member_IDs.ne("").astype(int)
    d["Covalent_lasso"] = d.Covalent_lasso_member_IDs.ne("").astype(int)

    keywords = pd.read_csv(args.keywords, sep="\t", dtype=str)
    if not {"Entry", "Organism (ID)", "Length", "Keyword ID"}.issubset(keywords.columns):
        raise ValueError("Unexpected frozen keyword snapshot schema")
    if keywords.Entry.isna().any() or not keywords.Entry.is_unique:
        raise ValueError("Keyword snapshot has duplicate or missing accessions")
    keywords = keywords.rename(columns={"Entry": "UniProt_ID", "Length": "Current_UniProt_length",
                                        "Organism (ID)": "Current_taxonomy_ID", "Organism": "Current_organism",
                                        "Keyword ID": "UniProt_keyword_IDs", "Keywords": "UniProt_keywords"})
    d = d.merge(keywords, on="UniProt_ID", how="left", validate="one_to_one")
    d["Keyword_source_covered"] = d.Current_taxonomy_ID.notna()
    kw = d.UniProt_keyword_IDs.fillna("").str.split(";").map(lambda v: {s.strip() for s in v})
    for name, keyword in [("Membrane_keyword", "KW-0472"), ("Secreted_keyword", "KW-0964")]:
        d[name] = pd.array([int(keyword in x) if covered else pd.NA
                            for x, covered in zip(kw, d.Keyword_source_covered)], dtype="Int64")
    d["Current_UniProt_length"] = pd.to_numeric(d.Current_UniProt_length, errors="coerce")
    d["Length_differs_from_current_UniProt"] = (
        d.Protein_length.notna() & d.Current_UniProt_length.notna()
        & d.Protein_length.ne(d.Current_UniProt_length)).astype(int)

    rules = [
        ("Retain unambiguous single-accession turnover records", ~d.Grouped_accession),
        ("Require exact native-annotation accession match", d.NCLE.notna()),
        ("Retain supplied mean pLDDT >=70", d.Mean_pLDDT.ge(70)),
        ("Exclude supplied knotted proteins", d.Knotted.eq(0)),
        ("Exclude supplied covalent-lasso proteins", d.Covalent_lasso.eq(0)),
        ("Require coverage in frozen UniProt keyword snapshot", d.Keyword_source_covered),
        ("Exclude UniProt Membrane keyword KW-0472", d.Membrane_keyword.eq(0)),
        ("Exclude UniProt Secreted keyword KW-0964", d.Secreted_keyword.eq(0)),
    ]
    flow, retained = [], pd.Series(True, index=d.index)
    d["Disposition"] = "Included"

    def capture(step, removed):
        sub = d.loc[retained]
        flow.append({"Step": step, "Remaining": len(sub), "Removed_at_step": removed,
                     "NCLE_1": int(sub.NCLE.eq(1).sum()), "NCLE_0": int(sub.NCLE.eq(0).sum()),
                     "NCLE_missing": int(sub.NCLE.isna().sum()),
                     **{f"NED_{str(c).replace('.', '_')}": int(sub.R2.lt(c).sum()) for c in args.cutoffs}})

    capture("All turnover source records", 0)
    for step, keep in rules:
        keep = keep.fillna(False)
        removed = retained & ~keep
        d.loc[removed, "Disposition"] = step
        retained &= keep
        capture(step, int(removed.sum()))
    d["Included"] = retained.astype(int)
    for threshold in args.cutoffs:
        d[f"NED_{str(threshold).replace('.', '_')}"] = d.R2.lt(threshold).astype(int)
    d["NED"] = d.R2.lt(args.primary_cutoff).astype(int)
    d["Degradation_profile"] = np.where(d.NED.eq(1), "NED", "ED")
    cohort = d.loc[retained].copy()
    if not cohort.UniProt_ID.is_unique or cohort.NCLE.isna().any():
        raise ValueError("Invalid final cohort accession/NCLE integrity")

    all_results, counts, diagnostics = [], [], []
    for threshold in args.cutoffs:
        sub = cohort.copy()
        sub["NED"] = sub.R2.lt(threshold).astype(int)
        outcome = f"NED (R2 < {threshold:.2f}) vs ED (R2 >= {threshold:.2f})"
        analysis = "Primary" if threshold == args.primary_cutoff else "Cutoff sensitivity"
        raw, table_rows = crude(sub)
        all_results.append({"Analysis": analysis, "Outcome": outcome, "R2_cutoff": threshold, "Method": "Crude", "Term": "NCLE", **raw})
        for row in table_rows:
            counts.append({"Analysis": analysis, "Outcome": outcome, "R2_cutoff": threshold,
                           **row, "NED_fraction": row["NED"] / row["Total"]})
        rows, diagnostic = fit_model1(sub)
        if diagnostic["Status"] != "ok":
            raise ValueError(f"Model could not be estimated at R2={threshold}: {diagnostic}")
        all_results.extend({"Analysis": analysis, "Outcome": outcome, "R2_cutoff": threshold, "Method": "Model 1", **row} for row in rows)
        diagnostics.append({"Analysis": analysis, "R2_cutoff": threshold, "Events": int(sub.NED.sum()),
                            "Non_events": int(sub.NED.eq(0).sum()),
                            "Exactly_at_cutoff": int(sub.R2.eq(threshold).sum()),
                            "Smallest_2x2_cell": min(int(sub.loc[sub.NCLE.eq(n), "NED"].eq(y).sum())
                                                     for n in [0, 1] for y in [0, 1]),
                            **diagnostic})
    cohort["z_Protein_length"] = (cohort.Protein_length - cohort.Protein_length.mean()) / cohort.Protein_length.std(ddof=1)

    native.to_csv(args.output_dir / "native_annotation_export.csv", index=False)
    d.to_csv(args.output_dir / "all_turnover_record_audit.csv", index=False)
    d.loc[d.Grouped_accession].to_csv(args.output_dir / "grouped_accession_review.csv", index=False)
    cohort.to_csv(args.output_dir / "eligible_yeast_cohort.csv", index=False)
    evidence.to_csv(args.output_dir / "covalent_lasso_evidence.csv", index=False)
    pd.DataFrame(flow).to_csv(args.output_dir / "cohort_flow.csv", index=False)
    pd.DataFrame(all_results).to_csv(args.output_dir / "association_results.csv", index=False)
    pd.DataFrame(counts).to_csv(args.output_dir / "contingency_counts.csv", index=False)
    pd.DataFrame(diagnostics).to_csv(args.output_dir / "model_diagnostics.csv", index=False)
    manifest = {
        "Status": "Quyen-approved operational R2-based ED/NED classification",
        "Species": "Saccharomyces cerevisiae", "Final_N": len(cohort),
        "Model": "NED ~ NCLE + z(Protein_length)",
        "Primary_R2_cutoff": args.primary_cutoff,
        "R2_cutoffs": args.cutoffs,
        "Classification": "NED=1 if R2 < cutoff; ED=0 if R2 >= cutoff; exactly the cutoff is ED",
        "Interpretation": "Poor single-exponential fit is an operational NED label, not confirmation of bi-exponential kinetics or misfolding",
        "Length_source": "Supplied length_AF, preserved even if current UniProt length differs",
        "Quality_source": "Supplied plddt interpreted as protein mean; values not recomputed",
        "Length_z_definition": "Within eligible cohort, sample SD ddof=1",
        "NCLE_coding": "Supplied entanglement_AF Yes=1, No=0; separate knot/covalent exclusions",
        "Membrane_definition": "Presence of UniProt keyword KW-0472",
        "Secreted_definition": "Presence of UniProt keyword KW-0964",
        "Keywords_source": "Frozen live UniProt taxonomy_id:4932 snapshot; all retained IDs covered",
        "No_disorder_filter": True, "No_complex_covariate": True,
        "R2_outcomes": [f"R2<{c:.2f}" for c in args.cutoffs],
        "Multiple_testing_adjustment": "None; primary 0.90 and all three cutoff-sensitivity results reported without selecting by P-value",
        "Native_source_rows": len(native), "Native_source_NCLE_counts": native.NCLE.value_counts().to_dict(),
        "Native_pLDDT_ge70": int(native.Mean_pLDDT.ge(70).sum()),
        "Current_UniProt_length_mismatches_in_cohort": int(cohort.Length_differs_from_current_UniProt.sum()),
        "Nonnumeric_corrected_rates_in_cohort": int(pd.to_numeric(cohort["Degradation rates (min-1)"], errors="coerce").isna().sum()),
        "Input_sources": [fingerprint(p) for p in [args.turnover, args.native, args.knots, args.alphalasso, args.keywords]],
        "Runtime": {"Python": platform.python_version(), "pandas": pd.__version__, "numpy": np.__version__, "scipy": scipy.__version__},
        "Model_diagnostics": diagnostics,
    }
    (args.output_dir / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2, default=str) + "\n")
    # Compact typed payload for the workbook builder; preserve full audit in CSV.
    display = ["UniProt_ID", "NCLE", "Mean_pLDDT", "Protein_length", "R2", "Degradation_profile", "NED",
               *[f"NED_{str(c).replace('.', '_')}" for c in sorted(args.cutoffs)],
               "Knotted", "Covalent_lasso", "Membrane_keyword", "Secreted_keyword", "Keyword_source_covered",
               "Grouped_accession", "Included", "Disposition", "Source_excel_row", "Current_UniProt_length",
               "Length_differs_from_current_UniProt", "Gene name", "Protein name"]
    payload = {"results": all_results, "counts": counts, "flow": flow,
               "cohort": json.loads(cohort[display + ["z_Protein_length"]].to_json(orient="records")),
               "audit": json.loads(d[display].to_json(orient="records")), "manifest": manifest}
    (args.output_dir / "workbook_payload.json").write_text(json.dumps(payload, indent=2, default=str) + "\n")
    print(pd.DataFrame(flow).to_string(index=False))
    print(pd.DataFrame(all_results).loc[lambda x: x.Term.eq("NCLE")].to_string(index=False))


if __name__ == "__main__":
    main()
