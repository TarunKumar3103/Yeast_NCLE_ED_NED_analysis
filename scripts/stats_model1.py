"""Numerical routines for crude and length-adjusted R2-defined NED association.

Adapted from the independently validated earlier Model 1 runner.
Quyen-defined NED=1 is R2 below the cutoff; ED=0 is R2 at or above it.
"""
from __future__ import annotations
import numpy as np
import pandas as pd
from scipy import optimize, special, stats

def crude(cohort: pd.DataFrame) -> tuple[dict, list[dict]]:
    a = int(((cohort["NCLE"] == 1) & (cohort["NED"] == 1)).sum())
    b = int(((cohort["NCLE"] == 1) & (cohort["NED"] == 0)).sum())
    c = int(((cohort["NCLE"] == 0) & (cohort["NED"] == 1)).sum())
    d = int(((cohort["NCLE"] == 0) & (cohort["NED"] == 0)).sum())
    table = np.array([[a, b], [c, d]], dtype=int)
    rows = [
        {"NCLE": 1, "NED": a, "ED": b, "Total": a + b},
        {"NCLE": 0, "NED": c, "ED": d, "Total": c + d},
    ]
    if len(cohort) == 0 or (table.sum(axis=0) == 0).any() or (table.sum(axis=1) == 0).any():
        return {"N": len(cohort), "Status": "not_estimable_missing_comparison_group"}, rows
    estimate, p = stats.fisher_exact(table, alternative="two-sided")
    if np.all(table > 0):
        se = np.sqrt(np.sum(1.0 / table))
        lo, hi = np.exp(np.log(estimate) + stats.norm.ppf(0.975) * np.array([-se, se]))
        ci_method = "Woolf log-OR Wald"
    else:
        interval = stats.contingency.odds_ratio(table, kind="conditional").confidence_interval(0.95)
        lo, hi = interval.low, interval.high
        ci_method = "Exact conditional interval (uncorrected sample OR reported)"
    return {
        "N": len(cohort), "OR": float(estimate), "CI_lower": float(lo), "CI_upper": float(hi),
        "P_value": float(p), "P_method": "Two-sided Fisher exact", "CI_method": ci_method,
        "Status": "ok" if np.all(table > 0) else "zero_cell",
    }, rows


def fit_model1(cohort: pd.DataFrame) -> tuple[list[dict], dict]:
    n = len(cohort)
    if n == 0 or cohort["NCLE"].nunique() != 2 or cohort["NED"].nunique() != 2:
        return [], {"N": n, "Status": "not_estimable_missing_comparison_group"}
    mean = float(cohort["Protein_length"].mean())
    sd = float(cohort["Protein_length"].std(ddof=1))
    if not np.isfinite(sd) or sd <= 0:
        return [], {"N": n, "Status": "not_estimable_constant_length"}
    cohort["z_Protein_length"] = (cohort["Protein_length"] - mean) / sd
    x = np.column_stack([np.ones(n), cohort["NCLE"], cohort["z_Protein_length"]])
    y = cohort["NED"].to_numpy(dtype=float)
    if np.linalg.matrix_rank(x) != 3:
        return [], {"N": n, "Status": "not_estimable_rank_deficient"}
    signed = (2 * y - 1)[:, None] * x
    separation = optimize.linprog(
        np.zeros(3), A_ub=np.vstack([-signed, -signed.sum(axis=0)]),
        b_ub=np.r_[np.zeros(n), -1.0], bounds=[(None, None)] * 3, method="highs",
    )
    if separation.status == 0:
        return [], {"N": n, "Status": "not_estimable_complete_or_quasi_separation"}
    if separation.status != 2:
        return [], {"N": n, "Status": "separation_diagnostic_failed", "Message": separation.message}

    def loss(beta):
        eta = x @ beta
        return float(np.sum(np.logaddexp(0, eta) - y * eta))

    def gradient(beta):
        return x.T @ (special.expit(x @ beta) - y)

    def information(beta):
        p = special.expit(x @ beta)
        return x.T @ ((p * (1 - p))[:, None] * x)

    initial = np.array([special.logit(y.mean()), 0.0, 0.0])
    fitted = optimize.minimize(loss, initial, method="trust-exact", jac=gradient, hess=information,
                               options={"gtol": 1e-8, "maxiter": 1000})
    score_max = float(np.max(np.abs(gradient(fitted.x))))
    if not np.isfinite(fitted.fun) or score_max > 1e-6 or np.max(np.abs(fitted.x)) > 30:
        return [], {"N": n, "Status": "logistic_fit_not_reliable", "Message": str(fitted.message),
                    "Max_abs_score": score_max}
    info = information(fitted.x)
    covariance = np.linalg.inv(info)
    errors = np.sqrt(np.diag(covariance))
    if not np.isfinite(errors).all():
        return [], {"N": n, "Status": "nonfinite_standard_errors"}
    rows = []
    for term, beta, se in zip(["Intercept", "NCLE", "z_Protein_length"], fitted.x, errors):
        z = beta / se
        radius = stats.norm.ppf(0.975) * se
        rows.append({"Term": term, "N": n, "Coefficient": float(beta), "SE": float(se),
                     "OR": float(np.exp(beta)), "CI_lower": float(np.exp(beta - radius)),
                     "CI_upper": float(np.exp(beta + radius)), "P_value": float(2 * stats.norm.sf(abs(z))),
                     "P_method": "Two-sided Wald z", "CI_method": "Logistic Wald", "Status": "ok"})
    return rows, {"N": n, "Status": "ok", "Mean_length": mean, "SD_length_ddof1": sd,
                  "Log_likelihood": -float(fitted.fun), "Max_abs_score": score_max,
                  "Information_condition_number": float(np.linalg.cond(info)), "Iterations": int(fitted.nit)}


