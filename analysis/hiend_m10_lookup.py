# -*- coding: utf-8 -*-
"""M10 Geneformer lookup and calls (analysis plan fixed before running)
Positive = >= 1 signature gene with Goal_end_FDR < 0.05 (template criterion); genes not tested count as non-positive.
Sensitivity: signatures with test coverage >= 50% only. Real vs random: two-sided Fisher; random positive rate with Clopper-Pearson 95% CI, lower limit > 0.075 -> excess type I error.
Signatures = N11_out/hiend_m8/signature_genes_all_arms.tsv (same set as proxy analysis / M8 / M3 / M12).
"""
import os
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact, beta

B = os.environ.get("HCC_AUDIT_ROOT", ".")
O = os.path.join(B, "N11_out", "hiend_m10")
st = pd.read_csv(os.path.join(O, "gf", "stats_delete_all_n1000", "all_genes_delete_all_n1000_stats.csv"))
tested = set(st.loc[st["Goal_end_FDR"].notna(), "Gene_name"].astype(str))
pos_genes = set(st.loc[st["Goal_end_FDR"] < 0.05, "Gene_name"].astype(str))

# Note: the literal "null" in column kind is read as missing by pandas by default -> groupby silently drops all random signatures; keep_default_na=False is required
G = pd.read_csv(os.path.join(B, "N11_out", "hiend_m8", "signature_genes_all_arms.tsv"), sep="\t", keep_default_na=False, na_values=[""])
assert G["kind"].isin(["real", "null"]).all() and G["kind"].eq("null").any(), G["kind"].value_counts(dropna=False)
G["sid"] = G["arm"] + "|" + G["src"] + "|" + G["kind"] + "|" + G["rep"].astype(str)
S = G.groupby(["sid", "arm", "src", "kind", "rep"]).agg(
    k=("gene", "nunique"),
    n_tested=("gene", lambda x: len(set(x) & tested)),
    n_pos=("gene", lambda x: len(set(x) & pos_genes))).reset_index()
S["coverage"] = S["n_tested"] / S["k"]
S["M10_pos"] = (S["n_pos"] >= 1).astype(int)
assert len(S) == G["sid"].nunique(), (len(S), G["sid"].nunique())
S.to_csv(os.path.join(O, "m10_per_signature.tsv"), sep="\t", index=False)

full_ids = set(pd.read_csv(os.path.join(B, "N11_out", "ring3_full", "published_coverage.tsv"), sep="\t").query("full == True")["sig_id"])


def cp(x, n):
    lo = 0.0 if x == 0 else beta.ppf(0.025, x, n - x + 1)
    hi = 1.0 if x == n else beta.ppf(0.975, x + 1, n - x)
    return lo, hi


def comp(re, nu, lab, filt):
    xr, nr, xn, nn = int(re["M10_pos"].sum()), len(re), int(nu["M10_pos"].sum()), len(nu)
    lo, hi = cp(xn, nn)
    orr, p = fisher_exact([[xr, nr - xr], [xn, nn - xn]]) if nr and nn else (np.nan, np.nan)
    return dict(comparison=lab, filter=filt, n_real=nr, rate_real=xr / nr if nr else np.nan, n_null=nn, FPR_null=xn / nn if nn else np.nan,
                FPR_CP_lo=lo, FPR_CP_hi=hi, typeI_inflated=lo > 0.075, fisher_OR=orr, fisher_P=p,
                coverage_median_real=re["coverage"].median(), coverage_median_null=nu["coverage"].median(),
                k_median_real=re["k"].median(), k_median_null=nu["k"].median())


rows = []
for filt in ["all", "coverage>=0.5"]:
    T = S if filt == "all" else S[S["coverage"] >= 0.5]
    A_r = T[(T.arm == "A") & (T.kind == "real")]
    rows.append(comp(A_r, T[(T.arm == "A") & (T.kind == "null")], "A_vs_fitsig_k", filt))
    rows.append(comp(A_r, T[T.arm == "Ak"], "A_vs_Ak_sizematched", filt))
    rows.append(comp(T[(T.arm == "B") & (T.kind == "real")], T[(T.arm == "B") & (T.kind == "null")], "B_vs_fitsig_k", filt))
    rows.append(comp(A_r[A_r.src.isin(full_ids)], T[(T.arm == "Ak") & T.src.isin(full_ids)], "Afull_vs_Ak_sizematched", filt))
D = pd.DataFrame(rows)
D.to_csv(os.path.join(O, "m10_decision.tsv"), sep="\t", index=False)
gl = pd.DataFrame([dict(genes_with_result=len(tested), genes_FDR05=len(pos_genes), frac=len(pos_genes) / len(tested),
                        FDR05_shift_pos=int(((st["Goal_end_FDR"] < 0.05) & (st["Shift_to_goal_end"] > 0)).sum()),
                        FDR05_shift_neg=int(((st["Goal_end_FDR"] < 0.05) & (st["Shift_to_goal_end"] < 0)).sum()),
                        N_Detections_median=float(st["N_Detections"].median()))])
gl.to_csv(os.path.join(O, "m10_gene_level.tsv"), sep="\t", index=False)
pd.set_option("display.width", 250)
print(gl.to_string(index=False)); print(D.to_string(index=False)); print("DONE")
