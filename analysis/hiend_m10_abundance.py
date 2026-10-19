# -*- coding: utf-8 -*-
"""M10 Geneformer: detection-abundance bias check (analysis plan fixed before running)"""
import os
import gzip
import random
import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.stats import mannwhitneyu, spearmanr, fisher_exact, beta

B = os.environ.get("HCC_AUDIT_ROOT", ".")
GF = os.path.join(B, "N11_out", "hiend_m10", "gf"); O = os.path.join(B, "N11_out", "hiend_m10")
st = pd.read_csv(os.path.join(GF, "stats_delete_all_n1000", "all_genes_delete_all_n1000_stats.csv"))
st = st[st["Goal_end_FDR"].notna()].copy()
st["sig"] = st["Goal_end_FDR"] < 0.05
st = st.drop_duplicates("Gene_name")

# Auxiliary abundance: all Disease cells
feats = pd.read_csv(os.path.join(GF, "features.tsv"), sep="\t")["gene_symbol"].astype(str).values
meta = pd.read_csv(os.path.join(GF, "cell_meta.tsv"), sep="\t", index_col=0)
with gzip.open(os.path.join(GF, "counts.mtx.gz"), "rb") as fh:
    m = mmread(fh).tocsr()                       # genes x cells
dis = np.where(meta["Type"].values == "Disease")[0]
md = m[:, dis]; lib = np.asarray(md.sum(axis=0)).ravel()
cpm = md.multiply(1e6 / lib).tocsr(); cpm.data = np.log1p(cpm.data)
ab = pd.DataFrame({"Gene_name": feats, "mean_log1p_cpm": np.asarray(cpm.mean(axis=1)).ravel(),
                   "detect_rate": np.asarray((md > 0).mean(axis=1)).ravel()}).drop_duplicates("Gene_name")
st = st.merge(ab, on="Gene_name", how="left")

# ---- Check 1 ----
rows = []
for col in ["N_Detections", "mean_log1p_cpm", "detect_rate"]:
    a = st.loc[st.sig, col].dropna(); b = st.loc[~st.sig, col].dropna()
    rows.append(dict(measure=col, n_sig=len(a), n_nonsig=len(b), median_sig=a.median(), median_nonsig=b.median(),
                     MWU_P=mannwhitneyu(a, b, alternative="two-sided").pvalue))
r1 = pd.DataFrame(rows)
rho, prho = spearmanr(-np.log10(st["Goal_end_FDR"].clip(lower=1e-300)), st["N_Detections"])
r1["spearman_rho_neglog10FDR_vs_NDet"] = rho; r1["spearman_P"] = prho
r1.to_csv(os.path.join(O, "m10_abund_check1.tsv"), sep="\t", index=False)
st["NDet_decile"] = pd.qcut(st["N_Detections"].rank(method="first"), 10, labels=False) + 1
dec = st.groupby("NDet_decile").agg(n=("sig", "size"), n_sig=("sig", "sum"), frac_sig=("sig", "mean"),
                                    NDet_min=("N_Detections", "min"), NDet_max=("N_Detections", "max")).reset_index()
dec.to_csv(os.path.join(O, "m10_abund_deciles.tsv"), sep="\t", index=False)
st.to_csv(os.path.join(O, "m10_abund_gene_table.tsv"), sep="\t", index=False)

# ---- Check 2 ----
G = pd.read_csv(os.path.join(B, "N11_out", "hiend_m8", "signature_genes_all_arms.tsv"), sep="\t", keep_default_na=False, na_values=[""])
assert G["kind"].isin(["real", "null"]).all() and G["kind"].eq("null").any()
U = [g for g in open(os.path.join(B, "N11_out", "hiend_m8", "U_pool.txt")).read().split()]
dec_of = dict(zip(st["Gene_name"], st["NDet_decile"])); sigset = set(st.loc[st.sig, "Gene_name"])
pools = {d: sorted(st.loc[st.NDet_decile == d, "Gene_name"]) for d in range(1, 11)}
untested = sorted(set(U) - set(dec_of))
pub = G[(G.arm == "A") & (G.kind == "real")].groupby("src")["gene"].apply(lambda x: sorted(set(x)))
res = []
for i, (sid, genes) in enumerate(pub.items(), start=1):
    real_pos = int(any(g in sigset for g in genes))
    layers = [dec_of.get(g, 0) for g in genes]
    res.append(dict(src=sid, type="real", rep=0, pos=real_pos, k=len(genes), n_tested=sum(l > 0 for l in layers)))
    for rep in range(1, 201):
        rng = random.Random(int(10e6 + i * 1e3 + rep)); draw = []
        for l in sorted(set(layers)):
            n = layers.count(l); pool = pools[l] if l > 0 else untested
            draw += rng.sample(pool, min(n, len(pool)))
        res.append(dict(src=sid, type="abund_matched", rep=rep, pos=int(any(g in sigset for g in draw)), k=len(draw),
                        n_tested=sum(g in dec_of for g in draw)))
R = pd.DataFrame(res); R.to_csv(os.path.join(O, "m10_abund_matched_rows.tsv"), sep="\t", index=False)
re_, nu = R[R.type == "real"], R[R.type == "abund_matched"]
xr, nr, xn, nn = int(re_.pos.sum()), len(re_), int(nu.pos.sum()), len(nu)
lo = 0.0 if xn == 0 else beta.ppf(0.025, xn, nn - xn + 1); hi = 1.0 if xn == nn else beta.ppf(0.975, xn + 1, nn - xn)
orr, p = fisher_exact([[xr, nr - xr], [xn, nn - xn]])
c1 = r1.set_index("measure").loc["N_Detections"]
verdict = ("explained_by_abundance" if (c1.MWU_P < 0.05 and c1.median_sig > c1.median_nonsig and p >= 0.05)
           else "difference_remains" if p < 0.05 else "abundance_claim_withdrawn")
r2 = pd.DataFrame([dict(n_real=nr, rate_real=xr / nr, n_null=nn, rate_null=xn / nn, null_CP_lo=lo, null_CP_hi=hi, fisher_OR=orr,
                        fisher_P=p, verdict=verdict)])
r2.to_csv(os.path.join(O, "m10_abund_check2.tsv"), sep="\t", index=False)
pd.set_option("display.width", 250)
print(r1.to_string(index=False)); print(dec.to_string(index=False)); print(r2.to_string(index=False)); print("DONE")
