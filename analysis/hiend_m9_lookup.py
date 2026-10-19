# -*- coding: utf-8 -*-
"""M9 scTenifoldKnk 查表判读（预注册 N20_m9_knk_prereg.md，结果出来后不改）
基因级：去掉靶基因后 diffRegulation p.adj<0.05 的基因数 >=1 → 该基因模板阳性
签名级：签名中网络内基因 >=1 个模板阳性 → 签名阳性；网络内无基因的签名按非阳性计
敏感性：只算网络内至少 1 个基因的签名。真实 vs 随机 Fisher 双侧；随机阳性率 Clopper–Pearson 95% CI，下界 > 0.075 → I 类错误失控
"""
import os
import glob
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact, beta

B = os.environ.get("HCC_AUDIT_ROOT", ".")
O = os.path.join(B, "N11_out", "hiend_m9_hvg")
KD = os.path.join(O, "ko_formal")
kol = pd.read_csv(os.path.join(O, "ko_list_formal.tsv"), sep="\t")
rows = []
for g in kol["gene"]:
    f = os.path.join(KD, f"{g}.tsv")
    if not os.path.exists(f):
        rows.append(dict(gene=g, done=False, n_tested=np.nan, n_sig=np.nan)); continue
    d = pd.read_csv(f, sep="\t", keep_default_na=False, na_values=["", "NA"])
    d = d[d["gene"] != g]
    rows.append(dict(gene=g, done=True, n_tested=len(d), n_sig=int((d["p.adj"] < 0.05).sum())))
GL = pd.DataFrame(rows)
GL["gene_pos"] = (GL["n_sig"] >= 1).astype(int)
GL.to_csv(os.path.join(O, "m9_gene_level.tsv"), sep="\t", index=False)
n_missing = int((~GL["done"]).sum())
assert n_missing == 0, f"{n_missing} KO genes not finished"
pos = set(GL.loc[GL.gene_pos == 1, "gene"]); W = set(open(os.path.join(O, "genes_after_qc.txt")).read().split())

G = pd.read_csv(os.path.join(B, "N11_out", "hiend_m8", "signature_genes_all_arms.tsv"), sep="\t", keep_default_na=False, na_values=[""])
assert G["kind"].isin(["real", "null"]).all() and G["kind"].eq("null").any()
G["sid"] = G["arm"] + "|" + G["src"] + "|" + G["kind"] + "|" + G["rep"].astype(str)
S = G.groupby(["sid", "arm", "src", "kind", "rep"]).agg(
    k=("gene", "nunique"), n_in_net=("gene", lambda x: len(set(x) & W)), n_pos=("gene", lambda x: len(set(x) & pos))).reset_index()
assert len(S) == G["sid"].nunique()
S["coverage"] = S["n_in_net"] / S["k"]
S["M9_pos"] = (S["n_pos"] >= 1).astype(int)
S.to_csv(os.path.join(O, "m9_per_signature.tsv"), sep="\t", index=False)
full_ids = set(pd.read_csv(os.path.join(B, "N11_out", "ring3_full", "published_coverage.tsv"), sep="\t").query("full == True")["sig_id"])


def cp(x, n):
    lo = 0.0 if x == 0 else beta.ppf(0.025, x, n - x + 1)
    hi = 1.0 if x == n else beta.ppf(0.975, x + 1, n - x)
    return lo, hi


def comp(re, nu, lab, filt):
    xr, nr, xn, nn = int(re["M9_pos"].sum()), len(re), int(nu["M9_pos"].sum()), len(nu)
    lo, hi = cp(xn, nn); orr, p = fisher_exact([[xr, nr - xr], [xn, nn - xn]])
    return dict(comparison=lab, filter=filt, n_real=nr, rate_real=xr / nr, n_null=nn, FPR_null=xn / nn, FPR_CP_lo=lo, FPR_CP_hi=hi,
                typeI_inflated=lo > 0.075, fisher_OR=orr, fisher_P=p,
                coverage_median_real=re["coverage"].median(), coverage_median_null=nu["coverage"].median())


out = []
for filt in ["all", "n_in_net>=1"]:
    T = S if filt == "all" else S[S["n_in_net"] >= 1]
    A_r = T[(T.arm == "A") & (T.kind == "real")]
    out.append(comp(A_r, T[(T.arm == "A") & (T.kind == "null")], "A_vs_fitsig_k", filt))
    out.append(comp(A_r, T[T.arm == "Ak"], "A_vs_Ak_sizematched", filt))
    out.append(comp(T[(T.arm == "B") & (T.kind == "real")], T[(T.arm == "B") & (T.kind == "null")], "B_vs_fitsig_k", filt))
    out.append(comp(A_r[A_r.src.isin(full_ids)], T[(T.arm == "Ak") & T.src.isin(full_ids)], "Afull_vs_Ak_sizematched", filt))
D = pd.DataFrame(out); D.to_csv(os.path.join(O, "m9_decision.tsv"), sep="\t", index=False)
gsum = pd.DataFrame([dict(genes_knocked=len(GL), genes_pos=len(pos), frac_pos=len(pos) / len(GL),
                          n_sig_median=GL["n_sig"].median(), n_sig_max=GL["n_sig"].max())])
gsum.to_csv(os.path.join(O, "m9_gene_summary.tsv"), sep="\t", index=False)
pd.set_option("display.width", 250); print(gsum.to_string(index=False)); print(D.to_string(index=False)); print("DONE")
