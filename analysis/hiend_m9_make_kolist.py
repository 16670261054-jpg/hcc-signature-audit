# -*- coding: utf-8 -*-
"""M9 KO 基因清单（预注册 N20 + 变更 2026-10-17）
用法：python hiend_m9_make_kolist.py <网络目录名> pilot|formal
  pilot ：已发表签名基因中随机 10 个 + U 池随机 10 个（都须在网络内），种子 20261017
  formal：所有签名中落在网络内的全部基因，去重
"""
import os
import sys
import random
import pandas as pd
B = os.path.join(os.environ.get("HCC_AUDIT_ROOT", "."), "N11_out") + os.sep
net, mode = sys.argv[1], sys.argv[2]
G = pd.read_csv(B + "hiend_m8/signature_genes_all_arms.tsv", sep="\t", keep_default_na=False, na_values=[""])
assert G["kind"].isin(["real", "null"]).all() and G["kind"].eq("null").any()
W = set(open(B + net + "/genes_after_qc.txt").read().split())
if mode == "formal":
    ko = pd.DataFrame({"gene": sorted(set(G.loc[G.gene.isin(W), "gene"]))})
    ko.to_csv(B + net + "/ko_list_formal.tsv", sep="\t", index=False)
else:
    rng = random.Random(20261017)
    pub = sorted(set(G.loc[(G.arm == "A") & (G.kind == "real") & G.gene.isin(W), "gene"]))
    U = [g for g in open(B + "hiend_m8/U_pool.txt").read().split() if g in W]
    p = rng.sample(pub, 10); r = rng.sample(sorted(set(U) - set(p)), 10)
    ko = pd.DataFrame({"gene": p + r, "kind": ["published"] * 10 + ["random"] * 10})
    ko.to_csv(B + net + "/ko_list_pilot.tsv", sep="\t", index=False)
print(net, mode, "KO genes:", len(ko), "| network genes:", len(W))
