import json, re, os, math, csv, collections
B = json.load(open("raw/build.json")); C = json.load(open("raw/cand.json")); A = json.load(open("raw/authors.json"))
S = json.load(open("raw/scan.json")) if os.path.exists("raw/scan.json") else {}
sigs = B["sigs"]; long = B["long"]; EXC = B["exclude"]

# ---- coef vs HR consistency for table sources (rows: gene|coef|HR|...)
for s in sigs:
    if not s["source"].startswith("Table"): continue
    bad = []
    for r in long:
        if r["PMID"] != s["PMID"]: continue
        cells = [c.strip() for c in r["evidence_quote"].split("|")]
        nums = []
        for c in cells[1:]:
            try: nums.append(float(c.replace("−", "-").replace(",", "")))
            except Exception: pass
        if len(nums) >= 2 and nums[1] > 0 and r["coef"] != "NA":
            if abs(math.exp(float(r["coef"])) - nums[1]) / nums[1] > 0.05: bad.append(r["gene_original"])
    if bad:
        msg = "表中exp(coef)与HR不符的基因:" + ",".join(bad)
        if msg not in s["note"]: s["note"] = (s["note"] + "; " + msg).strip("; ")

# ---- author-group dedup: >=3 shared authors, or same last author + >=2 shared
def norm(a): return a.lower().replace("-", " ").strip()
grp = {}; gid = 0
order = sorted(sigs, key=lambda s: (s["year"], s["PMID"]))
for i, s in enumerate(order):
    a1 = set(map(norm, A.get(s["PMID"], [])))
    for t in order[:i]:
        a2 = set(map(norm, A.get(t["PMID"], [])))
        sh = len(a1 & a2)
        same_last = A.get(s["PMID"]) and A.get(t["PMID"]) and norm(A[s["PMID"]][-1]) == norm(A[t["PMID"]][-1])
        fam1 = set(s["mechanism"].split("+")); fam2 = set(t["mechanism"].split("+"))
        if (sh >= 3 or (same_last and sh >= 2)) and (fam1 & fam2):
            grp[s["PMID"]] = grp[t["PMID"]]; break
    else:
        gid += 1; grp[s["PMID"]] = f"G{gid:03d}"
gcount = collections.Counter(grp.values())
seen = set()
for s in order:
    g = grp[s["PMID"]]; s["author_group"] = g
    if gcount[g] > 1:
        s["group_primary"] = "是" if g not in seen else "否(同组重复,不计)"
    else:
        s["group_primary"] = "是"
    seen.add(g)

# ---- reliability
for s in sigs:
    flags = []
    if s["complete"] != "是": flags.append("系数不全")
    if s["unmapped"]: flags.append("基因名无法映射HGNC:" + s["unmapped"])
    if re.search(r"推断|不对应|不符|不稳定|笔误|拼写|破折号", s["note"]): flags.append("见备注")
    s["reliability"] = "需复核(" + ";".join(flags) + ")" if flags else "正常"

D = "."
def w(fn, rows, cols):
    with open(os.path.join(D, fn), "w", encoding="utf8", newline="") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows: f.write("\t".join(str(r.get(c, "")).replace("\t", " ").replace("\n", " ") for c in cols) + "\n")

colsL = ["sig_id", "PMID", "PMCID", "first_author", "year", "journal", "mechanism", "outcome", "training_cohort", "gene_original",
         "gene_hgnc", "hgnc_status", "coef", "coef_source", "evidence_quote"]
w("signatures_long.tsv", long, colsL)
for s in sigs:
    s["机制"] = s["mechanism"]; s["基因数"] = s["n"]; s["是否有完整系数"] = s["complete"]; s["有无湿实验验证"] = s["wet"]
    s["备注"] = s["note"]
colsI = ["sig_id", "PMID", "PMCID", "first_author", "last_author", "year", "journal", "机制", "outcome", "training", "基因数",
         "是否有完整系数", "有无湿实验验证", "author_group", "group_primary", "reliability", "source", "备注", "title"]
w("signatures_index.tsv", sigs, colsI)

# ---- gene_only
inc = {s["PMID"] for s in sigs}
go = []
for p, m in C.items():
    if p in inc: continue
    if p in EXC:
        reason = "排除:" + EXC[p]
    elif not m.get("pmcid"):
        reason = "无PMC开放全文，未能读到正文/表格"
    elif S.get(p, {}).get("nobody") or (os.path.exists(f"raw/txt/{p}.txt") and open(f"raw/txt/{p}.txt", encoding="utf8").read().startswith("##NOBODY")):
        reason = "PMC仅有摘要/出版商不允许全文下载"
    else:
        reason = "全文未检出逐基因系数(公式只写Σcoef×expr，系数可能仅在图/补充材料)，或非签名研究"
    go.append(dict(PMID=p, PMCID=m.get("pmcid", ""), first_author=m["first_author"], year=m["year"], journal=m["journal"],
                   mechanism=__import__("build_mech").mech(m["title"]) if False else "", reason=reason, title=m["title"]))
w("gene_only.tsv", go, ["PMID", "PMCID", "first_author", "year", "journal", "reason", "title"])
json.dump(sigs, open("raw/sigs_final.json", "w"), ensure_ascii=False, indent=0)
prim = [s for s in sigs if s["group_primary"] == "是"]
print("sigs", len(sigs), "primary", len(prim), "complete&primary", sum(1 for s in prim if s["complete"] == "是"))
print(collections.Counter(s["mechanism"] for s in prim if s["complete"] == "是"))
print("gene_only", len(go), collections.Counter(r["reason"][:12] for r in go))
print("reliability", collections.Counter(s["reliability"][:3] for s in sigs))
print("wet", collections.Counter(s["wet"][:1] for s in sigs))
print("outcome", collections.Counter(s["outcome"] for s in sigs))
for g, n in gcount.items():
    if n > 1: print(g, [(s["PMID"], s["first_author"], s["last_author"], s["mechanism"]) for s in sigs if s["author_group"] == g])
