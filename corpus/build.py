import json, re, os, sys, collections
sys.path.insert(0, "scripts")
from hgnc_util import hgnc
C = json.load(open("raw/cand.json")); R = json.load(open("raw/parsed2.json")); A = json.load(open("raw/authors.json"))


def txt(p):
    fn = f"raw/txt/{p}.txt"
    return open(fn, encoding="utf8").read() if os.path.exists(fn) else ""


def tnum(c):
    c = c.strip().replace("−", "-").replace("–", "-").replace(" ", "")
    if c.startswith("."): c = "0" + c
    if c.startswith("-."): c = "-0" + c[1:]
    return float(c)


def table_sig(p, hdr):
    L = txt(p).split("\n"); out = []
    for i, l in enumerate(L):
        if l.startswith("TR") and hdr in l:
            cols = [c.strip() for c in l[3:].split("|")]
            ci = [k for k, c in enumerate(cols) if re.match(r"^(coef|coefficient)$", c, re.I)][0]
            gi = [k for k, c in enumerate(cols) if re.search(r"symbol", c, re.I)]
            gi = gi[0] if gi else 0
            for l2 in L[i + 1:]:
                if not l2.startswith("TR"): break
                cc = [c.strip() for c in l2[3:].split("|")]
                try: out.append((cc[gi], tnum(cc[ci]), l2[3:][:200]))
                except Exception: pass
            return out
    return out


TABLES = {"34746309": ("Id | Coef", "Table(Id/Coef/HR)"), "35293027": ("Gene | Coef", "Table(Gene/Coef)"),
          "35309514": ("Gene | Coef", "Table(Gene/Coef)"), "35782053": ("Ensembl ID | Symbol ID", "Table(Ensembl/Symbol/Coef)"),
          "35813858": ("Gene | Coefficient", "Table 1"), "36212155": ("ID | Coef", "Table(ID/Coef/HR)"),
          "36703728": ("Gene | Coefficient", "Table 1"), "38259960": ("Gene symbol | Description | Coefficient", "Table 1"),
          "39525011": ("GS | Coef", "Table 1"), "39251909": (" | Coef | Exp(coef)", "Table 2"), "37014321": ("Gene | Coef | aHR", "Table(Gene/Coef/aHR)")}
MANUAL = {
    "33351066": [("VPS35", 0.0326), ("VPS26A", 0.1966), ("PRKCD", 0.0213), ("BIRC5", 0.2453), ("HMOX1", 0.1205), ("VEGFA", 0.1554), ("WAC", 0.4405), ("FEZ2", 0.3583)],
    "36611208": [("CASP3", 0.0097), ("IRAK1", 0.0065), ("MAPK1", 0.0053), ("MAPK3", 0.0083), ("YWHAB", 0.0018)],
    "36845733": [("GLS", 0.53), ("DLAT", 0.42), ("CDKN2A", 0.505), ("FDX1", -0.135), ("LIAS", -0.491)],
}
MANUAL_NOTE = {"36611208": "原文写作'0.97 × 10–2'等，系数为按原文科学计数直接换算(0.97e-2=0.0097)，非估算",
               "33351066": "公式为 expression(GENE) 写法，人工抄录",
               "36845733": "原文'0.53*GLS+0.42*DLAT+ 0.505*CDKN2A-0.135*FDX1-0.491*LIAS'，人工抄录"}
EXTRA = {"37502362": [("ANP32B", 0.561)]}
INCOMPLETE = {"35123387": "IL18 系数原文排版缺损为'+ 179 * IL18exp.)'，无法确定(0.179或-0.179)",
              "35651286": "CASP8 系数原文为'×10.0134952699440289'，疑排版错误，无法确定"}
EXCLUDE = {"34497633": "基因对(gene-pair)签名，非单基因 Σcoef×expr",
           "36516498": "ICD 基因对(gene-pair)签名",
           "38495502": "签名含 lncRNA SNHG4，芯片队列无法完整复现",
           "42524197": "有12基因RiskScore但自动抽取未找到系数表，未人工核", "41611937": "NETosis签名含snoRNA SNORA31(非mRNA)，芯片队列难以复现"}
SIGNFLAG = {"34869350": "原文用破折号连接各项，CBS/CD2/PON1/PPARGC1A 的负号由'—'/'−'推断",
            "36276933": "CD40LG 前为'—'(破折号)，按负号记录",
            "36212155": "表中Coef与HR不对应(Coef 0.075对HR 1.62)，HR疑为单因素结果",
            "37370041": "CLN3 系数10.05(HR 23130)，模型极不稳定",
            "37502362": "公式写成 riskScore=e^(...)；RAB328 非 HGNC 符号(疑笔误)；ANP32B 项写作'ANP32B*0.561'，人工补录",
            "38649860": "CYP26B 非完整符号(疑CYP26B1)；签名以m1A+自噬命名",
            "35659304": "原文'CSAP8'，按CASP8拼写错误处理",
            "35966333": "多组学签名，此处只取转录层分数(Scoregene expression)",
            "38259960": "正文公式为3位小数，表1为全精度，采用表1；正文'CPS11'为CPS1之误",
            "34803509": "文中有LASSO(12基因)与多因素Cox(10基因)两个公式，取LASSO公式(文中称其区分更明显)",
            "35847985": "另有ICGC重拟合的4基因公式，取TCGA训练公式",
            "36920176": "另有ICGC重拟合公式，取TCGA训练公式",
            "37653738": "文中第二公式为lncRNA模型已舍弃，取CRG mRNA签名",
            "35422802": "同文另有4基因RFS模型(未收)；此处为OS模型",
            "37462769": "系数极小(1e-4量级)，可能基于未log的表达值",
            "38157255": "签名以铜代谢(copper metabolism)命名而非cuproptosis，标题称SLC27A5调控cuproptosis；归入cuproptosis需谨慎"}


def mech(t):
    t = t.lower(); m = []
    for k, v in [("cuproptos", "cuproptosis"), ("cuprotos", "cuproptosis"), ("entosis", "entosis"), ("cell-in-cell", "entosis"), ("copper metabolism", "cuproptosis"), ("disulfidptos", "disulfidptosis"), ("ferroptos", "ferroptosis"),
                 ("pyroptos", "pyroptosis"), ("necroptos", "necroptosis"), ("panoptos", "PANoptosis"), ("anoikis", "anoikis"),
                 ("autophag", "autophagy"), ("immunogenic cell death", "ICD"), ("parthanatos", "parthanatos"), ("netosis", "NETosis"),
                 ("neutrophil extracellular", "NETosis"), ("necso", "NECSO")]:
        if k in t and v not in m: m.append(v)
    if re.search(r"programmed cell death|program cell death|cell death patterns|cell-death patterns|cell death-related|non-apoptotic|rcd panel|cell death genes", t):
        if not m: m = ["PCD_multi"]
        elif m != ["ICD"]: m = ["PCD_multi"] + m
    return "+".join(m) if m else "other"


WET = re.compile(r"qRT-PCR|RT-qPCR|quantitative real-time PCR|real-time quantitative PCR|western blot|transfect|siRNA|shRNA|CCK-8|Transwell|xenograft", re.I)


def wet(t):
    m = WET.findall(t)
    return ("有(自动:" + ",".join(sorted(set(x.lower() for x in m)))[:80] + ")") if m else "未检出(自动)"


def outcome(p):
    ti = C[p]["title"].lower()
    return "RFS" if re.search(r"recurrence|relapse|disease-free|\brfs\b", ti) else "OS"


def training(t):
    for m in re.finditer(r"train(?:ing)?\s+(?:set|cohort|dataset|group)", t, re.I):
        w = t[max(0, m.start() - 250):m.end() + 250]
        h = re.findall(r"TCGA|ICGC|LIRI-JP|GSE\d{4,6}|E-TABM-36", w)
        if h: return h[0] + "(自动)"
    h = re.findall(r"TCGA|ICGC|GSE\d{4,6}", t)
    return (h[0] + "(自动,无training表述)") if h else "未判读"


LNC = re.compile(r"^(A[CLP]\d{6}\.\d|LINC\d|.*-AS\d|SNHG\d|MIR\d|NRAV$|LUCAT)")
sigs = []; long = []


def snippet(line, g, c=None):
    best = None
    for m in re.finditer(re.escape(g), line):
        w = line[max(0, m.start() - 45):m.end() + 45]
        if c and c.strip().lstrip("-−–").strip() in w: best = m.start(); break
        if best is None and re.search(r"[×*∗]", w): best = m.start()
    i = line.find(g) if best is None else best
    i = max(0, i - 60); return line[i:i + 200]


def add(p, pairs, source):
    t = txt(p); genes = []
    for g, v, q in pairs:
        h, st = hgnc(g)
        if g == "CSAP8": h, st = "CASP8", "manual(typo)"
        genes.append((g, h, st, v, q))
    if any(LNC.match(g) for g, *_ in genes):
        EXCLUDE.setdefault(p, "签名含 lncRNA"); return
    sid = f"S{len(sigs) + 1:03d}"; au = A.get(p, [])
    s = dict(sig_id=sid, PMID=p, PMCID=C[p]["pmcid"], first_author=C[p]["first_author"], last_author=au[-1] if au else "",
             year=C[p]["year"], journal=C[p]["journal"], mechanism=mech(C[p]["title"]), outcome=outcome(p), training=training(t),
             n=len(genes), complete="否" if p in INCOMPLETE else "是", wet=wet(t),
             note="; ".join(x for x in [INCOMPLETE.get(p, ""), SIGNFLAG.get(p, ""), MANUAL_NOTE.get(p, "")] if x),
             title=C[p]["title"], unmapped=",".join(g for g, h, *_ in genes if not h), source=source)
    sigs.append(s)
    for g, h, st, v, q in genes:
        bad = (p == "35123387" and g.startswith("IL18")) or (p == "35651286" and g == "CASP8")
        long.append(dict(sig_id=sid, PMID=p, PMCID=s["PMCID"], first_author=s["first_author"], year=s["year"], journal=s["journal"],
                         mechanism=s["mechanism"], outcome=s["outcome"], training_cohort=s["training"], gene_original=g, gene_hgnc=h,
                         hgnc_status=st, coef="NA" if bad else repr(v), coef_source=source, evidence_quote=q[:200].replace("\t", " ")))


for p in C:
    if p in EXCLUDE: continue
    if not ("2021" <= C[p]["year"] <= "2026"): EXCLUDE[p] = "年份不在2021-2026"; continue
    if p in MANUAL:
        t = txt(p); line = ""
        for l in t.split("\n"):
            if MANUAL[p][0][0] in l and MANUAL[p][1][0] in l and re.search(r"[×*∗]", l): line = l; break
        add(p, [(g, v, snippet(line, g)) for g, v in MANUAL[p]], "正文公式(人工抄录)"); continue
    if p in TABLES:
        rows = table_sig(p, TABLES[p][0])
        if rows: add(p, rows, TABLES[p][1]); continue
    if p in R:
        line = R[p]["line"]
        prs = [(g, v, snippet(line, g, c)) for g, v, s, c in R[p]["pairs"] if g != "HCC"]
        prs += [(g, v, snippet(line, g)) for g, v in EXTRA.get(p, [])]
        add(p, prs, "正文公式(Results/Methods段)")

json.dump({"sigs": sigs, "long": long, "exclude": EXCLUDE}, open("raw/build.json", "w"), ensure_ascii=False, indent=0)
print(len(sigs), len(long))
print(collections.Counter(s["mechanism"] for s in sigs))
for s in sigs:
    if s["unmapped"]: print("UNMAPPED", s["PMID"], s["unmapped"])
print({k: v for k, v in EXCLUDE.items()})
