import os,re,json,csv
C=json.load(open("raw/cand.json"))
# HGNC
appr={};prev={};alias={}
with open("ref/hgnc.txt",encoding="utf8") as f:
    r=csv.DictReader(f,delimiter="\t")
    for d in r:
        if d["status"]!="Approved": continue
        s=d["symbol"]; appr[s.upper()]=s
        for x in d["prev_symbol"].strip('"').split("|"):
            if x: prev.setdefault(x.upper(),set()).add(s)
        for x in d["alias_symbol"].strip('"').split("|"):
            if x: alias.setdefault(x.upper(),set()).add(s)
def hgnc(g):
    u=g.upper()
    if u in appr: return appr[u],"approved"
    if u in prev and len(prev[u])==1: return list(prev[u])[0],"prev_symbol"
    if u in alias and len(alias[u])==1: return list(alias[u])[0],"alias"
    return "","not_found"
def clean(g):
    g=g.strip("\"“”'.(")
    for suf in ("EXP","Exp","exp"):
        if g.endswith(suf) and g[:-len(suf)].upper() in appr: return g[:-len(suf)]
    for pre in ("EX","Exp","exp"):
        if g.startswith(pre) and g[len(pre):].upper() in appr and g.upper() not in appr: return g[len(pre):]
    return g
num=r"[-−–]?\s*\d*\.\d+|[-−–]?\s*\d+"
SG=r"([+\-−–—]?)"
MUL=r"\s*[×\*✕·∗]\s*"
G=r"[\"“]?([A-Z][A-Za-z0-9\-\.]{1,14})[\"”]?"
pre=r"(?:the\s+)?(?:gene\s+)?(?:mRNA\s+)?(?:expression\s+(?:level\s+|value\s+)?of\s+|Exp(?:ression)?\s*\(?\s*|expr\s*\(?|exp\s*\(?\s*)?"
post=r"\s*\)?\s*(?:exp\.?|Exp\.?|EXP|expression(?:\s+level)?(?:\s*\(exp\.\))?|relative expression)?\s*\)?"
pA=re.compile(SG+r"\s*[\(\[]?\s*\(?\s*[\"“]?("+num+r")[\"”]?\s*\)?"+MUL+r"[\(\[]?\s*"+pre+G+post)
pB=re.compile(SG+r"\s*[\(\[]?\s*"+pre+G+post+MUL+r"[\(\[]?\s*[\"“]?\(?\s*("+num+r")\s*\)?")
stop={"COEF","COEFFICIENT","RISK","RISKSCORE","EXPRESSION","EXP","THE","SCORE","HR","CI","OS","TCGA","LIHC","ICGC","GEO","SUM","EXPR","COEFI","EXPI","PCD","CRG","FRG","AUC","RS","PSRS","CRRS","FRRS"}
def tonum(sign,c):
    c=c.replace("−","-").replace("–","-").replace(" ","")
    v=float(c)
    if sign in ("-","−","–","—") : v=-abs(v) if v>0 else v  # operator minus
    return v
f1=re.compile(r"\d\.\d+[\"”]?\s*\)?\s*[×\*✕·∗]")
f2=re.compile(r"[×\*✕∗]\s*[\"“]?\(?\s*[-−–]?\d*\.\d+")
res={}
for fn in sorted(os.listdir("raw/txt")):
    pmid=fn[:-4]
    if pmid not in C: continue
    L=open("raw/txt/"+fn,encoding="utf8").read().split("\n")
    for l in L:
        if not l.startswith(("P","F")): continue
        if len(f1.findall(l))<2 and len(f2.findall(l))<2: continue
        a=[(s,g,c) for s,c,g in pA.findall(l)]
        b=[(s,g,c) for s,g,c in pB.findall(l)]
        def good(ps): return [x for x in ps if clean(x[1]).upper() not in stop and not re.match(r"^\d",x[1])]
        a=good(a); b=good(b)
        # prefer the set with more HGNC-recognised genes
        sa=sum(1 for x in a if hgnc(clean(x[1]))[0]); sb=sum(1 for x in b if hgnc(clean(x[1]))[0])
        pairs=a if sa>=sb else b
        if len(pairs)>=2 and pmid not in res:
            res[pmid]={"line":l[2:],"pairs":[(clean(g),tonum(s,c),s,c) for s,g,c in pairs]}
json.dump(res,open("raw/parsed2.json","w"),ensure_ascii=False,indent=0)
print(len(res))
