import os,re,json
C=json.load(open("raw/cand.json"))
num=r"[-−–]?\s*\d+(?:\.\d+)?(?:\s*[eE×x]\s*10\^?\s*[-−–]?\d+|[eE][-−–+]?\d+)?"
gene=r"[A-Z][A-Za-z0-9\-\.]{1,14}"
MUL=r"\s*[×\*✕·]\s*"
pA=re.compile(r"([+\-−–]?)\s*\(?\s*\(?("+num+r")\s*\)?"+MUL+r"\(?\s*(?:the\s+)?(?:expression\s+(?:level\s+)?of\s+|Exp(?:ression)?\s*\(?|expr\s*\(?)?("+gene+r")")
pB=re.compile(r"([+\-−–]?)\s*\(?\s*(?:Exp(?:ression)?\s*\(?|expr\s*\(?)?("+gene+r")\)?\s*(?:expression|exp|Exp)?\s*(?:level)?\s*\)?"+MUL+r"\(?\s*("+num+r")\s*\)?")
def fnum(s):
    s=s.replace("−","-").replace("–","-").replace(" ","")
    return s
f1=re.compile(r"\d\.\d+\s*\)?\s*[×\*✕·]\s*\(?\s*(?:expression\s*(?:level\s*)?of\s*|Exp\s*)?\(?[A-Z]")
f2=re.compile(r"[A-Z][A-Z0-9\-]{1,12}\)?\s*[×\*✕]\s*\(?\s*[-−–]?\d+\.\d+")
stop={"Coef","Coefficient","Risk","RiskScore","Expression","Exp","The","Score","HR","CI","OS","TCGA","LIHC","ICGC","GEO","Sum","Expr","Coefi","Coefficienti","Expi","PCD","CRG","FRG"}
out=[]
for fn in sorted(os.listdir("raw/txt")):
    pmid=fn[:-4]; L=open("raw/txt/"+fn,encoding="utf8").read().split("\n")
    for l in L:
        if not l.startswith(("P","F")): continue
        a=f1.findall(l); b=f2.findall(l)
        if len(a)<2 and len(b)<2: continue
        if len(a)>=len(b):
            pairs=[(s,g,c) for s,c,g in pA.findall(l)]
        else:
            pairs=[(s,g,c) for s,g,c in pB.findall(l)]
        pairs=[(s,g.rstrip(".-"),c) for s,g,c in pairs if g.rstrip(".-") not in stop]
        if len(pairs)>=2:
            out.append({"pmid":pmid,"line":l[2:],"pairs":pairs})
json.dump(out,open("raw/parsed.json","w"),ensure_ascii=False,indent=0)
print(len(out),len(set(o["pmid"] for o in out)))
with open("raw/review.txt","w",encoding="utf8") as f:
    for o in out:
        t=C[o["pmid"]]["title"][:110]
        f.write(f"=={o['pmid']} | {t}\n")
        # snippet around first pair
        l=o["line"]; i=max(0,l.find(o["pairs"][0][1])-120)
        f.write("  RAW: "+l[i:i+700]+"\n")
        f.write("  PARSED: "+"; ".join(f"{g}:{s}{fnum(c)}" for s,g,c in o["pairs"])+"\n")
