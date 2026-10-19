import os,re,json
C=json.load(open("raw/cand.json"))
# formula: number followed by multiplication then gene, or gene × number
f1=re.compile(r"[-−–]?\s*\(?\s*[-−–]?\d+\.\d+\s*\)?\s*[×\*✕·x]\s*\(?\s*(?:expression\s*(?:level\s*)?of\s*|Exp\s*)?\(?[A-Z][A-Z0-9\-]{1,12}",)
f2=re.compile(r"[A-Z][A-Z0-9\-]{1,12}\)?\s*[×\*✕]\s*\(?\s*[-−–]?\d+\.\d+")
res={}
for fn in os.listdir("raw/txt"):
    pmid=fn[:-4]; L=open("raw/txt/"+fn,encoding="utf8").read().split("\n")
    form=[l for l in L if l.startswith(("P","F")) and (len(f1.findall(l))>=2 or len(f2.findall(l))>=2)]
    tab=[]
    for i,l in enumerate(L):
        if l.startswith("TR") and re.search(r"coef|Coef|COEF|β|Beta|beta",l):
            tab.append(i)
    res[pmid]={"formula":len(form),"tabcoef":len(tab),"nobody":L[0]=="##NOBODY" if L else True}
nf=sum(1 for v in res.values() if v["formula"]); nt=sum(1 for v in res.values() if v["tabcoef"] and not v["formula"])
print("formula",nf,"table-only",nt,"none",sum(1 for v in res.values() if not v["formula"] and not v["tabcoef"]))
json.dump(res,open("raw/scan.json","w"))
