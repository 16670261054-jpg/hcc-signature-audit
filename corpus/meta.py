import json, subprocess, time, re, sys
def curl(url,data=None):
    cmd=["curl","-s","--ssl-no-revoke","-m","90",url]
    if data: cmd=["curl","-s","--ssl-no-revoke","-m","90","-d",data,url]
    for i in range(4):
        r=subprocess.run(cmd,capture_output=True)
        if r.stdout.strip(): return r.stdout.decode("utf8","replace")
        time.sleep(2)
    return ""
P=json.load(open("raw/pmids.json"))
ids=list(P)
meta={}
for k in range(0,len(ids),150):
    b=ids[k:k+150]
    s=curl("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi","db=pubmed&retmode=json&id="+",".join(b))
    j=json.loads(s)["result"]
    for i in b:
        r=j.get(i,{})
        meta[i]={"title":r.get("title",""),"journal":r.get("source",""),"year":r.get("pubdate","")[:4],
          "first_author":(r.get("authors") or [{}])[0].get("name",""),"pmcid":"",
          "pubtype":";".join(r.get("pubtype",[])),"mechs":P[i]}
        for a in r.get("articleids",[]):
            if a["idtype"]=="pmc": meta[i]["pmcid"]=a["value"]
    time.sleep(0.4)
json.dump(meta,open("raw/meta.json","w"),ensure_ascii=False)
bad=re.compile(r"lncRNA|long non-?coding|miRNA|microRNA|circRNA|non-coding|lncRNAs",re.I)
hcc=re.compile(r"hepatocellular|HCC|liver cancer|hepatic carcinoma|liver hepatocellular",re.I)
keep={i:m for i,m in meta.items() if not bad.search(m["title"]) and hcc.search(m["title"])}
print(len(meta),len(keep),sum(1 for m in keep.values() if m["pmcid"]))
json.dump(keep,open("raw/keep.json","w"),ensure_ascii=False)
