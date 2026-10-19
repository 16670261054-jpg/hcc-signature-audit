import json,subprocess,time
C=json.load(open("raw/cand.json")); ids=list(C); A={}
for k in range(0,len(ids),150):
    b=ids[k:k+150]
    for t in range(4):
        r=subprocess.run(["curl","-s","--ssl-no-revoke","-m","90","-d","db=pubmed&retmode=json&id="+",".join(b),"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"],capture_output=True)
        if r.stdout.strip(): break
        time.sleep(2)
    j=json.loads(r.stdout)["result"]
    for i in b: A[i]=[a["name"] for a in j[i].get("authors",[])]
json.dump(A,open("raw/authors.json","w"),ensure_ascii=False); print(len(A))
