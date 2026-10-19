import json, subprocess, time, urllib.parse, sys, os
OUT=sys.argv[1]
def curl(url):
    for i in range(4):
        r=subprocess.run(["curl","-s","--ssl-no-revoke","-m","60",url],capture_output=True)
        if r.stdout.strip(): return r.stdout.decode("utf8","replace")
        time.sleep(2)
    return ""
mechs={
 "cuproptosis":"(cuproptosis OR cuproptosis-related OR copper-induced cell death)",
 "disulfidptosis":"(disulfidptosis)",
 "ferroptosis":"(ferroptosis)",
 "pyroptosis":"(pyroptosis)",
 "necroptosis":"(necroptosis)",
 "PANoptosis":"(PANoptosis)",
 "anoikis":"(anoikis)",
 "autophagy":"(autophagy OR autophagy-related)",
 "NETosis":"(NETosis OR \"neutrophil extracellular trap\")",
 "parthanatos":"(parthanatos)",
 "PCD":"(\"programmed cell death\" OR \"cell death\" OR \"PCD-related\")",
 "immunogenic_cell_death":"(\"immunogenic cell death\")",
 "cellular_senescence_skip":None,
}
tail=" AND (hepatocellular carcinoma OR HCC) AND (signature OR \"risk model\" OR \"prognostic model\" OR \"risk score\") AND 2021:2026[dp] NOT lncRNA[ti] NOT review[pt]"
log=[]; allp={}
for m,q in mechs.items():
    if not q: continue
    term=q+tail
    url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmode=json&retmax=300&term="+urllib.parse.quote(term)
    s=curl(url)
    try: j=json.loads(s)["esearchresult"]; ids=j["idlist"]; cnt=j["count"]
    except Exception as e: ids=[]; cnt="FAIL"
    log.append((m,term,cnt,len(ids)))
    for i in ids: allp.setdefault(i,[]).append(m)
    time.sleep(0.4)
with open(os.path.join(OUT,"search_log.tsv"),"w",encoding="utf8") as f:
    f.write("mechanism\tquery\tcount\tretrieved\tdate\n")
    for r in log: f.write("\t".join(map(str,r))+"\t2026-10-15\n")
json.dump(allp,open(os.path.join(OUT,"raw","pmids.json"),"w"))
print(log); print(len(allp))
