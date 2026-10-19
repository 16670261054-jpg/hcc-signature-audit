import json, subprocess, time, os
from concurrent.futures import ThreadPoolExecutor
C=json.load(open("raw/cand.json"))
def get(pmid):
    pmc=C[pmid]["pmcid"]
    if not pmc: return
    fn=f"raw/bioc/{pmid}.xml"
    if os.path.exists(fn) and os.path.getsize(fn)>5000: return
    for i in range(3):
        r=subprocess.run(["curl","-s","--ssl-no-revoke","-m","120",f"https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_xml/{pmc}/unicode"],capture_output=True)
        if len(r.stdout)>5000 and b"<collection" in r.stdout[:2000]:
            open(fn,"wb").write(r.stdout); return
        time.sleep(3)
    # fallback efetch pmc
    r=subprocess.run(["curl","-s","--ssl-no-revoke","-m","120",f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id={pmc}"],capture_output=True)
    if len(r.stdout)>5000: open(f"raw/bioc/{pmid}.jats.xml","wb").write(r.stdout)
with ThreadPoolExecutor(4) as ex: list(ex.map(get,list(C)))
print("done",len(os.listdir("raw/bioc")))
