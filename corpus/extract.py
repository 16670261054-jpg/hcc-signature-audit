import json, re, os, sys, html
from xml.etree import ElementTree as ET
C=json.load(open("raw/cand.json"))
pat=re.compile(r"(risk\s*score|riskscore|risk_score|risk value|risk index|score\s*=|formula|coefficient|coef)",re.I)
out={}
for fn in sorted(os.listdir("raw/bioc")):
    pmid=fn.split(".")[0]
    try: root=ET.parse("raw/bioc/"+fn).getroot()
    except Exception as e: out[pmid]={"err":str(e)}; continue
    hits=[]; tables=[]
    for p in root.iter("passage"):
        typ=""; 
        for inf in p.findall("infon"):
            if inf.get("key")=="type": typ=inf.text or ""
        t=(p.findtext("text") or "")
        if typ.startswith("table") :
            x=p.find("infon[@key='xml']")
            tx=t
            if re.search(r"coef|Coef|β|beta|Beta",tx+(typ)) or re.search(r"coef",(x.text if x is not None else ""),re.I):
                tables.append((typ,tx[:3000]))
            continue
        if typ in ("ref",) : continue
        for sent in re.split(r"(?<=[.;])\s+(?=[A-Z])",t):
            if pat.search(sent) and re.search(r"[×\*x]\s*[A-Z]|[A-Z0-9]{2,}\s*[×\*x(]|\d\.\d{2,}",sent):
                hits.append(sent[:1500])
    out[pmid]={"hits":hits,"tables":tables}
json.dump(out,open("raw/extract.json","w"),ensure_ascii=False,indent=0)
n=sum(1 for v in out.values() if v.get("hits") or v.get("tables"))
print(len(out),n)
