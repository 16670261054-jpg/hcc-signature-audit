import os, re
from lxml import etree
def norm(s): return re.sub(r"\s+"," ",s or "").strip()
for fn in sorted(os.listdir("raw/bioc")):
    pmid=fn.split(".")[0]
    try: root=etree.parse("raw/bioc/"+fn, etree.XMLParser(recover=True,huge_tree=True)).getroot()
    except Exception as e: print("ERR",fn,e); continue
    lines=[]
    if fn.endswith(".jats.xml"):
        body=root.find(".//body")
        if body is None: lines.append("##NOBODY")
        for el in root.iter():
            tag=el.tag if isinstance(el.tag,str) else ""
            if tag=="p" and not any(a.tag in ("table-wrap","ref-list","fn") for a in el.iterancestors()):
                lines.append("P\t"+norm("".join(el.itertext())))
            elif tag in ("disp-formula",):
                lines.append("F\t"+norm("".join(el.itertext())))
            elif tag=="table-wrap":
                lab=norm("".join(el.findtext("label") or ""))+" "+norm("".join(el.find("caption").itertext()) if el.find("caption") is not None else "")
                lines.append("T\t"+lab[:300])
                for tr in el.iter("tr"):
                    lines.append("TR\t"+" | ".join(norm("".join(c.itertext())) for c in tr if isinstance(c.tag,str)))
            elif tag=="supplementary-material":
                lines.append("S\t"+norm("".join(el.itertext()))[:300]+" "+str([v for k,v in el.attrib.items() if 'href' in k]))
    else:
        for p in root.iter("passage"):
            typ=p.findtext("infon[@key='type']") or ""
            if typ=="ref": continue
            if typ=="table":
                x=p.findtext("infon[@key='xml']")
                lines.append("T\t"+norm(p.findtext("text"))[:200])
                if x:
                    try:
                        t=etree.fromstring(x.encode(),etree.XMLParser(recover=True))
                        for tr in t.iter("tr"):
                            lines.append("TR\t"+" | ".join(norm("".join(c.itertext())) for c in tr if isinstance(c.tag,str)))
                    except Exception: lines.append("TR\t"+norm(p.findtext("text")))
                else: lines.append("TR\t"+norm(p.findtext("text")))
            else:
                lines.append(("P" if typ!="table_caption" else "T")+"\t"+norm(p.findtext("text")))
    open(f"raw/txt/{pmid}.txt","w",encoding="utf8").write("\n".join(lines))
print(len(os.listdir("raw/txt")))
