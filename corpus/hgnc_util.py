import csv
appr = {}; prev = {}; alias = {}
with open("ref/hgnc.txt", encoding="utf8") as f:
    for d in csv.DictReader(f, delimiter="\t"):
        if d["status"] != "Approved":
            continue
        s = d["symbol"]; appr[s.upper()] = s
        for x in d["prev_symbol"].strip('"').split("|"):
            if x: prev.setdefault(x.upper(), set()).add(s)
        for x in d["alias_symbol"].strip('"').split("|"):
            if x: alias.setdefault(x.upper(), set()).add(s)


def hgnc(g):
    u = g.upper()
    if u in appr:
        return appr[u], ("approved" if appr[u] == g else "approved(case_fixed)")
    if u in prev and len(prev[u]) == 1:
        return list(prev[u])[0], "prev_symbol"
    if u in alias and len(alias[u]) == 1:
        return list(alias[u])[0], "alias"
    return "", "not_found"
