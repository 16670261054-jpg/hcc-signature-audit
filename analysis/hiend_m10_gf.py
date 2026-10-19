#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""M10 Geneformer genome-wide in silico deletion (audited module)

Audited configuration = the template Geneformer perturbation pipeline (parameters and order reproduced):
  Pretrained V1-10M; TranscriptomeTokenizer(model_version="V1", custom_attr_name_dict={"Type": "Type"});
  EmbExtractor(emb_mode="cell", cell_emb_style="mean_pool", emb_layer=-1);
  InSilicoPerturber(perturb_type="delete", genes_to_perturb="all", start=Disease -> goal=Control, max_ncells=1000, emb_layer=-1);
  InSilicoPerturberStats(mode="goal_state_shift") -> Goal_end_FDR; template positive = >= 1 signature gene with Goal_end_FDR < 0.05
Mapping: Disease = tumour hepatocytes, Control = adjacent normal hepatocytes (GSE149614).
Per-gene results do not depend on signature membership -> run once for all genes, then look up signatures.

Steps: build / tokenize / check / embs / perturb / stats
  perturb reads GF_MAXN (50 for a test run, 1000 = template value) and GF_TAG (output subdirectory)
Windows: nproc=1, executable code inside __main__.
"""
import os
import sys
import gzip
import time
import pickle
import warnings
warnings.filterwarnings("ignore")

B = os.environ.get("HCC_AUDIT_ROOT", ".")
GFD = os.path.join(B, "N11_out", "hiend_m10", "gf")
H5AD_DIR = os.path.join(GFD, "h5ad")
TOKD = os.path.join(GFD, "tokenized")
DATASET = os.path.join(TOKD, "hcc_hep.dataset")
MODEL_PATH = os.environ.get("GENEFORMER_MODEL_DIR", "models/Geneformer-V1-10M")
STATES = {"state_key": "Type", "start_state": "Disease", "goal_state": "Control", "alt_states": []}
EMB_LAYER = -1
EMB_DIR = os.path.join(GFD, "state_embs")
STATE_PKL = os.path.join(EMB_DIR, "hcc_hep.pkl")
FBS = int(os.environ.get("GF_FBS", "32"))      # affects memory/speed only, not results (template 100)
MAXN = int(os.environ.get("GF_MAXN", "50"))
PTYPE = os.environ.get("GF_PTYPE", "delete")   # delete (template) / overexpress (extension)
TAG = os.environ.get("GF_TAG", f"{PTYPE}_all_n{MAXN}")
ISP_DIR = os.path.join(GFD, f"isp_{TAG}")
STATS_DIR = os.path.join(GFD, f"stats_{TAG}")


def say(*a):
    print(f"[{time.strftime('%F %T')}]", *a, flush=True)


def build():
    import numpy as np
    import pandas as pd
    import scanpy as sc
    import scipy.sparse as sp
    from scipy.io import mmread
    from geneformer import ENSEMBL_DICTIONARY_FILE_30M, TOKEN_DICTIONARY_FILE_30M
    feats = pd.read_csv(os.path.join(GFD, "features.tsv"), sep="\t")
    meta = pd.read_csv(os.path.join(GFD, "cell_meta.tsv"), sep="\t", index_col=0)
    with gzip.open(os.path.join(GFD, "counts.mtx.gz"), "rb") as fh:
        m = mmread(fh).T.tocsr()
    assert m.shape == (len(meta), len(feats)), m.shape
    name2id = pickle.load(open(ENSEMBL_DICTIONARY_FILE_30M, "rb"))
    token_dict = pickle.load(open(TOKEN_DICTIONARY_FILE_30M, "rb"))
    probe = list(name2id.items())[:3]
    say("name2id probe:", probe)
    if all(str(k).startswith("ENSG") for k, _ in probe):
        name2id = {v: k for k, v in name2id.items()}
    sym = feats["gene_symbol"].astype(str).values
    ensg = np.array([name2id.get(s, "") for s in sym], dtype=object)
    keep = np.array([(e != "") and (e in token_dict) for e in ensg])
    say(f"genes {len(sym)} -> mapped & in V1 vocab {keep.sum()}")
    m, ensg, sym = m[:, keep], ensg[keep], sym[keep]
    codes, uniq = pd.factorize(pd.Series(ensg))
    if len(uniq) < len(ensg):
        agg = sp.csr_matrix((np.ones(len(codes)), (np.arange(len(codes)), codes)), shape=(len(codes), len(uniq)))
        m = (m @ agg).tocsr()
        first = pd.Series(range(len(sym))).groupby(codes).first()
        sym, ensg = sym[first.values], np.asarray(uniq, dtype=object)
        say(f"collapsed duplicate ENSG -> {len(uniq)}")
    var = pd.DataFrame({"gene_symbol": sym, "ensembl_id": ensg}, index=pd.Index(ensg))
    obs = meta.copy()
    obs["n_counts"] = obs["nCount"].astype(float).values   # whole-transcriptome library size
    obs["Type"] = obs["Type"].astype(str)
    os.makedirs(H5AD_DIR, exist_ok=True)
    ad = sc.AnnData(X=m, obs=obs, var=var)
    say("Type counts:", dict(ad.obs["Type"].value_counts()))
    ad.write_h5ad(os.path.join(H5AD_DIR, "hcc_hep.h5ad"))
    pd.DataFrame({"gene_symbol": sym, "ensembl_id": ensg}).to_csv(os.path.join(GFD, "gene_anno.tsv"), sep="\t", index=False)
    say("OK build")


def tokenize():
    from geneformer import TranscriptomeTokenizer
    os.makedirs(TOKD, exist_ok=True)
    tk = TranscriptomeTokenizer(custom_attr_name_dict={"Type": "Type", "donor": "donor"},
                                nproc=1, chunk_size=512, model_version="V1", collapse_gene_ids=True)
    tk.tokenize_data(H5AD_DIR, TOKD, "hcc_hep", file_format="h5ad")
    say("OK tokenize ->", DATASET)


def check():
    import collections
    import numpy as np
    from datasets import load_from_disk
    ds = load_from_disk(DATASET)
    say(ds)
    say(collections.Counter(ds["Type"]))
    ln = np.array(ds["length"])
    say(f"length median={np.median(ln):.0f} min={ln.min()} max={ln.max()}")


def embs():
    from geneformer import EmbExtractor
    os.makedirs(EMB_DIR, exist_ok=True)
    if os.path.exists(STATE_PKL):
        say("state_embs exists, skip"); return
    ee = EmbExtractor(model_type="Pretrained", num_classes=0, emb_mode="cell", cell_emb_style="mean_pool",
                      filter_data=None, max_ncells=1000, emb_layer=EMB_LAYER, summary_stat="exact_mean",
                      forward_batch_size=FBS, nproc=1, model_version="V1")
    ee.get_state_embs(STATES, MODEL_PATH, DATASET, EMB_DIR, "hcc_hep")
    say("OK embs ->", STATE_PKL)


def perturb():
    from geneformer import InSilicoPerturber
    state_embs = pickle.load(open(STATE_PKL, "rb"))
    say("state keys:", list(state_embs.keys()), "MAXN", MAXN, "FBS", FBS, "->", ISP_DIR)
    os.makedirs(ISP_DIR, exist_ok=True)
    t0 = time.time()
    isp = InSilicoPerturber(perturb_type=PTYPE, genes_to_perturb="all", model_type="Pretrained", num_classes=0,
                            emb_mode="cell", cell_emb_style="mean_pool", cell_states_to_model=STATES,
                            state_embs_dict=state_embs, max_ncells=MAXN, emb_layer=EMB_LAYER,
                            forward_batch_size=FBS, nproc=1, model_version="V1")
    isp.perturb_data(MODEL_PATH, DATASET, ISP_DIR, f"all_genes_{TAG}")
    say(f"OK perturb {time.time() - t0:.0f} s")


def stats():
    from geneformer import InSilicoPerturberStats
    os.makedirs(STATS_DIR, exist_ok=True)
    st = InSilicoPerturberStats(mode="goal_state_shift", genes_perturbed="all", combos=0, anchor_gene=None,
                                cell_states_to_model=STATES, model_version="V1")
    st.get_stats(ISP_DIR, None, STATS_DIR, f"all_genes_{TAG}_stats")
    say("OK stats ->", STATS_DIR)


if __name__ == "__main__":
    for s in sys.argv[1:]:
        say("== step", s)
        globals()[s]()
