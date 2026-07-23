#!/usr/bin/env python3
"""
PEP (Proximity-aware ENSG Pipeline) - core proximity module
Identifies low-annotated ENSG loci and their nearest protein-coding gene.

Inputs:
  - DESeq2 results table (gene_id, gene_name, gene_type, baseMean, log2FoldChange, lfcSE, padj)
  - gene coordinate table (gene_id, chr, start, end, strand, tss) from featureCounts GTF

Universe definition (canonical, this study):
  low-annotated = gene_name missing OR ENSG-prefixed OR clone-accession pattern
  expressed     = baseMean > EXPR_THRESHOLD (default 10)
  main chromosomes only (chr1-22, X, Y)

Distance metric: TSS-to-TSS absolute distance to the nearest expressed
protein-coding gene on the SAME chromosome (strand-agnostic).
"""
import pandas as pd, numpy as np, re, argparse

MAIN_CHR = ["chr%d" % i for i in range(1, 23)] + ["chrX", "chrY"]

def is_low_annotated(name):
    if pd.isna(name) or str(name).strip() == "": return True
    n = str(name)
    if n.startswith("ENSG"): return True
    if re.match(r"^[A-Z]{2}\d{6}\.\d+$", n): return True   # clone accessions
    if re.match(r"^A[CLP]\d", n): return True                 # AC/AL/AP clones
    if re.match(r"^Z\d{5}", n): return True
    return False

def build_universe(de, coords, expr=10):
    de = de.copy()
    de["low_name"] = de.gene_name.map(is_low_annotated)
    u = de[(de.low_name) & (de.baseMean > expr)].merge(
        coords[["gene_id","chr","start","end","strand","tss"]], on="gene_id", how="left")
    return u[u.chr.isin(MAIN_CHR)].dropna(subset=["tss"]).reset_index(drop=True)

def protein_coding_ref(de, coords, expr=10):
    pc = de[(de.gene_type=="protein_coding") & (de.baseMean > expr)].merge(
        coords[["gene_id","chr","tss","strand"]], on="gene_id", how="left")
    pc = pc[pc.chr.isin(MAIN_CHR)].dropna(subset=["tss"]).copy()
    pc["tss"] = pc.tss.astype(int)
    return pc

def nearest_pc(universe, pc):
    pc_by_chr = {c: g.reset_index(drop=True) for c, g in pc.groupby("chr")}
    out = []
    for _, row in universe.iterrows():
        sub = pc_by_chr.get(row.chr)
        if sub is None or len(sub)==0:
            out.append((np.nan, np.nan, np.nan)); continue
        i = np.abs(sub.tss.values - int(row.tss)).argmin()
        r = sub.iloc[i]
        out.append((r.gene_name, int(abs(r.tss-int(row.tss))), r.gene_id))
    res = pd.DataFrame(out, columns=["nearest_pc_gene","dist_tss","nearest_pc_id"])
    u = pd.concat([universe.reset_index(drop=True), res], axis=1)
    u["dist_cat"] = pd.cut(u.dist_tss, [-1,10000,50000,100000,np.inf],
                           labels=["<10kb","10-50kb","50-100kb",">100kb"])
    return u

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--de", required=True); ap.add_argument("--coords", required=True)
    ap.add_argument("--expr", type=float, default=10); ap.add_argument("--out", default="ensg_nearest_gene.csv")
    a = ap.parse_args()
    de = pd.read_csv(a.de); coords = pd.read_csv(a.coords)
    u = build_universe(de, coords, a.expr)
    pc = protein_coding_ref(de, coords, a.expr)
    u = nearest_pc(u, pc)
    u.to_csv(a.out, index=False)
    print(f"universe={len(u)}  within50kb={100*(u.dist_tss<=50000).mean():.1f}%")
