#!/usr/bin/env python3
"""
PEP permutation module - spatial non-randomness tests.

Implements the two null models used to test proximity and chromosome/ZNF
enrichment claims (Reviewer R1-2, R2-6):

  Null A (chromosome-conditional): shuffle each ENSG locus to a random position
          on its OWN chromosome -> tests general proximity to protein-coding genes.
  Null B (gene-density-matched):   redistribute ENSG loci across chromosomes with
          probability proportional to each chromosome's protein-coding gene count
          -> tests whether per-chromosome enrichment (e.g. chr19) exceeds the
          expectation set by gene density.

All p-values are empirical two-sided / one-sided as documented per test.
Default N_PERM = 10000, seed fixed for reproducibility.
"""
import pandas as pd, numpy as np, argparse

MAIN_CHR = ["chr%d" % i for i in range(1, 23)] + ["chrX", "chrY"]


def chrom_sizes(coords):
    """Max coordinate per chromosome as a proxy for chromosome length."""
    return coords.groupby("chr").end.max().to_dict()


def nearest_dist_random(rng, n_by_chr, pc_tss_by_chr, sizes):
    """Null A: random TSS on same chromosome; return nearest-PC distances."""
    dists = []
    for c, n in n_by_chr.items():
        pcs = pc_tss_by_chr.get(c)
        if pcs is None or len(pcs) == 0 or n == 0:
            continue
        pos = rng.integers(0, int(sizes[c]) + 1, size=n)
        idx = np.searchsorted(pcs, pos)
        idx = np.clip(idx, 1, len(pcs) - 1)
        left = np.abs(pos - pcs[idx - 1]); right = np.abs(pcs[idx] - pos)
        dists.append(np.minimum(left, right))
    return np.concatenate(dists) if dists else np.array([])


def run_proximity_null(universe, pc, n_perm=10000, seed=0, within=50000):
    """Null A: is observed fraction-within-`within` and median distance
    beyond chromosome-conditional random expectation?"""
    rng = np.random.default_rng(seed)
    sizes = {c: universe[universe.chr == c].end.max() if (universe.chr == c).any()
             else 0 for c in MAIN_CHR}
    # better: use pc coords for sizes
    for c in MAIN_CHR:
        m = pc[pc.chr == c]
        if len(m): sizes[c] = max(sizes[c], m.tss.max())
    n_by_chr = universe.chr.value_counts().to_dict()
    pc_tss_by_chr = {c: np.sort(g.tss.values.astype(int))
                     for c, g in pc.groupby("chr")}
    obs_within = float((universe.dist_tss <= within).mean())
    obs_median = float(universe.dist_tss.median())
    null_within = np.empty(n_perm); null_median = np.empty(n_perm)
    for i in range(n_perm):
        d = nearest_dist_random(rng, n_by_chr, pc_tss_by_chr, sizes)
        null_within[i] = (d <= within).mean()
        null_median[i] = np.median(d)
    p_within = (1 + np.sum(null_within >= obs_within)) / (n_perm + 1)
    p_median = (1 + np.sum(null_median <= obs_median)) / (n_perm + 1)
    return dict(obs_within=obs_within, null_within_mean=null_within.mean(),
                p_within=p_within, obs_median=obs_median,
                null_median_mean=null_median.mean(), p_median=p_median,
                null_within=null_within, null_median=null_median)


def run_density_null(universe, pc, n_perm=10000, seed=0):
    """Null B: per-chromosome ENSG fraction vs gene-density-matched expectation."""
    rng = np.random.default_rng(seed)
    pc_count = pc.chr.value_counts().reindex(MAIN_CHR).fillna(0)
    weights = (pc_count / pc_count.sum()).values
    n = len(universe)
    obs_frac = universe.chr.value_counts().reindex(MAIN_CHR).fillna(0).values / n
    null = np.zeros((n_perm, len(MAIN_CHR)))
    for i in range(n_perm):
        draw = rng.multinomial(n, weights)
        null[i] = draw / n
    exp_frac = weights
    p = np.array([(1 + np.sum(np.abs(null[:, j] - exp_frac[j]) >=
                              abs(obs_frac[j] - exp_frac[j]))) / (n_perm + 1)
                  for j in range(len(MAIN_CHR))])
    out = pd.DataFrame(dict(chr=MAIN_CHR, observed=obs_frac, expected=exp_frac,
                            log2_OE=np.log2((obs_frac + 1e-9) / (exp_frac + 1e-9)),
                            p_value=p))
    return out.sort_values("log2_OE", ascending=False).reset_index(drop=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--universe", required=True, help="ensg_nearest_gene.csv from pep_core")
    ap.add_argument("--de", required=True); ap.add_argument("--coords", required=True)
    ap.add_argument("--nperm", type=int, default=10000); ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()
    from pep_core import protein_coding_ref
    u = pd.read_csv(a.universe)
    de = pd.read_csv(a.de); coords = pd.read_csv(a.coords)
    pc = protein_coding_ref(de, coords)
    prox = run_proximity_null(u, pc, a.nperm, a.seed)
    print(f"[Null A proximity] within50kb obs={prox['obs_within']:.3f} "
          f"null={prox['null_within_mean']:.3f} p={prox['p_within']:.4f}")
    dens = run_density_null(u, pc, a.nperm, a.seed)
    dens.to_csv("perm_chr_enrichment.csv", index=False)
    print("[Null B density] top/bottom chromosomes written to perm_chr_enrichment.csv")
    print(dens[dens.chr == "chr19"].to_string(index=False))
