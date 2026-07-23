# PEP — Proximity-aware ENSG Pipeline

Computational framework accompanying:

> Aksoy ZB. *Proximity-Aware Transcriptomics Reveals Non-Random Genomic
> Clustering of Low-Annotated ncRNAs Near Protein-Coding Genes in SMA iPSCs.*
> Human Molecular Genetics (2026).

PEP interprets **low-annotated ncRNA loci** — transcripts represented only by an
Ensembl gene identifier (ENSG) with no curated symbol — by their genomic
position relative to protein-coding genes, rather than discarding them as noise.
From a DESeq2 result table and a per-gene coordinate table it builds the set of
expressed low-annotated loci, assigns each its nearest protein-coding gene, and
tests whether the observed proximity and per-chromosome distribution exceed
chance under two explicit null models.

## Pipeline

**`scripts/pep_core.py`** — builds the analysis universe and proximity table.
- *Low-annotated* locus: `gene_name` missing, `ENSG`-prefixed, or a clone
  accession (AC/AL/AP/Z or `[A-Z]{2}\d{6}\.\d+`).
- *Expressed*: `baseMean > EXPR_THRESHOLD` (default 10).
- Main chromosomes only (chr1–22, X, Y).
- Distance metric: absolute TSS-to-TSS distance to the nearest expressed
  protein-coding gene on the same chromosome (strand-agnostic).
- Output: `ensg_nearest_gene.csv` (locus, nearest PC gene, distance, distance
  category <10 kb / 10–50 kb / 50–100 kb / >100 kb).

**`scripts/pep_permutation.py`** — two null models (N = 10,000, fixed seed).
- *Null A — chromosome-conditional*: each ENSG locus is shuffled to a random
  position on its own chromosome; tests whether the observed fraction within
  50 kb and the median nearest-gene distance exceed chance.
- *Null B — gene-density-matched*: ENSG loci are redistributed across
  chromosomes with probability proportional to each chromosome's protein-coding
  gene count; tests whether per-chromosome enrichment (e.g. chr19 KRAB-ZNF)
  exceeds the expectation set by gene density. Output: `perm_chr_enrichment.csv`.

**`run_pipeline.sh`** — runs both steps end-to-end into `results/`.

## Installation

```bash
git clone https://github.com/zeynepbusraaa/PEP.git
cd PEP
conda env create -f environment.yml   # or: pip install -r requirements.txt
conda activate pep
```

Only `pandas` and `numpy` are required.

## Usage

Reproduce the paper's results directly from the bundled example inputs:

```bash
bash run_pipeline.sh
```

or run the steps manually / on your own data:

```bash
python scripts/pep_core.py \
    --de data_example/deseq2_results.csv \
    --coords data_example/gene_coords.csv \
    --expr 10 --out results/ensg_nearest_gene.csv

cd scripts
PYTHONPATH=. python pep_permutation.py \
    --universe ../results/ensg_nearest_gene.csv \
    --de ../data_example/deseq2_results.csv \
    --coords ../data_example/gene_coords.csv \
    --nperm 10000 --seed 0
```

Environment overrides for `run_pipeline.sh`: `DE`, `COORDS`, `EXPR`, `NPERM`, `SEED`.

## Expected output (bundled example, seed 0)

```
universe=691  within50kb=76.8%
[Null A proximity] within50kb obs=0.768 null≈0.42 p<0.005
[Null B density]   chr19 log2(O/E)≈-0.05  p≈0.56  (not enriched)
```

The chromosome-conditional test confirms genome-wide proximity to
protein-coding genes; Null B shows the apparent chr19 signal does not exceed
gene-density expectation — the transparency result reported in the response
letter.

## Reproducibility

All randomised procedures take an explicit `--seed`. Gene identity, biotype and
coordinates are read from the supplied tables (derived from the pinned
quantification GTF), so results are reproducible given the same inputs.

## Repository layout

```
PEP/
├── run_pipeline.sh
├── scripts/
│   ├── pep_core.py            # universe + nearest protein-coding gene
│   └── pep_permutation.py     # Null A (proximity) + Null B (density)
├── data_example/
│   ├── deseq2_results.csv     # DESeq2 table used in the paper
│   ├── gene_coords.csv        # per-gene coordinates
│   └── README.md
├── results/                   # created at runtime
├── environment.yml  requirements.txt
├── CITATION.cff  LICENSE  README.md
```

## License

MIT — see [LICENSE](LICENSE).

## Citation

Please cite the article and the archived release (see [CITATION.cff](CITATION.cff));
the Zenodo DOI is minted on release.
