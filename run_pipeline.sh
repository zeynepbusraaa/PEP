#!/usr/bin/env bash
# PEP full pipeline runner. Edit the input paths below, then: bash run_pipeline.sh
set -euo pipefail

DE=${DE:-data_example/deseq2_results.csv}
COORDS=${COORDS:-data_example/gene_coords.csv}
EXPR=${EXPR:-10}
NPERM=${NPERM:-10000}
SEED=${SEED:-0}
mkdir -p results

echo "[1/2] Building ENSG universe + nearest protein-coding gene..."
python scripts/pep_core.py --de "$DE" --coords "$COORDS" --expr "$EXPR" \
    --out results/ensg_nearest_gene.csv

echo "[2/2] Running permutation tests (N=$NPERM, seed=$SEED)..."
cd scripts
PYTHONPATH=. python pep_permutation.py \
    --universe ../results/ensg_nearest_gene.csv \
    --de "../$DE" --coords "../$COORDS" --nperm "$NPERM" --seed "$SEED"
mv -f perm_chr_enrichment.csv ../results/ 2>/dev/null || true
cd ..
echo "Done. Outputs in results/"
