#!/bin/zsh

set -e

BASE="work/codon/soluble_control_batches"
OUTBASE="work/topology/deeptmhmm/soluble_batches"

for i in $(seq -w 2 13)
do
    echo "=== START batch${i} ==="

    rm -rf biolib_results

    biolib run 'DTU/DeepTMHMM:1.0.24' \
        --fasta "${BASE}/soluble_batch_${i}.faa"

    mkdir -p "${OUTBASE}/batch${i}"

    cp biolib_results/deeptmhmm_results.md \
        "${OUTBASE}/batch${i}/"

    cp biolib_results/predicted_topologies.3line \
        "${OUTBASE}/batch${i}/"

    cp biolib_results/TMRs.gff3 \
        "${OUTBASE}/batch${i}/"

    n=$(grep -c '^>' "${OUTBASE}/batch${i}/predicted_topologies.3line")

    echo "batch${i} sequences=${n}"
    echo "=== DONE batch${i} ==="
    echo
done

echo "ALL BATCHES FINISHED"
