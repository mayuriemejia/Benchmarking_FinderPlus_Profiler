#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# build_master_truthset.sh
# Concatena datasets/ds{1-4}/results/sim{1-10}/ds{N}_sim{K}_truthset.tsv
# en un único truthset.tsv (mismo header, una sola vez).
#
# El dataset 'control' no tiene truthset (no se le insertó nada) — no
# necesita entrada aquí; eval_benchmark.py ya maneja dataset=='control'
# como truth_genes vacío directamente.
#
# Uso:
#   bash build_master_truthset.sh \
#       /media/sequentia/isilon/students/visitor6/monkey_d_benchmark/datasets \
#       truthset_master.tsv
# ============================================================================

BASE_DIR="${1:?Uso: build_master_truthset.sh <base_dir> <output.tsv>}"
OUT_TSV="${2:?Uso: build_master_truthset.sh <base_dir> <output.tsv>}"

n_files=0
header_written=0

> "$OUT_TSV"

for tsv in "$BASE_DIR"/{ds1,ds2,ds3,ds4}/results/sim*/*_truthset.tsv; do
    [[ -f "$tsv" ]] || continue

    if [[ $header_written -eq 0 ]]; then
        head -n1 "$tsv" > "$OUT_TSV"
        header_written=1
    fi
    tail -n +2 "$tsv" >> "$OUT_TSV"
    n_files=$((n_files + 1))
done

n_expected=40
echo "Truthset maestro generado: $OUT_TSV"
echo "  $n_files archivos concatenados (esperados: $n_expected)"
echo "  $(($(wc -l < "$OUT_TSV") - 1)) filas de inserción totales"

if [[ $n_files -ne $n_expected ]]; then
    echo "  [WARN] El total no coincide con los $n_expected esperados (4 datasets x 10 sims)." >&2
    echo "  Revisa si insert_resistance_genes_v2.py terminó todos los datasets." >&2
fi
