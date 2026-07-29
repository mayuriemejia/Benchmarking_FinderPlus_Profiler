#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# generate_samplesheet.sh
# Genera el samplesheet.tsv (sample_id, r1, r2, organism) para el Paso 4
# a partir de los R1/R2 ya producidos por insilicoseq/main.nf.
#
# Espera la estructura: {base_dir}/{ds1..ds4,control}/reads/sim{1..10}/
#                            {dataset}_sim{N}_R1.fastq.gz
#                            {dataset}_sim{N}_R2.fastq.gz
#
# organism se deja SIEMPRE vacío: son metagenomas espikados con varias
# especies mezcladas, no tiene sentido un único --organism para AMRFinderPlus.
#
# Uso:
#   bash generate_samplesheet.sh \
#       /media/sequentia/isilon/students/visitor6/monkey_d_benchmark/datasets \
#       samplesheet_full.tsv
# ============================================================================

BASE_DIR="${1:?Uso: generate_samplesheet.sh <base_dir> <output.tsv>}"
OUT_TSV="${2:?Uso: generate_samplesheet.sh <base_dir> <output.tsv>}"

echo -e "sample_id\tr1\tr2\torganism" > "$OUT_TSV"

n_found=0
n_missing_pair=0

for r1 in "$BASE_DIR"/{ds1,ds2,ds3,ds4,control}/reads/sim*/*_R1.fastq.gz; do
    [[ -f "$r1" ]] || continue

    r2="${r1/_R1.fastq.gz/_R2.fastq.gz}"
    sample_id=$(basename "$r1" _R1.fastq.gz)

    if [[ ! -f "$r2" ]]; then
        echo "  [WARN] Falta R2 para $sample_id (esperado: $r2), se omite" >&2
        n_missing_pair=$((n_missing_pair + 1))
        continue
    fi

    echo -e "${sample_id}\t${r1}\t${r2}\t" >> "$OUT_TSV"
    n_found=$((n_found + 1))
done

n_expected=50
echo "Samplesheet generado: $OUT_TSV"
echo "  $n_found samples encontrados (esperados: $n_expected)"
[[ $n_missing_pair -gt 0 ]] && echo "  $n_missing_pair samples con R1 sin su R2 correspondiente (omitidos)"

if [[ $n_found -ne $n_expected ]]; then
    echo "  [WARN] El total no coincide con los $n_expected esperados (4 datasets x 10 sims + control x 10)." >&2
    echo "  Revisa si insilicoseq/main.nf terminó todos los jobs antes de continuar." >&2
fi
