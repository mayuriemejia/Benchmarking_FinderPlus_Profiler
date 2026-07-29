#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

# ─── Argumentos opcionales ───────────────────────────────
# Uso: ./run_pipeline_docker.sh [ruta/samplesheet.tsv]
# Si no se pasa argumento, usa samplesheet.tsv por defecto
SAMPLESHEET="${1:-samplesheet.tsv}"

if [ ! -f "${SAMPLESHEET}" ]; then
    echo "ERROR: No se encuentra el samplesheet: ${SAMPLESHEET}" >&2
    exit 1
fi

# ─── Verificar Docker ────────────────────────────────────
ensure_docker() {
    if docker info >/dev/null 2>&1; then
        return 0
    fi
    echo "Docker no está accesible. Intentando iniciarlo..." >&2
    if command -v systemctl >/dev/null 2>&1; then
        sudo systemctl start docker >/dev/null 2>&1 || true
    fi
    if docker info >/dev/null 2>&1; then return 0; fi
    if command -v service >/dev/null 2>&1; then
        sudo service docker start >/dev/null 2>&1 || true
    fi
    if docker info >/dev/null 2>&1; then return 0; fi
    echo "ERROR: No se puede iniciar Docker." >&2
    exit 1
}

ensure_docker

echo "Lanzando AMRProfiler pipeline con samplesheet: ${SAMPLESHEET}"

exec nextflow run main.nf \
    -profile docker \
    --samplesheet "${SAMPLESHEET}" \
    "$@"
