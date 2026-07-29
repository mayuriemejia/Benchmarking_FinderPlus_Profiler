#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

# ─── Parseo de argumentos ────────────────────────────────
# Uso:
#   ./run_pipeline_docker.sh samplesheet.tsv
#   ./run_pipeline_docker.sh --samplesheet samplesheet.tsv
SAMPLESHEET=""
EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --samplesheet)
            SAMPLESHEET="$2"
            shift 2
            ;;
        --samplesheet=*)
            SAMPLESHEET="${1#*=}"
            shift
            ;;
        -profile|--profile)
            # Ignoramos -profile aquí; ya está hardcodeado abajo
            shift 2
            ;;
        *)
            # Primer argumento posicional sin flag = samplesheet
            if [[ -z "$SAMPLESHEET" && "$1" != -* ]]; then
                SAMPLESHEET="$1"
            else
                EXTRA_ARGS+=("$1")
            fi
            shift
            ;;
    esac
done

# ─── Valor por defecto ───────────────────────────────────
SAMPLESHEET="${SAMPLESHEET:-samplesheet.tsv}"

# ─── Verificar que el samplesheet existe ─────────────────
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

echo "Lanzando pipeline con samplesheet: ${SAMPLESHEET}"

exec nextflow run main.nf \
    -profile docker \
    --samplesheet "${SAMPLESHEET}" \
    "${EXTRA_ARGS[@]+"${EXTRA_ARGS[@]}"}"
