#!/bin/bash

################################################################################
# clean_genomes_dual.sh  (v2.1 — lock + tmpdir único en AMRFinderPlus)
#
# Limpieza y enmascaramiento de genomas AMR para benchmark monkey_d.
# Compatible con: AMRFinderPlus (modo combinado -n -p -g) + AMRProfiler (dual)
#
# Cambios v2 -> v2.1:
#   - Lock de ejecución (flock) para evitar que dos invocaciones del script
#     corran en paralelo sobre el mismo OUTPUT_DIR (causa confirmada de
#     "Could not read file ...sprot.tmp.1.blast" por escrituras concurrentes
#     sobre las mismas rutas compartidas)
#   - AMRFinderPlus ahora escribe su TSV en un tmpdir único por proceso
#     (mktemp) en vez de directamente en $fasta_dir, que es compartido
#     entre genomas y podía colisionar con otra ejecución simultánea
#
# Mejora v1 -> v2 (sin cambios):
#   - Anota cada genoma con Prokka antes de correr AMRFinderPlus
#   - Pasa proteínas (.faa) + GFF limpio + nucleótido (.fna) a AMRFinderPlus
#   - Activa el pipeline completo BLASTP + HMM, recuperando genes intrínsecos
#     divergentes (ampC, catA, arr, fosB, etc.) que el modo solo-nucleótido
#     (-n) no detectaba por caer bajo el umbral de BLASTX sin modelos HMM
#   - El GFF de Prokka se limpia eliminando la sección ##FASTA embebida antes
#     de pasarlo a AMRFinderPlus (--annotation_format prokka)
#
# Uso:
#   cd genomas/
#   bash scripts_clean/clean_genomes_dual.sh
#
# Requisitos: Docker disponible en PATH
################################################################################

set -euo pipefail

# ============================================================================
# CONFIGURACIÓN
# ============================================================================

GENOMES_DIR="$(pwd)/genomas_originales"
OUTPUT_DIR="$(pwd)/control_no_amr"
AMRP_DB="/media/sequentia/isilon/students/visitor6/monkey_d_benchmark/amrprofiler_nf/amrprofiler_db/amrprofiler-main"
FLANK=50

BEDTOOLS_IMG="quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_2"
SEQKIT_IMG="quay.io/biocontainers/seqkit:2.10.0--h9ee0642_0"
AMRFINDER_IMG="staphb/ncbi-amrfinderplus:latest"
AMRPROFILER_IMG="amrprofiler:1.0.0"
PROKKA_IMG="staphb/prokka:latest"

GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

log_info()  { echo -e "${GREEN}▶${NC} $1" >&2; }
log_warn()  { echo -e "${YELLOW}⚠${NC} $1" >&2; }
log_error() { echo -e "${RED}✗${NC} $1" >&2; }

# ============================================================================
# LOCK DE EJECUCIÓN — evita corridas paralelas accidentales sobre OUTPUT_DIR
# ============================================================================

acquire_lock() {
    local lockfile="$OUTPUT_DIR/.pipeline.lock"
    mkdir -p "$OUTPUT_DIR"
    exec 200>"$lockfile"
    if ! flock -n 200; then
        log_error "Ya hay una ejecución de este pipeline en curso (lock: $lockfile)."
        log_error "Si estás seguro de que no hay ninguna corriendo, bórralo a mano: rm -f $lockfile"
        exit 1
    fi
    log_info "Lock adquirido: $lockfile"
}

# ============================================================================
# MAPEO accesión → especie (desde assembly_data_report.jsonl)
# ============================================================================

build_species_map() {
    local jsonl="$GENOMES_DIR/assembly_data_report.jsonl"
    if [[ ! -f "$jsonl" ]]; then
        log_warn "assembly_data_report.jsonl no encontrado, AMRProfiler usará 'Unknown'"
        return 0
    fi
    python3 -c "
import json, sys
with open('$jsonl') as f:
    for line in f:
        try:
            d = json.loads(line)
            acc = d.get('accession','')
            org = d.get('organism',{}).get('organismName','Unknown')
            sp = ' '.join(org.split()[:2])
            print(acc + '\t' + sp)
        except: pass
" > "$OUTPUT_DIR/stats/species_map.tsv"
    log_info "Mapa de especies generado: $(wc -l < "$OUTPUT_DIR/stats/species_map.tsv") entradas"
}

get_species() {
    local acc="$1"
    local sp
    sp=$(grep "^${acc}" "$OUTPUT_DIR/stats/species_map.tsv" 2>/dev/null | cut -f2 | head -1)
    echo "${sp:-Unknown}"
}

# ============================================================================
# FUNCIONES AUXILIARES
# ============================================================================

setup_directories() {
    mkdir -p "$OUTPUT_DIR"/{bed,masked,logs,tsv,stats,amrprofiler,prokka}
    log_info "Directorios creados"
}

validate_fasta() {
    local fasta="$1"
    [[ -f "$fasta" ]] && grep -q "^>" "$fasta"
}

clean_fasta_headers() {
    sed 's/^>\([^ ]*\).*/>\1/' "$1" > "$2"
}

# ============================================================================
# COMBINAR BED: AMRFinderPlus + AMRProfiler → BED final unificado
# ============================================================================

combine_beds() {
    local acc="$1"
    local amrf_bed="$OUTPUT_DIR/bed/${acc}_raw_amrf.bed"
    local amrp_bed="$OUTPUT_DIR/bed/${acc}_raw_amrp.bed"
    local combined="$OUTPUT_DIR/bed/${acc}_raw.bed"

    > "$combined"
    [[ -s "$amrf_bed" ]] && cat "$amrf_bed" >> "$combined"
    [[ -s "$amrp_bed" ]] && cat "$amrp_bed" >> "$combined"

    if [[ -s "$combined" ]]; then
        sort -k1,1 -k2,2n "$combined" -o "$combined"
    fi
}

# ============================================================================
# PASO 0: Anotación con Prokka (necesaria para AMRFinderPlus modo HMM)
# ============================================================================

run_prokka_annotation() {
    local acc="$1"
    local fasta_dir="$2"
    local fasta_name="$3"
    local prokka_out="$OUTPUT_DIR/prokka/${acc}"
    local faa_out="$prokka_out/${acc}.faa"
    local gff_clean="$prokka_out/${acc}_clean.gff"

    if [[ -s "$faa_out" && -s "$gff_clean" ]]; then
        log_info "  │  Prokka ya disponible para ${acc}, reutilizando"
        echo "${prokka_out}"
        return 0
    fi

    mkdir -p "$prokka_out"

    # Ejecutar Prokka en /tmp (disco local) para evitar cuelgues de tbl2asn sobre NFS
    local tmp_out
    tmp_out=$(mktemp -d "/tmp/prokka_${acc}_XXXXXX")
    chmod 777 "$tmp_out"

    log_info "  ├─ [0/2] Anotando con Prokka (en /tmp local, timeout 600s)..."

    if ! timeout 600 docker run --rm --user "$(id -u):$(id -g)" \
        -v "$fasta_dir:/data" \
        -v "$tmp_out:/out" \
        "$PROKKA_IMG" \
        prokka --outdir /out --prefix "$acc" --force --quiet \
        "/data/$fasta_name" \
        2>>"$OUTPUT_DIR/logs/errors.log"; then
        log_warn "  │  Prokka falló o excedió timeout para $acc — AMRFinderPlus correrá solo en modo -n"
        rm -rf "$tmp_out"
        echo ""
        return 1
    fi

    # Copiar resultados de /tmp de vuelta a NFS (ya son propiedad del usuario actual gracias a --user)
    cp "$tmp_out"/* "$prokka_out"/ 2>/dev/null
    rm -rf "$tmp_out"

    # Limpiar sección ##FASTA embebida del GFF (incompatible con AMRFinderPlus)
    if [[ -f "$prokka_out/${acc}.gff" ]]; then
        sed '/^##FASTA/,$d' "$prokka_out/${acc}.gff" > "$gff_clean"
        local n_cds
        n_cds=$(grep -c $'\tCDS\t' "$gff_clean" 2>/dev/null || echo 0)
        log_info "  │  Prokka OK → ${n_cds} CDS anotados"
    else
        log_warn "  │  Prokka no generó .gff para ${acc}"
        echo ""
        return 1
    fi

    echo "${prokka_out}"
    return 0
}

# ============================================================================
# PASO 1: AMRFinderPlus (modo combinado -n -p -g con Prokka)
#         v2.1: salida en tmpdir único por proceso, evita colisiones si dos
#         ejecuciones del pipeline llegan a correr en paralelo por error
# ============================================================================

run_amrfinderplus() {
    local acc="$1"
    local fasta_dir="$2"
    local fasta_name="$3"
    local out_tsv="$OUTPUT_DIR/tsv/${acc}_amrfinder.tsv"
    local raw_bed="$OUTPUT_DIR/bed/${acc}_raw_amrf.bed"

    log_info "  ├─ [1/2] AMRFinderPlus..."

    # Intentar anotación Prokka para modo combinado (BLASTP + HMM)
    local prokka_dir
    prokka_dir=$(run_prokka_annotation "$acc" "$fasta_dir" "$fasta_name") || prokka_dir=""

    # tmpdir único por proceso para la salida de amrfinder (evita colisión
    # con otra ejecución concurrente escribiendo sobre el mismo fasta_dir)
    local tmp_amrf
    tmp_amrf=$(mktemp -d "/tmp/amrfinder_${acc}_XXXXXX")
    chmod 777 "$tmp_amrf"

    local amrfinder_args
    local extra_mount=()

    if [[ -n "$prokka_dir" && \
          -s "$prokka_dir/${acc}.faa" && \
          -s "$prokka_dir/${acc}_clean.gff" ]]; then
        # Modo combinado: nucleótido + proteína + GFF (máxima sensibilidad)
        extra_mount=(-v "$prokka_dir:/protein")
        amrfinder_args=(
            -n "/data/$fasta_name"
            -p "/protein/${acc}.faa"
            -g "/protein/${acc}_clean.gff"
            --annotation_format prokka
            --plus
            -o "/out/${acc}_amrfinder.tsv"
        )
        log_info "  │  Modo: nucleótido + proteína + HMM (Prokka)"
    else
        # Fallback: solo nucleótido
        amrfinder_args=(
            -n "/data/$fasta_name"
            --plus
            -o "/out/${acc}_amrfinder.tsv"
        )
        log_warn "  │  Modo: solo nucleótido (Prokka no disponible)"
    fi

    if ! docker run --rm --user "$(id -u):$(id -g)" \
        -v "$fasta_dir:/data" \
        -v "$tmp_amrf:/out" \
        "${extra_mount[@]}" \
        "$AMRFINDER_IMG" \
        amrfinder "${amrfinder_args[@]}" \
        2>>"$OUTPUT_DIR/logs/errors.log"; then
        log_error "  │  AMRFinderPlus falló para ${acc}"
        rm -rf "$tmp_amrf"
        return 1
    fi

    cp "$tmp_amrf/${acc}_amrfinder.tsv" "$out_tsv"
    rm -rf "$tmp_amrf"

    # TSV → BED (columnas: $2=contig, $3=start, $4=stop, 1-based → 0-based)
    awk 'BEGIN{OFS="\t"} NR>1 && $2!="" {
        s=$3; e=$4
        if (s>e) { t=s; s=e; e=t }
        print $2, s-1, e
    }' "$out_tsv" | sort -k1,1 -k2,2n > "$raw_bed"

    local n_hits
    n_hits=$(wc -l < "$raw_bed")
    log_info "  │  AMRFinderPlus: ${n_hits} hits"
    return 0
}

# ============================================================================
# PASO 2: AMRProfiler
# ============================================================================

run_amrprofiler() {
    local acc="$1"
    local fasta="$2"
    local species="$3"
    local raw_bed="$OUTPUT_DIR/bed/${acc}_raw_amrp.bed"
    local amrp_out="$OUTPUT_DIR/amrprofiler/${acc}"
    local csv_out="$amrp_out/${acc}_final_results_tool1.csv"

    log_info "  ├─ [2/2] AMRProfiler (${species})..."

    mkdir -p "$amrp_out"
    chmod 777 "$amrp_out"

    # Limpiar headers antes de pasar a AMRProfiler
    local fasta_clean="$amrp_out/${acc}_input.fasta"
    local fasta_name
    fasta_name=$(basename "$fasta_clean")
    clean_fasta_headers "$fasta" "$fasta_clean"

    if ! docker run --rm --user "$(id -u):$(id -g)" \
        -v "$AMRP_DB:/amrprofiler_db" \
        -v "$amrp_out:/workdir" \
        -w "/workdir" \
        "$AMRPROFILER_IMG" \
        python /amrprofiler_db/amrprofiler.py \
            "/workdir/$fasta_name" \
            "$species" \
            "/amrprofiler_db/" \
            --threads 2 \
            --identity_threshold 70 \
            --coverage_threshold 70 \
        2>>"$OUTPUT_DIR/logs/errors.log"; then
        log_warn "  │  AMRProfiler falló o sin resultados para ${acc}"
        touch "$raw_bed"
        return 0
    fi

    [[ -f "$amrp_out/final_results_tool1.csv" ]] && \
        mv "$amrp_out/final_results_tool1.csv" "$csv_out"

    if [[ -f "$csv_out" ]] && [[ $(wc -l < "$csv_out") -gt 1 ]]; then
        awk -F'\t' 'BEGIN{OFS="\t"} NR>1 && $1!="" && $13~/^[0-9]/ {
            s=$13+0; e=$14+0
            if (s>e) { t=s; s=e; e=t }
            print $1, s-1, e
        }' "$csv_out" | sort -k1,1 -k2,2n > "$raw_bed"
        local n_hits
        n_hits=$(wc -l < "$raw_bed")
        log_info "  │  AMRProfiler: ${n_hits} hits"
    else
        log_warn "  │  AMRProfiler: sin hits"
        touch "$raw_bed"
    fi

    return 0
}

# ============================================================================
# PASO 3: Enmascaramiento con BED combinado
# ============================================================================

mask_genome() {
    local acc="$1"
    local fasta_dir="$2"
    local fasta_name="$3"
    local fasta="$fasta_dir/$fasta_name"
    local raw_bed="$OUTPUT_DIR/bed/${acc}_raw.bed"

    if [[ ! -s "$raw_bed" ]]; then
        log_warn "  │  Sin hits en ninguna herramienta — copiando sin enmascarar"
        clean_fasta_headers "$fasta" "$OUTPUT_DIR/masked/${acc}_clean.fasta"
        echo "$acc    0   0   sin_amr" >> "$OUTPUT_DIR/stats/summary.txt"
        return 0
    fi

    log_info "  ├─ Extrayendo longitudes de contigs..."
    local genome_file="$OUTPUT_DIR/bed/${acc}.genome"
    docker run --rm --user "$(id -u):$(id -g)" \
        -v "$fasta_dir:/data" \
        "$SEQKIT_IMG" seqkit fx2tab --name --only-id --length "/data/$fasta_name" \
        > "$genome_file" 2>>"$OUTPUT_DIR/logs/errors.log"

    log_info "  ├─ Slop (${FLANK} bp) + merge..."
    local slop_bed="$OUTPUT_DIR/bed/${acc}_slop.bed"
    local final_bed="$OUTPUT_DIR/bed/${acc}_final.bed"

    docker run --rm --user "$(id -u):$(id -g)" \
        -v "$OUTPUT_DIR/bed:/data" \
        "$BEDTOOLS_IMG" bedtools slop \
        -i "/data/${acc}_raw.bed" \
        -g "/data/${acc}.genome" \
        -b "$FLANK" \
        > "$slop_bed" 2>>"$OUTPUT_DIR/logs/errors.log"

    docker run --rm --user "$(id -u):$(id -g)" \
        -v "$OUTPUT_DIR/bed:/data" \
        "$BEDTOOLS_IMG" bedtools merge \
        -i "/data/${acc}_slop.bed" \
        > "$final_bed" 2>>"$OUTPUT_DIR/logs/errors.log"

    log_info "  ├─ Enmascarando FASTA con N's..."
    local masked_fasta="$OUTPUT_DIR/masked/${acc}_masked.fasta"
    docker run --rm --user "$(id -u):$(id -g)" \
        -v "$fasta_dir:/fasta" \
        -v "$OUTPUT_DIR/bed:/bed" \
        -v "$OUTPUT_DIR/masked:/masked" \
        "$BEDTOOLS_IMG" bedtools maskfasta \
        -fi "/fasta/$fasta_name" \
        -bed "/bed/${acc}_final.bed" \
        -fo "/masked/${acc}_masked.fasta" \
        2>>"$OUTPUT_DIR/logs/errors.log"

    log_info "  ├─ Limpiando headers..."
    clean_fasta_headers "$masked_fasta" "$OUTPUT_DIR/masked/${acc}_clean.fasta"
    rm -f "$masked_fasta"

    local n_hits n_bp
    n_hits=$(wc -l < "$final_bed")
    n_bp=$(awk '{sum+=$3-$2} END {print sum+0}' "$final_bed")
    log_info "  └─ ✓ ${n_hits} región(es) enmascarada(s) | ${n_bp} bp"

    echo "$acc    $n_hits $n_bp   $(date +%s)" >> "$OUTPUT_DIR/stats/summary.txt"
}

# ============================================================================
# PROCESAR UN GENOMA
# ============================================================================

process_genome() {
    local genome_dir="$1"
    local acc
    acc=$(basename "$genome_dir")
    local fasta
    fasta=$(find "$genome_dir" -name "*_genomic.fna" | head -1)

    if ! validate_fasta "$fasta"; then
        log_warn "${acc}: FASTA no encontrado, omitido"
        return 1
    fi

    local fasta_dir
    fasta_dir=$(dirname "$fasta")
    local fasta_name
    fasta_name=$(basename "$fasta")
    local species
    species=$(get_species "$acc")

    log_info "${acc} (${species})"

    run_amrfinderplus "$acc" "$fasta_dir" "$fasta_name" || return 1
    run_amrprofiler   "$acc" "$fasta" "$species"
    combine_beds      "$acc"
    mask_genome       "$acc" "$fasta_dir" "$fasta_name"
}

# ============================================================================
# MAIN
# ============================================================================

main() {
    echo ""
    echo "╔═══════════════════════════════════════════════════════════════════╗"
    echo "║  Limpieza AMR dual: AMRFinderPlus (Prokka+HMM) + AMRProfiler     ║"
    echo "║  v2.1 — lock de ejecución + tmpdir único en AMRFinderPlus         ║"
    echo "╚═══════════════════════════════════════════════════════════════════╝"
    echo ""

    [[ ! -d "$GENOMES_DIR" ]] && {
        log_error "GENOMES_DIR no encontrado: $GENOMES_DIR"
        exit 1
    }

    setup_directories
    acquire_lock

    echo "genome        amr_regions     amr_bp  timestamp" \
        > "$OUTPUT_DIR/stats/summary.txt"

    build_species_map

    local total=0 success=0 failed=0

    for genome_dir in "$GENOMES_DIR"/GCF_*/; do
        [[ ! -d "$genome_dir" ]] && continue
        total=$((total + 1))
        if process_genome "$genome_dir"; then
            success=$((success + 1))
        else
            failed=$((failed + 1))
        fi
    done

    echo ""
    echo "══════════════════════════════════════════════════════════════════════"
    log_info "Combinando ${success} genomas limpios..."
    cat "$OUTPUT_DIR"/masked/*_clean.fasta \
        > "$OUTPUT_DIR/metagenome_clean_combined.fasta"

    local n_seqs total_mb
    n_seqs=$(grep -c "^>" "$OUTPUT_DIR/metagenome_clean_combined.fasta" || echo 0)
    total_mb=$(( $(wc -c < "$OUTPUT_DIR/metagenome_clean_combined.fasta") / 1024 / 1024 ))

    echo ""
    echo "╔═══════════════════════════════════════════════════════════════════╗"
    echo "║  RESUMEN FINAL                                                    ║"
    echo "╚═══════════════════════════════════════════════════════════════════╝"
    echo ""
    echo "  Genomas procesados: ${success} / ${total}  (fallidos: ${failed})"
    echo "  Metagenoma:         ${n_seqs} secuencias, ~${total_mb} MB"
    echo ""
    echo "  📁 Salidas:"
    echo "     ├─ ${OUTPUT_DIR}/masked/                 ← genomas individuales"
    echo "     ├─ ${OUTPUT_DIR}/metagenome_clean_combined.fasta"
    echo "     ├─ ${OUTPUT_DIR}/tsv/                    ← resultados AMRFinderPlus"
    echo "     ├─ ${OUTPUT_DIR}/amrprofiler/            ← resultados AMRProfiler"
    echo "     ├─ ${OUTPUT_DIR}/prokka/                 ← anotaciones Prokka"
    echo "     └─ ${OUTPUT_DIR}/stats/summary.txt"
    echo ""
    echo "✅ Proceso completado"
}

# Guardar log con timestamp
LOG_FILE="$(dirname "$OUTPUT_DIR")/run_$(date +%Y%m%d_%H%M).log"
main "$@" 2>&1 | tee "$LOG_FILE"
