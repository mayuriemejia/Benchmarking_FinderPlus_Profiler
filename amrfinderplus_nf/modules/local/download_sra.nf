/*
  (•˕ •マ.ᐟ
  Módulo: DOWNLOAD_SRA
  ================================================================================
  Si los FASTQ ya existen en local  -> los usa directamente (sin descarga).
  Si no existen (r1/r2 vacíos/null) -> los descarga desde SRA con fasterq-dump.
  ================================================================================
*/

process DOWNLOAD_SRA {

    tag "${sample_id}"

    container params.img_sra

    cpus params.threads

    publishDir "${params.outdir}/${sample_id}/01_raw_reads", mode: 'copy'

    input:
    tuple val(sample_id), val(r1), val(r2)

    output:
    tuple val(sample_id), path("${sample_id}_1.fastq.gz"), path("${sample_id}_2.fastq.gz"), emit: reads

    script:
    def use_local = (r1 != null && r1.toString() != 'null' && r1.toString() != '' &&
                     r2 != null && r2.toString() != 'null' && r2.toString() != '')

    if (use_local) {
        """
        echo "[INFO] Reads locales detectados para ${sample_id}, omitiendo descarga SRA."

        # ── R1 ──────────────────────────────────────────────────────────────────
        src_r1=\$(readlink -f "${r1}" 2>/dev/null || echo "${r1}")
        dst_r1="${sample_id}_1.fastq.gz"

        if [[ "\$src_r1" == *.gz ]]; then
            if ln "\$src_r1" "\$dst_r1" 2>/dev/null; then
                echo "[INFO] Hard-link creado para R1."
            else
                cp -L "\$src_r1" "\$dst_r1"
                echo "[INFO] Copia realizada para R1."
            fi
        else
            echo "[INFO] R1 sin comprimir, comprimiendo..."
            gzip -c "\$src_r1" > "\$dst_r1"
        fi

        # ── R2 ──────────────────────────────────────────────────────────────────
        src_r2=\$(readlink -f "${r2}" 2>/dev/null || echo "${r2}")
        dst_r2="${sample_id}_2.fastq.gz"

        if [[ "\$src_r2" == *.gz ]]; then
            if ln "\$src_r2" "\$dst_r2" 2>/dev/null; then
                echo "[INFO] Hard-link creado para R2."
            else
                cp -L "\$src_r2" "\$dst_r2"
                echo "[INFO] Copia realizada para R2."
            fi
        else
            echo "[INFO] R2 sin comprimir, comprimiendo..."
            gzip -c "\$src_r2" > "\$dst_r2"
        fi

        echo "[INFO] Reads listos: \$dst_r1  \$dst_r2"
        """
    } else {
        """
        echo "[INFO] No se encontraron reads locales para ${sample_id}. Descargando desde SRA..."

        fasterq-dump ${sample_id} \\
            --split-files \\
            --threads ${task.cpus} \\
            --outdir .

        if [[ ! -f "${sample_id}_1.fastq" || ! -f "${sample_id}_2.fastq" ]]; then
            echo "ERROR: fasterq-dump no generó los archivos esperados para ${sample_id}" >&2
            exit 1
        fi

        echo "[INFO] Comprimiendo reads descargados..."
        gzip "${sample_id}_1.fastq"
        gzip "${sample_id}_2.fastq"

        echo "[INFO] Descarga y compresión completadas para ${sample_id}."
        """
    }
}
