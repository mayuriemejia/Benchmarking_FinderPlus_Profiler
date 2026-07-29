/*
  ================================================================================
  Módulo: DOWNLOAD
  ================================================================================
  Si r1/r2 locales -> los usa directamente (sin descarga).
  Si vacíos/null   -> descarga desde SRA con fasterq-dump.
  ================================================================================
*/

process DOWNLOAD {

    tag "${sample_id}"

    publishDir "${params.outdir}/${sample_id}/01_raw_reads", mode: 'copy'

    input:
    tuple val(sample_id), val(r1), val(r2)

    output:
    tuple val(sample_id), path("${sample_id}_1.fastq.gz"), path("${sample_id}_2.fastq.gz"), emit: reads

    script:
    def use_local = (r1 != null && r1.toString() != 'null' && r1.toString().trim() != '' &&
                     r2 != null && r2.toString() != 'null' && r2.toString().trim() != '')

    if (use_local) {
        """
        echo "[INFO] Reads locales detectados para ${sample_id}, omitiendo descarga SRA."

        src_r1=\$(readlink -f "${r1}" 2>/dev/null || echo "${r1}")
        if [[ "\$src_r1" == *.gz ]]; then
            ln "\$src_r1" "${sample_id}_1.fastq.gz" 2>/dev/null || cp -L "\$src_r1" "${sample_id}_1.fastq.gz"
        else
            gzip -c "\$src_r1" > "${sample_id}_1.fastq.gz"
        fi

        src_r2=\$(readlink -f "${r2}" 2>/dev/null || echo "${r2}")
        if [[ "\$src_r2" == *.gz ]]; then
            ln "\$src_r2" "${sample_id}_2.fastq.gz" 2>/dev/null || cp -L "\$src_r2" "${sample_id}_2.fastq.gz"
        else
            gzip -c "\$src_r2" > "${sample_id}_2.fastq.gz"
        fi

        echo "[INFO] Reads listos: ${sample_id}_1.fastq.gz  ${sample_id}_2.fastq.gz"
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

        gzip "${sample_id}_1.fastq"
        gzip "${sample_id}_2.fastq"

        echo "[INFO] Descarga y compresión completadas para ${sample_id}."
        """
    }
}
