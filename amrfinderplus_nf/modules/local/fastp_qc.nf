/*
  (•˕ •マ.ᐟ
  Módulo: FASTP_QC
  ================================================================================
  Control de calidad, eliminación de adaptadores y trimming de lecturas crudas.
  ================================================================================
*/

process FASTP_QC {

    tag "${sample}"

    container params.img_fastp

    cpus params.threads

    publishDir "${params.outdir}/${sample}/02_qc", mode: 'copy'

    input:
    tuple val(sample), path(r1), path(r2)

    output:
    tuple val(sample), path("${sample}_1.trimmed.fastq.gz"), path("${sample}_2.trimmed.fastq.gz"), emit: trimmed_reads
    
    path "${sample}_fastp.html", emit: html_report
    path "${sample}_fastp.json", emit: json_report

    script:
    """
    echo "[INFO] Iniciando control de calidad (fastp) para: ${sample}"
    
    fastp \\
        -i "${r1}" \\
        -I "${r2}" \\
        -o "${sample}_1.trimmed.fastq.gz" \\
        -O "${sample}_2.trimmed.fastq.gz" \\
        -h "${sample}_fastp.html" \\
        -j "${sample}_fastp.json" \\
        --thread ${task.cpus} \\
        --detect_adapter_for_pe \\
        --qualified_quality_phred 20 \\
        --length_required 50
    """
}
