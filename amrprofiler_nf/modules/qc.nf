/*
  ================================================================================
  Módulo: QC (Quality Control)
  Mayurie Mejía | Trainee/Rookie | Sequentia Biotech
  ================================================================================
  Control de calidad, eliminación de adaptadores y trimming con fastp.
  ================================================================================
*/

process QC {
    
    tag "${sample_id}"
    
    publishDir "${params.outdir}/${sample_id}/02_qc", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}_1.trimmed.fastq.gz"), path("${sample_id}_2.trimmed.fastq.gz"), emit: reads
    path "${sample_id}_fastp.{html,json}", emit: report

    script:
    """
    echo "[INFO] Iniciando control de calidad (fastp) para: ${sample_id}"
    
    fastp \\
        -i "${r1}" -I "${r2}" \\
        -o "${sample_id}_1.trimmed.fastq.gz" \\
        -O "${sample_id}_2.trimmed.fastq.gz" \\
        -h "${sample_id}_fastp.html" \\
        -j "${sample_id}_fastp.json" \\
        --thread ${task.cpus} \\
        --detect_adapter_for_pe \\
        --qualified_quality_phred 20 \\
        --length_required 50
    """
}
