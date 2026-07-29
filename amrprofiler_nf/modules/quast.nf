/*
  ================================================================================
  Módulo: QUAST (Quality Assessment of Assembly)
  Mayurie Mejía | Trainee/Rookie | Sequentia Biotech
  ================================================================================
  Evaluación métrica del ensamblaje (N50, longitud, contigs, etc.).
  ================================================================================
*/

process QUAST {
    
    tag "${sample}"
    
    publishDir "${params.outdir}/${sample}/03_assembly/quast_report", mode: 'copy'

    input:
    tuple val(sample), path(contigs)

    output:
    path "quast_report/", emit: report

    script:
    """
    echo "[INFO] Iniciando evaluación QUAST para: ${sample}"
    
    quast.py "${contigs}" \\
        -o quast_report \\
        --threads ${task.cpus} \\
        --min-contig 500
    """
}
