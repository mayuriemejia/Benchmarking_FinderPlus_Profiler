/*
  ================================================================================
  Módulo: AMRPROFILER
  ================================================================================
*/

process AMRPROFILER {
    tag "${sample}"
    container 'amrprofiler:1.0.0'

    publishDir "${params.outdir}/amrprofiler_results/${sample}", mode: 'copy'

    input:
    tuple val(sample), path(contigs), val(organism)
    path amrp_db

    output:
    tuple val(sample), path("${sample}_final_results_tool1.csv"),     optional: true, emit: genes
    tuple val(sample), path("${sample}_Core_mutations_results.csv"),  optional: true, emit: mutations
    tuple val(sample), path("${sample}_rRNA_mutations_results.csv"),  optional: true, emit: rrna
    tuple val(sample), path("${sample}_blast_results_*.csv"),         optional: true, emit: blast

    script:
    """
    echo "[INFO] Perfilando resistencia para: ${sample} (Organismo: ${organism})"

    python "${amrp_db}/amrprofiler.py" \\
        "${contigs}" \\
        "${organism}" \\
        "${amrp_db}/" \\
        --threads ${task.cpus ?: 2} \\
        --identity_threshold 70 \\
        --coverage_threshold 70

    # Renombrado seguro de outputs
    [ -f final_results_tool1.csv ]      && mv final_results_tool1.csv      "${sample}_final_results_tool1.csv"      || touch "${sample}_final_results_tool1.csv"
    [ -f Core_mutations_results.csv ]   && mv Core_mutations_results.csv   "${sample}_Core_mutations_results.csv"   || touch "${sample}_Core_mutations_results.csv"
    [ -f rRNA_mutations_results.csv ]   && mv rRNA_mutations_results.csv   "${sample}_rRNA_mutations_results.csv"   || touch "${sample}_rRNA_mutations_results.csv"
    [ -f blast_results_core.csv ]       && mv blast_results_core.csv       "${sample}_blast_results_core.csv"
    [ -f blast_results_rRNA.csv ]       && mv blast_results_rRNA.csv       "${sample}_blast_results_rRNA.csv"

    echo "[INFO] AMRProfiler finalizado exitosamente para ${sample}."
    """
}
