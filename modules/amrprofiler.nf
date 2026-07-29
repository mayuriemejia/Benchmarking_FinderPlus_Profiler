// (•˕ •マ.ᐟ
// Módulo: AMRPROFILER

process AMRPROFILER {
    tag "${sample}"
    container params.img_amrprofiler
    cpus params.threads
    publishDir "${params.outdir}/${sample}/05_amrprofiler", mode: 'copy'

    input:
    tuple val(sample), path(contigs), val(organism)
    path amrp_db

    output:
    tuple val(sample), path("${sample}_final_results_tool1.csv"),    optional: true, emit: genes
    tuple val(sample), path("${sample}_Core_mutations_results.csv"), optional: true, emit: mutations
    tuple val(sample), path("${sample}_rRNA_mutations_results.csv"), optional: true, emit: rrna
    tuple val(sample), path("${sample}_blast_results_*.csv"),        optional: true, emit: blast

    script:
    """
    chmod 777 .

    python "${amrp_db}/amrprofiler.py" \\
        "${contigs}" \\
        "${organism}" \\
        "${amrp_db}/" \\
        --threads            ${task.cpus} \\
        --identity_threshold 70 \\
        --coverage_threshold 70

    [ -f final_results_tool1.csv ]    && mv final_results_tool1.csv    "${sample}_final_results_tool1.csv"    || true
    [ -f Core_mutations_results.csv ] && mv Core_mutations_results.csv "${sample}_Core_mutations_results.csv" || true
    [ -f rRNA_mutations_results.csv ] && mv rRNA_mutations_results.csv "${sample}_rRNA_mutations_results.csv" || true
    [ -f blast_results_core.csv ]     && mv blast_results_core.csv     "${sample}_blast_results_core.csv"     || true
    [ -f blast_results_rRNA.csv ]     && mv blast_results_rRNA.csv     "${sample}_blast_results_rRNA.csv"     || true

    echo "[OK] AMRProfiler: ${sample}"
    """
}
