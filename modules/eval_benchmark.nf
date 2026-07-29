process EVAL_BENCHMARK {
    tag "${sample}"
    container params.img_python
    publishDir "${params.outdir}/${sample}/06_benchmark", mode: 'copy'
    input:
    tuple val(sample),
          path(amrfinder_tsv),
          path(amrprofiler_csv),
          path(contigs),
          path(faa)
    path truthset
    path eval_script
    output:
    tuple val(sample), path("${sample}_metrics.tsv"),        emit: metrics
    tuple val(sample), path("${sample}_fp_amrfinder.tsv"),   emit: fp_amrfinder
    tuple val(sample), path("${sample}_fp_amrprofiler.tsv"), emit: fp_amrprofiler
    tuple val(sample), path("${sample}_summary.txt"),        emit: summary
    script:
    def amrp = amrprofiler_csv ? amrprofiler_csv : 'null'
    def faa_arg = faa ? faa : 'null'
    """
    python3 ${eval_script} \
        "${sample}" \
        "${truthset}" \
        "${amrfinder_tsv}" \
        "${amrp}" \
        "${contigs}" \
        "${faa_arg}"
    """
}
