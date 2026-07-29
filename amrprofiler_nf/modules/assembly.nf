/*
  ================================================================================
  Módulo: ASSEMBLY (Modo Metagenómico)
  Mayurie Mejía | Trainee/Rookie | Sequentia Biotech
  ================================================================================
  Ensamblaje de novo con metaSPAdes, adaptado para comunidades mixtas.
  ================================================================================
*/

process ASSEMBLY {
    
    tag "${sample_id}"
    
    publishDir "${params.outdir}/${sample_id}/03_assembly", mode: 'copy'

    input:
    tuple val(sample_id), path(r1), path(r2)

    output:
    tuple val(sample_id), path("${sample_id}_contigs.fasta"), emit: contigs
    path "${sample_id}_assembly_stats.txt",                   emit: stats

    script:
    """
    echo "[INFO] Iniciando ensamblaje metaSPAdes para: ${sample_id}"
    
    spades.py \\
        --meta \\
        -1 "${r1}" -2 "${r2}" \\
        -o spades_out \\
        --threads ${task.cpus} \\
        --memory ${params.mem_gb} \\
        -k 21,33,55,77

    cp spades_out/contigs.fasta "${sample_id}_contigs.fasta"

    echo "Contigs: \$(grep -c '^>' "${sample_id}_contigs.fasta")"  > "${sample_id}_assembly_stats.txt"
    echo "Total bp: \$(awk '!/^>/{t+=length(\$0)} END{print t}' "${sample_id}_contigs.fasta")" >> "${sample_id}_assembly_stats.txt"
    """
}
