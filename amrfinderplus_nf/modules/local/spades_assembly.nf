/*
  (•˕ •マ.ᐟ
  Módulo: SPADES_ASSEMBLY (Modo Metagenómico)
  ================================================================================
  Ensamblaje de novo utilizando metaSPAdes.
  Adaptado para manejar coberturas y abundancias mixtas (5x, 10x, 20x, 30x).
  ================================================================================
*/

process SPADES_ASSEMBLY {

    tag "${sample}"

    container params.img_spades

    cpus params.threads
    memory "${params.mem_gb} GB"

    publishDir "${params.outdir}/${sample}/03_assembly", mode: 'copy'

    input:
    tuple val(sample), path(r1_trim), path(r2_trim)

    output:
    tuple val(sample), path("${sample}_contigs.fasta"), emit: assembly

    script:
    """
    echo "[INFO] Iniciando ensamblaje metaSPAdes para metagenoma: ${sample}"
    
    spades.py \\
        --meta \\
        -1 "${r1_trim}" \\
        -2 "${r2_trim}" \\
        -o ./spades_out \\
        --threads ${task.cpus} \\
        --memory ${params.mem_gb} \\
        -k 21,33,55,77

    cp ./spades_out/contigs.fasta ./${sample}_contigs.fasta
    """
}
