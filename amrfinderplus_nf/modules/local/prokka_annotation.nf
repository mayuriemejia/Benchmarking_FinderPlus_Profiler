/*
  (•˕ •マ.ᐟ
  Módulo: PROKKA_ANNOTATION
  ================================================================================
  Anotación rápida de genomas procariotas utilizando Prokka.

  v1.1: ahora también emite el .fna (contigs renombrados por Prokka con IDs
  compliant tipo "gnl|X|PREFIX_N", por el flag --compliant). Es necesario
  para AMRFinderPlus en modo -n: los contigs originales de SPAdes NO
  coinciden con los IDs que Prokka escribe en el .gff (--compliant los
  renombra), así que pasarle a -n los contigs crudos de SPAdes causa
  "GFF contig id ... is not in the DNA FASTA file". Hay que usar el .fna
  de Prokka, que sí tiene los mismos IDs que el .gff.
  ================================================================================
*/

process PROKKA_ANNOTATION {

    tag "${sample}"

    container params.img_prokka

    cpus params.threads

    publishDir "${params.outdir}/${sample}/04_prokka", mode: 'copy'

    input:
    tuple val(sample), path(contigs)

    output:
    tuple val(sample), path("*.faa"), path("*.gff"), path("*.fna"), emit: prokka_out

    script:
    """
    echo "[INFO] Iniciando anotación Prokka (Modo Metagenoma) para: ${sample}"

    prokka \\
        --outdir . \\
        --prefix "${sample}" \\
        --metagenome \\
        --cpus ${task.cpus} \\
        --centre X \\
        --compliant \\
        --force \\
        "${contigs}"
    """
}
