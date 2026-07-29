/*
  (•˕ •マ.ᐟ
  Módulo: AMRFINDER_PLUS
  ================================================================================
  Detección de genes de resistencia antimicrobiana (AMR) y factores de virulencia.

  v1.2: modo combinado COMPLETO — nucleótido (-n) + proteína (--protein) +
  GFF (--annotation_format prokka) + HMM (implícito con --plus).

  IMPORTANTE: -n usa el .fna de PROKKA (contigs renombrados, IDs compliant
  tipo "gnl|X|PREFIX_N"), NO los contigs crudos de SPAdes. AMRFinderPlus
  valida (gff_check) que los IDs de contig del .gff existan en el FASTA
  de -n; como Prokka corre con --compliant, sus IDs de salida no coinciden
  con los headers originales de SPAdes ("NODE_N_length_..."), así que pasar
  los contigs de SPAdes directamente rompe con "GFF contig id ... is not
  in the DNA FASTA file". v1.1 tenía este bug — corregido en v1.2.
  ================================================================================
*/
process AMRFINDER_PLUS {
    tag "${sample}"
    container params.img_amr
    cpus params.threads
    publishDir "${params.outdir}/${sample}/05_amrfinder", mode: 'copy'
    input:
    tuple val(sample), path(faa), path(gff), path(prokka_fna), val(organism)
    output:
    tuple val(sample), path("${sample}_amr_results.tsv"), emit: amr_report
    script:
    def org_flag = (organism && organism.trim() != '' && organism.toLowerCase() != 'null' && organism.toLowerCase() != 'metagenome') ? "--organism \"${organism}\"" : ""
    """
    echo "[INFO] Iniciando AMRFinderPlus para: ${sample}"
    if [ -n "${org_flag}" ]; then
        echo "[INFO] Modo específico de especie activado: ${organism}"
    else
        echo "[INFO] Modo general (metagenómico) activado."
    fi
    echo "[INFO] Modo combinado: nucleótido (Prokka .fna) + proteína + GFF + HMM"
    amrfinder \\
        -n                      "${prokka_fna}" \\
        --protein               "${faa}" \\
        --gff                   "${gff}" \\
        --annotation_format     prokka \\
        ${org_flag} \\
        --plus \\
        --output                "${sample}_amr_results.tsv" \\
        --threads               ${task.cpus}
    """
}
