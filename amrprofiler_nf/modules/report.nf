/*
  ================================================================================
  Módulo: REPORT
  Mayurie Mejía | Trainee/Rookie | Sequentia Biotech
  ================================================================================
  Genera un reporte resumen en formato texto plano para cada muestra.
  Maneja outputs opcionales de AMRProfiler (genes, mutaciones, rRNA).
  ================================================================================
*/
process REPORT {
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/05_report", mode: 'copy'

    input:
    tuple val(sample_id), path(contigs), val(organism), path(genes_csv), path(mut_csv), path(rrna_csv)

    output:
    tuple val(sample_id), path("report_${sample_id}.txt"), emit: final_report

    script:
    """
    echo "[INFO] Generando reporte final para: ${sample_id}"

    NUM_CONTIGS=\$(grep -c '^>' "${contigs}" 2>/dev/null || echo "0")
    TOTAL_BP=\$(awk '!/^>/{t+=length(\$0)} END{print t}' "${contigs}" 2>/dev/null || echo "0")

    # Manejo seguro de archivos opcionales
    if [[ -f "${genes_csv}" && -s "${genes_csv}" ]]; then
        TOTAL_GENES=\$(tail -n +2 "${genes_csv}" | wc -l)
    else
        TOTAL_GENES="N/A"
    fi

    if [[ -f "${mut_csv}" && -s "${mut_csv}" ]]; then
        TOTAL_MUT=\$(tail -n +2 "${mut_csv}" | wc -l)
    else
        TOTAL_MUT="N/A"
    fi

    if [[ -f "${rrna_csv}" && -s "${rrna_csv}" ]]; then
        TOTAL_RRNA=\$(tail -n +2 "${rrna_csv}" | wc -l)
    else
        TOTAL_RRNA="N/A"
    fi

    cat << HEREDOC > "report_${sample_id}.txt"
============================================================
 RESUMEN FINAL — AMRProfiler | ${sample_id}
============================================================
 Organismo               : ${organism}
 Contigs ensamblados     : \$NUM_CONTIGS
 Tamaño total (bp)       : \$TOTAL_BP
 ── Resultados AMRProfiler ──
 Genes AMR adquiridos    : \$TOTAL_GENES
 Mutaciones genes core   : \$TOTAL_MUT
 Mutaciones rRNA         : \$TOTAL_RRNA
============================================================
HEREDOC

    echo "[INFO] Reporte generado: report_${sample_id}.txt"
    """
}
