/*
  (•˕ •マ.ᐟ
  Módulo: AMR_EXCEL_REPORT
  ================================================================================
  Generación de reporte interactivo en formato Excel (openpyxl).
  ================================================================================
*/

process AMR_EXCEL_REPORT {

    tag "${sample}"

    container params.img_python

    publishDir "${params.outdir}/${sample}/06_report", mode: 'copy'

    input:
    tuple val(sample), path(results_tsv)
    path report_script

    output:
    tuple val(sample), path("${sample}_amr_report.xlsx"), emit: excel_report

    script:
    """
    echo "[INFO] Generando reporte Excel para: ${sample}"
    
    
    python "${report_script}" "${results_tsv}" "${sample}_amr_report.xlsx" "${sample}"
    """
}
