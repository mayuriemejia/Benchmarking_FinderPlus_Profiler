/*
  (•˕ •マ.ᐟ
  Subworkflow: AMR_PIPELINE
  ================================================================================
  Conecta módulos: SRA/Local -> QC -> Assembly -> Annotation -> AMR -> Report
  Acepta samplesheet.tsv con columnas: sample_id, r1, r2, organism

  NOTA: r1/r2 se pasan como val() al DOWNLOAD_SRA para permitir que sean
  nulos (descarga SRA) o rutas locales. DOWNLOAD_SRA normaliza la salida
  siempre a ${sample_id}_1.fastq.gz / ${sample_id}_2.fastq.gz.

  v1.2: -n usa el .fna de PROKKA (contigs renombrados con IDs compliant),
  NO los contigs crudos de SPAdes — evita el mismatch de IDs entre el .gff
  y el FASTA que AMRFinderPlus valida (gff_check). Por eso PROKKA_ANNOTATION.out.prokka_out
  ahora trae también el .fna, y ya no hace falta un join aparte con
  SPADES_ASSEMBLY.out.assembly para AMRFinderPlus.

  v1.1: ahora también emite `faa` (proteínas de Prokka), necesario para que
  EVAL_BENCHMARK pueda extraer la secuencia de los falsos positivos de
  AMRFinderPlus (antes no se pasaba como input real del proceso y esa
  columna siempre salía "N/A").
  ================================================================================
*/

include { DOWNLOAD_SRA }      from '../../modules/local/download_sra'
include { FASTP_QC }          from '../../modules/local/fastp_qc'
include { SPADES_ASSEMBLY }   from '../../modules/local/spades_assembly'
include { PROKKA_ANNOTATION } from '../../modules/local/prokka_annotation'
include { AMRFINDER_PLUS }    from '../../modules/local/amrfinder_plus'
include { AMR_EXCEL_REPORT }  from '../../modules/local/amr_excel_report'

workflow AMR_PIPELINE {

    take:
    ch_samplesheet

    main:

    // Ruta al script Python de reporte
    def report_script = file("${projectDir}/generate_amr_report.py", checkIfExists: true)

    // ── 1. Parseo del Samplesheet ────────────────────────────────────────────
    samples_ch = ch_samplesheet
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def sid = row.sample_id?.trim()
            def org = (row.organism && row.organism.trim() != '') ? row.organism.trim() : null
            def r1  = (row.r1 && row.r1.trim() != '')             ? row.r1.trim()       : null
            def r2  = (row.r2 && row.r2.trim() != '')             ? row.r2.trim()       : null

            if (!sid) error "ERROR: sample_id vacío en el samplesheet."

            tuple(sid, r1, r2, org)
        }

    // ── 2. Separación de canales ─────────────────────────────────────────────
    reads_ch    = samples_ch.map { sid, r1, r2, org -> tuple(sid, r1, r2) }
    organism_ch = samples_ch.map { sid, r1, r2, org -> tuple(sid, org)    }

    // ✻ 01 · Reads locales o descarga SRA
    DOWNLOAD_SRA(reads_ch)

    // ✻ 02 · Control de calidad (fastp)
    FASTP_QC(DOWNLOAD_SRA.out.reads)

    // ✻ 03 · Ensamblado de novo (metaSPAdes)
    SPADES_ASSEMBLY(FASTP_QC.out.trimmed_reads)

    // ✻ 04 · Anotación (Prokka, modo metagenoma) — emite faa, gff Y fna
    PROKKA_ANNOTATION(SPADES_ASSEMBLY.out.assembly)

    // ✻ 05 · AMRFinderPlus (modo combinado: -n con el .fna de Prokka + proteína + GFF)
    amr_in = PROKKA_ANNOTATION.out.prokka_out
        .join(organism_ch, by: 0)
    // tuple(sample, faa, gff, fna, organism)

    AMRFINDER_PLUS(amr_in)

    // ✻ 06 · Reporte Excel
    AMR_EXCEL_REPORT(AMRFINDER_PLUS.out.amr_report, report_script)

    // ── faa por separado, para EVAL_BENCHMARK (Paso 5) ────────────────────────
    faa_ch = PROKKA_ANNOTATION.out.prokka_out
        .map { sid, faa, gff, fna -> tuple(sid, faa) }

    emit:
    results_tsv = AMRFINDER_PLUS.out.amr_report
    faa         = faa_ch
}
