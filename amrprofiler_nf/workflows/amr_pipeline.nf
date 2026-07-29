/*
  ================================================================================
  Subworkflow: AMR_PIPELINE (AMRProfiler)
  ================================================================================
*/

include { DOWNLOAD    } from '../modules/download'
include { QC          } from '../modules/qc'
include { ASSEMBLY    } from '../modules/assembly'
include { QUAST       } from '../modules/quast'
include { AMRPROFILER } from '../modules/amrprofiler'
include { REPORT      } from '../modules/report'

workflow AMR_PIPELINE {

    take:
    ch_samplesheet

    main:
    db_ch = Channel.value(params.amrp_db)

    samples_ch = ch_samplesheet
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def sid     = (row.sample_id ?: row.sample)?.trim()
            def org     = (row.organism && row.organism.trim() != '') ? row.organism.trim() : (params.species ?: 'Salmonella enterica')
            def r1      = (row.r1      && row.r1.trim()      != '') ? row.r1.trim()      : null
            def r2      = (row.r2      && row.r2.trim()      != '') ? row.r2.trim()      : null
            def contigs = (row.contigs && row.contigs.trim() != '') ? file(row.contigs.trim(), checkIfExists: true) : null

            if (!sid) error "ERROR: sample_id vacío en el samplesheet."

            tuple(sid, r1, r2, org, contigs)
        }

    organism_ch = samples_ch.map { sid, r1, r2, org, ctg -> tuple(sid, org) }

    // ── Rama 1: contigs ya disponibles ───────────────────
    preassembled_ch = samples_ch
        .filter { sid, r1, r2, org, ctg -> ctg != null }
        .map    { sid, r1, r2, org, ctg -> tuple(sid, ctg, org) }

    // ── Rama 2: reads locales -> QC -> Assembly ───────────
    local_reads_ch = samples_ch
        .filter { sid, r1, r2, org, ctg -> ctg == null && r1 != null }
        .map    { sid, r1, r2, org, ctg -> tuple(sid, r1, r2) }

    // ── Rama 3: SRA -> Download -> QC -> Assembly ─────────
    to_download_ch = samples_ch
        .filter { sid, r1, r2, org, ctg -> ctg == null && r1 == null }
        .map    { sid, r1, r2, org, ctg -> tuple(sid, r1, r2) }

    downloaded_ch   = DOWNLOAD(to_download_ch)
    reads_for_qc_ch = local_reads_ch.mix(downloaded_ch)
    qc_ch           = QC(reads_for_qc_ch)
    assembly_ch     = ASSEMBLY(qc_ch.reads)

    assembled_ch = assembly_ch.contigs.join(organism_ch, by: 0)

    // ── Une ambas ramas ───────────────────────────────────
    contigs_ch = preassembled_ch.mix(assembled_ch)

    // ✻ QUAST
    QUAST(contigs_ch.map { sid, ctg, org -> tuple(sid, ctg) })

    // ✻ AMRProfiler
    amrp_out = AMRPROFILER(contigs_ch, db_ch)

    // ✻ Reporte — join por sample_id para alinear todos los outputs
    report_ch = contigs_ch
        .join(amrp_out.genes,     by: 0)
        .join(amrp_out.mutations, by: 0)
        .join(amrp_out.rrna,      by: 0)

    REPORT(report_ch)
}
