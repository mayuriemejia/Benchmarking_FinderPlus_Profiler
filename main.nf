#!/usr/bin/env nextflow
// (•˕ •マ.ᐟ
// AMR Benchmark — main.nf
// Mayurie Mejía | Sequentia Biotech
//
// Flujo:
//   reads → fastp → SPAdes → Prokka → AMRFinderPlus   (pipeline existente)
//                      ↓
//                   contigs ──────────→ AMRProfiler    (módulo propio)
//                                           ↓
//                   AMRFinder TSV + AMRProfiler CSV
//                                           ↓
//                        Evaluar vs truthset            (módulo propio)
//                                           ↓
//                   TP/FP/FN + Sensitivity/Precision/F1 + FP blast hits
//
// v1.1: se une también el canal `faa` (proteínas de Prokka) al join de
// EVAL_BENCHMARK — antes esa secuencia nunca llegaba al proceso y las
// columnas de secuencia en fp_amrfinder.tsv salían siempre "N/A".

nextflow.enable.dsl = 2

include { AMR_PIPELINE as AMRFINDER_PIPELINE } from "${params.amrfinder_nf_dir}/subworkflows/local/amr_pipeline"
include { AMRPROFILER    } from './modules/amrprofiler'
include { EVAL_BENCHMARK } from './modules/eval_benchmark'

workflow {

    if (!params.samplesheet) error "ERROR: especifica --samplesheet <ruta>"
    if (!params.truthset)    error "ERROR: especifica --truthset <ruta>"

    log.info """
    ╔══════════════════════════════════════════════════════════════╗
    ║   AMR Benchmark — AMRFinderPlus + AMRProfiler               ║
    ║   Sequentia Biotech                                          ║
    ╠══════════════════════════════════════════════════════════════╣
    ║  Samplesheet   : ${params.samplesheet}
    ║  Truthset      : ${params.truthset}
    ║  AMRFinder dir : ${params.amrfinder_nf_dir}
    ║  AMRProfiler DB: ${params.amrp_db}
    ║  Outdir        : ${params.outdir}
    ╚══════════════════════════════════════════════════════════════╝
    """.stripIndent()

    // ── Canal de muestras ─────────────────────────────────────────────────────
    ch_samplesheet = Channel.fromPath(params.samplesheet, checkIfExists: true)

    // Para pasar organism a AMRProfiler necesitamos leerlo aquí también
    organism_ch = ch_samplesheet
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def sid = row.sample_id?.trim()
            def org = row.organism?.trim() ?: 'Unknown'
            if (!sid) error "ERROR: sample_id vacío en samplesheet."
            tuple(sid, org)
        }

    // ════════════════════════════════════════════════════════════════════════
    // STEP 1: AMRFinderPlus pipeline completo (reutilizado sin modificar)
    //         reads → fastp → SPAdes → Prokka → AMRFinderPlus
    //         Escribe contigs en: outdir/sample/03_assembly/sample_contigs.fasta
    // ════════════════════════════════════════════════════════════════════════
    AMRFINDER_PIPELINE(ch_samplesheet)

    amrfinder_results = AMRFINDER_PIPELINE.out.results_tsv
    // tuple(sample, *_amr_results.tsv)

    faa_results = AMRFINDER_PIPELINE.out.faa
    // tuple(sample, *.faa)  — NUEVO, para EVAL_BENCHMARK

    // ════════════════════════════════════════════════════════════════════════
    // STEP 2: Recoger contigs del disco una vez que SPAdes ha terminado
    //         El subworkflow solo emite results_tsv, pero los contigs
    //         ya están escritos en publishDir por SPAdes
    // ════════════════════════════════════════════════════════════════════════
    contigs_ch = amrfinder_results
        .map { sid, tsv ->
            def contigs = file("${params.outdir}/${sid}/03_assembly/${sid}_contigs.fasta")
            tuple(sid, contigs)
        }
        .filter { sid, contigs ->
            if (!contigs.exists()) {
                log.warn "WARN: contigs no encontrados para ${sid} en ${contigs}"
                return false
            }
            return true
        }

    // ════════════════════════════════════════════════════════════════════════
    // STEP 3: AMRProfiler sobre los mismos contigs de SPAdes
    // ════════════════════════════════════════════════════════════════════════
    ch_amrp_db = Channel.value(file(params.amrp_db, checkIfExists: true))

    amrprofiler_in = contigs_ch
        .join(organism_ch, by: 0)
    // tuple(sample, contigs, organism)

    AMRPROFILER(amrprofiler_in, ch_amrp_db)

    amrprofiler_results = AMRPROFILER.out.genes
    // tuple(sample, *_final_results_tool1.csv)

    // ════════════════════════════════════════════════════════════════════════
    // STEP 4: Evaluación vs truthset
    //         → TP/FP/FN, Sensitivity, Precision, F1
    //         → FP con detalle de blast hit para auditoría
    // ════════════════════════════════════════════════════════════════════════
    ch_truthset = Channel.value(file(params.truthset, checkIfExists: true))
    ch_eval_script = Channel.value(file("${projectDir}/bin/eval_benchmark.py"))

    // Unir resultados de ambas herramientas + contigs + faa, por sample_id
    // remainder:true para no perder muestras donde AMRProfiler no detectó nada
    eval_in = amrfinder_results
        .join(amrprofiler_results, by: 0, remainder: true)
        .join(contigs_ch, by: 0)
        .join(faa_results, by: 0, remainder: true)
    // tuple(sample, amrfinder_tsv, amrprofiler_csv_or_null, contigs, faa_or_null)

    EVAL_BENCHMARK(eval_in, ch_truthset, ch_eval_script)
}
