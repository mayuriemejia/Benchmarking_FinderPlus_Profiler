#!/usr/bin/env nextflow

/*
    (•˕ •マ.ᐟ
    Ejercicio AMR | Salmonella
    Mayurie Mejía | Trainee/Rookie | Sequentia Biotech

    ================================================================================
    Pipeline AMR (Nextflow DSL2 · modular)
    Samplesheet -> [SRA/Local] -> fastp QC -> SPAdes -> Prokka -> AMRFinderPlus -> Excel
    ================================================================================
*/

nextflow.enable.dsl=2


include { AMR_PIPELINE } from './subworkflows/local/amr_pipeline'

workflow {

    // 1. Validación: el samplesheet es obligatorio
    if (!params.samplesheet) {
        error "ERROR: Debes especificar un samplesheet con: --samplesheet ruta/al/samplesheet.tsv"
    }

    ch_samplesheet = channel.fromPath(params.samplesheet, checkIfExists: true)

    AMR_PIPELINE( ch_samplesheet )
}
