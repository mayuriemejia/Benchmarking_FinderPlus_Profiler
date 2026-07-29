#!/usr/bin/env nextflow

/*
  (•˕ •マ.ᐟ
  AMRProfiler Pipeline — Salmonella enterica
  Mayurie Mejía | Trainee/Rookie | Sequentia Biotech
  ================================================================================
*/

nextflow.enable.dsl = 2

include { AMR_PIPELINE } from './workflows/amr_pipeline'

workflow {

    // 1. Validación: el samplesheet es obligatorio
    if (!params.samplesheet) {
        error "ERROR: Debes especificar un samplesheet con: --samplesheet ruta/al/samplesheet.tsv"
    }

    // 2. Registro en consola limpio y sin problemas de codificación
    log.info """
    ================================================================================
     AMRProfiler Pipeline
     Samplesheet : ${params.samplesheet}
     Especie     : ${params.species}
     Outdir      : ${params.outdir}
     DB          : ${params.amrp_db}
    ================================================================================
    """.stripIndent()

    // 3. Crear el canal a partir del Samplesheet
    ch_samplesheet = channel.fromPath(params.samplesheet, checkIfExists: true)

    // 4. Ejecutar el pipeline pasando el canal de entrada
    AMR_PIPELINE(ch_samplesheet)
}
