#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { BENCHMARK_SIMULATION } from './modules/insilicoSeq.nf'

workflow {
    log.info """
    ╔══════════════════════════════════════════════════════════╗
    ║     monkey_d_benchmark — InSilicoSeq Simulation         ║
    ║     Sequentia Biotech                                    ║
    ╚══════════════════════════════════════════════════════════╝
    Base dir       : ${params.base_dir}
    Control FASTA  : ${params.control_fasta}
    Datasets       : ds1, ds2, ds3, ds4, control
    Sims/dataset   : ${params.n_sims}
    Seeds          : ${params.rep_seeds}
    Reads/sim      : ${params.n_reads}
    Modelo ISS     : ${params.model}
    GC bias        : ${params.gc_bias}
    Compress       : ${params.compress}
    CPUs           : ${params.cpus}
    Jobs totales   : ${5 * params.n_sims}
    """.stripIndent()

    BENCHMARK_SIMULATION()
}
