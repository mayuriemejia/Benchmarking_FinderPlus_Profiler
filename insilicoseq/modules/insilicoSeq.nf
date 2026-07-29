// (•˕ •マ.ᐟ
// Module:   SIMULATE_READS_BENCHMARK
// Workflow: BENCHMARK_SIMULATION
//
// Diseño experimental completo:
//
//   UPSTREAM (insert_resistance_genes.py, ya ejecutado):
//     --seeds 11 22 33 44 55 66 77 88 99 110  →  sim1..sim10 por dataset
//     ds{1-4}/results/sim{K}/ds{N}_sim{K}_metagenome.fasta  (40 FASTAs)
//     control: sin inserción, usa masked FASTA directamente
//
//   ESTE PIPELINE (ISS reads generation):
//     1 ISS run por metagenoma, seed ISS pareado por sim_idx
//     Inputs  (auto-descubiertos): ds{1-4}/results/sim{1-10}/*_metagenome.fasta
//     Control: 10 runs con los mismos seeds sobre neg_ctrl_masked.fna
//     Total jobs ISS: 4 datasets × 10 sims + 1 control × 10 = 50 jobs
//
//   DOWNSTREAM (AMR detection):
//     50 outputs × 2 herramientas (AMRFinderPlus + AMRProfiler) = 100 runs
//     → 10 valores (Sensitivity, Specificity, F1) por dataset
//     → Wilcoxon pareado entre datasets (sim_idx como bloque)
//
//   Invariante de seeds ISS:
//     sim_idx K usa iss_seeds[K-1] en TODOS los datasets.
//     Seeds ISS distintos de los de inserción para evitar correlación.
//
//     sim_idx  │  1    2    3    4    5    6    7    8    9   10
//     ins seed │  11   22   33   44   55   66   77   88   99  110  ← insert_resistance_genes.py
//     ISS seed │ 200  400  600  800 1000 1200 1400 1600 1800 2000  ← este pipeline
//               ↑ idéntico en ds1, ds2, ds3, ds4 y control por columna

// ─── Params defaults ──────────────────────────────────────────────────────────
params.base_dir      = "/media/sequentia/isilon/students/visitor6/monkey_d_benchmark/datasets"
params.cpus          = 8
params.gc_bias       = true
params.compress      = true
params.n_reads = "6500000"
params.model         = "hiseq"
params.n_sims        = 10   // sim1..sim10

// Seeds ISS pareados por sim_idx — distintos de los seeds de inserción
params.rep_seeds = [200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000]

// ─── Process ──────────────────────────────────────────────────────────────────
process SIMULATE_READS_BENCHMARK {
    tag "ISS | ${dataset_id} | sim${sprintf('%02d', sim_idx)} | seed=${seed}"

    // Output dentro de la estructura existente del proyecto
    publishDir "${params.base_dir}/${dataset_id}/reads/sim${sim_idx}", mode: 'copy'

    container 'myinsilico:2.0'
    containerOptions '--entrypoint ""'

    cpus params.cpus

    input:
    tuple val(dataset_id), val(sim_idx), val(seed), path(metagenome_fasta)

    output:
    tuple val(dataset_id), val(sim_idx), val(seed),
          path("${dataset_id}_sim${sim_idx}_R1.fastq${params.compress ? '.gz' : ''}"),
          path("${dataset_id}_sim${sim_idx}_R2.fastq${params.compress ? '.gz' : ''}"),
          emit: reads

    script:
    def prefix        = "${dataset_id}_sim${sim_idx}"
    def gc_bias_flag  = params.gc_bias  ? '--gc_bias'  : ''
    def compress_flag = params.compress ? '--compress' : ''
    """
    iss generate \\
        --genomes ${metagenome_fasta} \\
        --model   ${params.model}     \\
        --n_reads ${params.n_reads}   \\
        --output  ${prefix}           \\
        --cpus    ${task.cpus}        \\
        --seed    ${seed}             \\
        ${gc_bias_flag}               \\
        ${compress_flag}

    echo "✓ ${dataset_id} | sim${sim_idx} | seed=${seed} completado"
    ls -lh ${prefix}*
    """
}

// ─── Workflow ─────────────────────────────────────────────────────────────────
workflow BENCHMARK_SIMULATION {
    main:

    // ── Canal experimental: auto-descubrir FASTAs de ds1–ds4 ─────────────────
    // Patrón: {base_dir}/{ds1..ds4}/results/sim{N}/{ds}_sim{N}_metagenome.fasta
    experimental_ch = Channel
        .fromPath("${params.base_dir}/{ds1,ds2,ds3,ds4}/results/sim*/*_metagenome.fasta")
        .map { fasta ->
            // Extraer dataset_id y sim_idx del path
            // .../ds1/results/sim3/ds1_sim3_metagenome.fasta
            //        ^^^           ^
            def path_str  = fasta.toString()
            def dataset   = (path_str =~ /\/(ds\d+)\/results/)[0][1]
            def sim_idx   = (path_str =~ /\/sim(\d+)\//)[0][1].toInteger()
            def seed      = params.rep_seeds[sim_idx - 1]
            tuple(dataset, sim_idx, seed, fasta)
        }

    // ── Canal control: 5 runs con los mismos seeds sobre el mismo FASTA ───────
    // El control no tiene inserción → el mismo masked FASTA para todas las sims
    control_ch = Channel
        .from(1..params.n_sims)
        .map { sim_idx ->
            def seed = params.rep_seeds[sim_idx - 1]
            tuple("control", sim_idx, seed, file(params.control_fasta, checkIfExists: true))
        }

    // ── Combinar y lanzar ─────────────────────────────────────────────────────
    // 4 datasets × 10 sims + 1 control × 10 runs = 50 jobs ISS
    all_ch = experimental_ch.mix(control_ch)

    SIMULATE_READS_BENCHMARK(all_ch)

    emit:
    reads = SIMULATE_READS_BENCHMARK.out.reads
}
