# AMRProfiler Pipeline — Nextflow
<!-- (•˕ •マ.ᐟ -->
<!-- Mayurie Mejía | Trainee/Rookie | Sequentia Biotech -->

Pipeline modular para detección de resistencia antimicrobiana. ✶
Diseñado para ser reproducible, limpio y fácil de seguir.

`SRA → fastp QC → SPAdes → QUAST → AMRProfiler → Reporte`

## Requisitos

- Nextflow ≥ 23.x
- Docker
- Imagen `amrprofiler:1.0.0` (ver abajo)
- Bases de datos AMRProfiler descargadas desde [Zenodo](https://zenodo.org/records/15674467)

## Setup (solo la primera vez)

```bash
# 1. Construir imagen Docker de AMRProfiler
docker build -t amrprofiler:1.0.0 docker/amrprofiler/

# 2. Descargar y descomprimir bases de datos
wget -O amrprofiler_db/amrprofiler_databases.zip \
    "https://zenodo.org/records/15674467/files/amrprofiler_databases.zip?download=1"
unzip amrprofiler_db/amrprofiler_databases.zip -d amrprofiler_db/
```

## Uso

```bash
# Correr con valores por defecto (SRR27334358 / Salmonella enterica)
nextflow run main.nf -profile docker

# Cambiar muestra o especie
nextflow run main.nf -profile docker \
    --sample SRR12345678 \
    --species "Escherichia coli"

# Reanudar desde caché
nextflow run main.nf -profile docker -resume
```

## Outputs

| Directorio            | Contenido                              |
|-----------------------|----------------------------------------|
| `results/01_raw_reads`| Lecturas crudas `.fastq.gz`            |
| `results/02_qc`       | Lecturas trimmed + reportes fastp      |
| `results/03_assembly` | `contigs.fasta` + reporte QUAST        |
| `results/04_amrprofiler` | CSVs genes AMR y mutaciones         |
| `results/05_report`   | Resumen final `.txt`                   |
| `results/trace`       | Timeline, DAG, trace y report HTML     |

## Trazabilidad

Generada automáticamente en `results/trace/`:
- `timeline.html` — línea de tiempo de procesos  ←
- `report.html`   — reporte de ejecución
- `trace.txt`     — log detallado con CPU/memoria por tarea
- `dag.svg`       — grafo del pipeline