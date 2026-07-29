# Pipeline AMR · Nextflow DSL2
```
 (•˕ •マ.ᐟ  Mayurie Mejía · Trainee/Rookie · Sequentia Biotech
```

***

## ✦ Flujo del pipeline

```
SRA → fastp QC → SPAdes → Prokka → AMRFinderPlus → Reporte Excel
```

***

## ✦ Estructura (buenas prácticas · módulos)

```
.
├── main.nf
├── nextflow.config
├── params.schema.json
├── subworkflows/
│   └── local/
│       └── amr_pipeline.nf
├── modules/
│   └── local/
│       ├── download_sra.nf
│       ├── fastp_qc.nf
│       ├── spades_assembly.nf
│       ├── prokka_annotation.nf
│       ├── amrfinder_db_update.nf
│       ├── amrfinder_plus.nf
│       └── amr_excel_report.nf
└── generate_amr_report.py
```

***

## ✦ Requisitos

- **Nextflow** (DSL2)
- **Docker Desktop** (perfil `-profile docker`)
- Internet (descarga SRA + actualización DB de AMRFinderPlus)

> Tip: En Windows suele ser más estable correrlo desde **WSL2** con Docker Desktop.

***

## ✦ Cómo ejecutarlo

```bash
nextflow run main.nf -profile docker \
  --sample SRR27334358 \
  --outdir amr_pipeline \
  --threads 4 \
  --mem_gb 8
```

> ← Ajusta `--sample`, `--threads` y `--mem_gb` según tu muestra y recursos disponibles.

***

## ✦ Parámetros útiles

- `--sample`: accession SRA (ej. `SRR27334358`)
- `--outdir`: carpeta de salida (default `amr_pipeline`)
- `--threads`: CPUs (default `4`)
- `--mem_gb`: memoria para SPAdes en GB (default `8`)
- `--organism`: string para AMRFinderPlus (default `Salmonella`)
- `--ncbi_dir`: (opcional) ruta a tu carpeta `.ncbi` para SRA Toolkit

***

## ✦ Outputs generados

```
amr_pipeline/
├── 01_raw_reads/         → lecturas crudas  (*_1.fastq.gz, *_2.fastq.gz)
├── 01b_qc/               → lecturas trimadas + reporte fastp (HTML/JSON)
├── 02_assembly/          → ensamblado (contigs.fasta)
├── 03_annotation/        → anotación Prokka (*.faa, *.gff, ...)
├── 04_amrfinder/         → genes de resistencia (*_amr_results.tsv)
│   └── db/               → snapshot de la base de datos AMRFinderPlus usada ✶
├── 05_report/            → reporte final (*_amr_report.xlsx)
└── pipeline_info/        → trazabilidad (trace / report / timeline)
```

***

## ✦ Nota · Windows

Si usas Docker en Windows, la configuración más estable es **WSL2 + Docker Desktop**,
ejecutando Nextflow desde dentro de WSL2.

Si necesitas configurar SRA Toolkit, usa `--ncbi_dir` para montar
tu directorio local `.ncbi` en el contenedor como `/root/.ncbi`. ←

