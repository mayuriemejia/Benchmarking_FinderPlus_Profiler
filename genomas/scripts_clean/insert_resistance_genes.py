#!/usr/bin/env python3
"""
insert_resistance_genes.py (v2)
PASO 2: Inserción de genes AMR en genomas limpios — diseño alineado con
        insilicoseq/modules/insilicoSeq.nf (seeds pareados por sim_idx)

Diseño:
  - Los 48 genomas limpios se reparten en 4 datasets de 12 "genomas base"
    cada uno, con composición de especies fija (quotas + relleno aleatorio),
    una sola vez, con --base-seed.
  - Cada dataset tiene un pool de genes AMR restringido (subconjunto de los
    40 disponibles), fijo, según la tabla de diseño experimental.
  - Para cada dataset se generan 10 simulaciones independientes (sim1..sim10),
    cada una con su propia semilla de inserción (11,22,...,110 — la MISMA
    semilla en la misma posición sim_idx para los 4 datasets, tal como espera
    insilicoSeq.nf). En cada sim, los 12 genomas base se re-cargan desde cero
    (limpios) y se les inserta una distribución nueva de genes de su pool.
  - El dataset "control" NO se genera aquí: ya existe como
    genomas/control_no_amr/metagenome_clean_combined.fasta (Paso 1, 48
    genomas completos, sin inserción).

Output (estructura esperada por insilicoSeq.nf):
  {output_dir}/ds1/results/sim1/ds1_sim1_metagenome.fasta
  {output_dir}/ds1/results/sim1/ds1_sim1_truthset.tsv
  ... (sim1..sim10, para ds1..ds4)
  {output_dir}/base_genome_assignment.tsv   (qué genoma cayó en qué dataset y por qué)

Uso:
  python3 insert_resistance_genes.py \
      --masked-dir genomas/control_no_amr/masked \
      --species-map genomas/control_no_amr/stats/species_map.tsv \
      --genes-dir genes_amr \
      --output-dir datasets \
      --base-seed 42
"""

import argparse
import random
import re
import sys
from pathlib import Path
from collections import defaultdict

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

# ============================================================================
# DISEÑO EXPERIMENTAL FIJO (según especificación del proyecto)
# ============================================================================

# sim_idx 1..10 -> semilla de inserción. La MISMA semilla en la misma
# posición se usa en los 4 datasets (pareado, igual que rep_seeds en ISS).
INSERTION_SEEDS = [11, 22, 33, 44, 55, 66, 77, 88, 99, 110]

DATASET_SPECS = {
    "ds1": {
        "quotas": {"Enterobacter": 5, "Pseudomonas": 3, "Acinetobacter": 2},
        "random_fill": 2,
        "genes": [
            "aac3II", "aac6Ib", "ant2Ia", "aph3Ia", "aph3Ib", "armA",
            "blaCTX-M-17", "blaGES-27", "blaIMP-8", "blaKPC-103", "blaNDM",
            "blaOXA-48", "blaPER-1", "blaSHV-12", "blaTEM-1", "blaVIM-20",
            "qnrA10", "qnrB", "qnrS", "rmtA", "rmtB", "sul1", "sul2",
        ],
    },
    "ds2": {
        "quotas": {"Klebsiella": 6, "Enterococcus": 3, "Staphylococcus": 2},
        "random_fill": 1,
        "genes": [
            "aac3II", "aac6Ib", "ant2Ia", "aph3Ia", "aph3Ib", "armA",
            "blaNDM", "blaOXA-48", "blaTEM-1", "dfrA1", "dfrA12", "ermB",
            "ermC", "mefA", "mphA", "mphE", "msrD", "rmtA", "rmtB",
            "sul1", "sul2",
        ],
    },
    "ds3": {
        "quotas": {"Acinetobacter": 6, "Pseudomonas": 3, "Enterobacter": 2},
        "random_fill": 1,
        "genes": [
            "aac3II", "aac6Ib", "blaKPC-103", "blaNDM", "blaOXA-48",
            "blaVIM-20", "dfrA12", "ermB", "ermC", "mcr-12", "mefA",
            "mphA", "mphE", "msrD", "qnrA10", "qnrS", "tetA", "tetB",
            "tetK", "tetL", "tetM", "tetO",
        ],
    },
    "ds4": {
        "quotas": {"Staphylococcus": 6, "Enterococcus": 3, "Klebsiella": 3},
        "random_fill": 0,
        "genes": [
            "armA", "blaIMP-8", "blaSHV-12", "blaTEM-1", "dfrA1", "dfrA12",
            "ermB", "ermC", "mcr-12", "mphA", "mphE", "qnrA10", "qnrB",
            "qnrS", "sul1", "sul2", "tetM", "tetO", "vanA", "vanB",
        ],
    },
}

GROUP_SIZE = 12  # genomas base por dataset (ds1-4)
CONTROL_SIZE = 12  # genomas para el dataset control (negativo, sin inserción)

# Correcciones manuales de género, aplicadas DESPUÉS de leer species_map.tsv.
# species_map.tsv se regenera en cada corrida de clean_genomes_dual.sh a partir
# de assembly_data_report.jsonl (taxonomía NCBI vigente), así que cualquier
# corrección debe vivir aquí para no perderse en la próxima regeneración.
#
#   GCF_001544195.1: NCBI lo reporta como "Tetragenococcus solitarius", pero
#   su biosample lo marca como "type strain of Enterococcus solitarius"
#   (reclasificación taxonómica posterior a su selección para este benchmark
#   ESKAPE). Se trata como Enterococcus para honrar la intención original.
SPECIES_OVERRIDES = {
    "GCF_001544195.1": "Enterococcus",
}
# NOTA: con 50 genomas totales, 4x12 (ds) + 12 (control) = 60 > 50, así que el
# control se elige de forma independiente y puede solapar con genomas ya
# usados en ds1-4 (misma secuencia limpia, sin ninguna inserción — sigue
# sirviendo como negativo real para ese genoma en concreto).


# ============================================================================
# FASTA helpers
# ============================================================================

def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def parse_fasta(path: Path) -> dict:
    contigs = {}
    header = None
    chunks = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    contigs[header] = "".join(chunks)
                header = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if header is not None:
            contigs[header] = "".join(chunks)
    return contigs


def write_fasta(fh, contig_id: str, seq: str, width: int = 70):
    fh.write(f">{contig_id}\n")
    for i in range(0, len(seq), width):
        fh.write(seq[i:i + width] + "\n")


def find_n_blocks(seq: str, min_len: int = 20):
    blocks = []
    for m in re.finditer(r"N+", seq, flags=re.IGNORECASE):
        length = m.end() - m.start()
        if length >= min_len:
            blocks.append((m.start(), m.end(), length))
    return blocks


def load_genes(genes_dir: Path) -> dict:
    """Devuelve {gene_name: {'seq':..., 'category':...}} para los 40 genes."""
    genes = {}
    for category_dir in sorted(genes_dir.iterdir()):
        if not category_dir.is_dir():
            continue
        for fasta_file in sorted(category_dir.glob("*.fasta")):
            recs = parse_fasta(fasta_file)
            if not recs:
                print(f"  [WARN] {fasta_file} vacío o ilegible, se omite", file=sys.stderr)
                continue
            gene_seq = next(iter(recs.values()))
            genes[fasta_file.stem] = {
                "name": fasta_file.stem,
                "category": category_dir.name,
                "seq": gene_seq.upper(),
            }
    return genes


# ============================================================================
# Reparto de los 48 genomas base en 4 datasets (quotas de especie + random)
# ============================================================================

def load_species_map(path: Path) -> dict:
    """accession -> género (primer token de la especie), con overrides aplicados"""
    species = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2 or not parts[0]:
                continue
            acc, sp = parts[0], parts[1]
            genus = sp.split()[0] if sp and sp != "Unknown" else "Unknown"
            species[acc] = genus

    for acc, genus_override in SPECIES_OVERRIDES.items():
        if acc in species and species[acc] != genus_override:
            print(f"  [INFO] Override de taxonomía: {acc} "
                  f"'{species[acc]}' -> '{genus_override}'", file=sys.stderr)
        species[acc] = genus_override

    return species


def assign_base_genomes(accessions: list, species_of: dict, seed: int) -> dict:
    rng = random.Random(seed)

    by_genus = defaultdict(list)
    for acc in accessions:
        by_genus[species_of.get(acc, "Unknown")].append(acc)
    for genus in by_genus:
        rng.shuffle(by_genus[genus])

    remaining = set(accessions)
    assignment = {ds: [] for ds in DATASET_SPECS}

    # 1. Cubrir quotas por género
    for ds_name, spec in DATASET_SPECS.items():
        for genus, count in spec["quotas"].items():
            pool = by_genus.get(genus, [])
            taken = 0
            while taken < count and pool:
                acc = pool.pop()
                if acc in remaining:
                    assignment[ds_name].append(acc)
                    remaining.discard(acc)
                    taken += 1
            if taken < count:
                print(f"  [WARN] {ds_name}: solo se encontraron {taken}/{count} "
                      f"genomas de género '{genus}'. Se rellenará con 'random'.",
                      file=sys.stderr)

    # 2. Relleno aleatorio (random_fill) + completar hasta GROUP_SIZE si faltó
    leftover = list(remaining)
    rng.shuffle(leftover)

    for ds_name, spec in DATASET_SPECS.items():
        need = GROUP_SIZE - len(assignment[ds_name])
        for _ in range(need):
            if not leftover:
                print(f"  [ERROR] No quedan genomas disponibles para completar "
                      f"{ds_name} (le faltan {need}).", file=sys.stderr)
                break
            assignment[ds_name].append(leftover.pop())

    for ds_name, accs in assignment.items():
        if len(accs) != GROUP_SIZE:
            print(f"  [WARN] {ds_name} terminó con {len(accs)} genomas "
                  f"(esperados {GROUP_SIZE})", file=sys.stderr)

    return assignment


# ============================================================================
# Inserción de genes (misma lógica robusta de coordenadas que v1)
# ============================================================================

def build_targets_for_genome(contigs: dict, n_genes: int, min_block_len: int):
    all_blocks = []
    for contig_id, seq in contigs.items():
        for start, end, length in find_n_blocks(seq, min_len=min_block_len):
            all_blocks.append({"contig": contig_id, "mode": "block",
                                "start": start, "end": end, "block_len": length})
    all_blocks.sort(key=lambda b: b["block_len"], reverse=True)

    targets = all_blocks[:n_genes]
    n_missing = n_genes - len(targets)
    if n_missing > 0:
        longest_contig = max(contigs.items(), key=lambda kv: len(kv[1]))[0]
        for _ in range(n_missing):
            targets.append({"contig": longest_contig, "mode": "append",
                             "start": None, "end": None, "block_len": None})
    return targets


def apply_insertions(contigs: dict, targets_with_genes: list):
    by_contig = defaultdict(list)
    for t in targets_with_genes:
        by_contig[t["contig"]].append(t)

    new_contigs = dict(contigs)
    truthset_records = []

    for contig_id, items in by_contig.items():
        seq = contigs[contig_id]
        block_items = sorted([i for i in items if i["mode"] == "block"],
                              key=lambda i: i["start"])
        append_items = [i for i in items if i["mode"] == "append"]

        pieces = []
        cursor = 0
        for item in block_items:
            gene_seq = item["gene"]["seq"]
            strand_seq = gene_seq if item["strand"] == "+" else revcomp(gene_seq)

            pieces.append(seq[cursor:item["start"]])
            insert_start_in_new = sum(len(p) for p in pieces)
            pieces.append(strand_seq)

            truthset_records.append({
                "contig": contig_id, "gene": item["gene"]["name"],
                "category": item["gene"]["category"], "mode": "block",
                "strand": item["strand"],
                "insert_start": insert_start_in_new,
                "insert_end": insert_start_in_new + len(strand_seq),
                "original_block_len": item["block_len"],
                "gene_len": len(strand_seq),
            })
            cursor = item["end"]
        pieces.append(seq[cursor:])
        new_seq = "".join(pieces)

        for item in append_items:
            gene_seq = item["gene"]["seq"]
            strand_seq = gene_seq if item["strand"] == "+" else revcomp(gene_seq)
            insert_start = len(new_seq)
            new_seq += strand_seq

            truthset_records.append({
                "contig": contig_id, "gene": item["gene"]["name"],
                "category": item["gene"]["category"], "mode": "append",
                "strand": item["strand"],
                "insert_start": insert_start,
                "insert_end": insert_start + len(strand_seq),
                "original_block_len": "NA", "gene_len": len(strand_seq),
            })

        new_contigs[contig_id] = new_seq

    return new_contigs, truthset_records


def select_control_genomes(accessions: list, size: int, seed: int) -> list:
    """Elige `size` genomas al azar para el control. Puede solapar con ds1-4."""
    rng = random.Random(seed)
    pool = accessions[:]
    if len(pool) < size:
        print(f"  [WARN] Solo hay {len(pool)} genomas disponibles, "
              f"menos de los {size} pedidos para control.", file=sys.stderr)
        return pool
    return rng.sample(pool, size)


def generate_control(control_genomes: list, masked_dir: Path, out_dir: Path):
    """Concatena los genomas de control SIN insertar nada (negativo puro)."""
    control_dir = out_dir / "control"
    control_dir.mkdir(parents=True, exist_ok=True)
    metagenome_path = control_dir / "control_metagenome.fasta"

    n_written = 0
    with open(metagenome_path, "w") as fasta_out:
        for acc in control_genomes:
            fasta_file = masked_dir / f"{acc}_clean.fasta"
            if not fasta_file.exists():
                print(f"  [WARN] No se encontró {fasta_file}, se omite {acc}", file=sys.stderr)
                continue
            contigs = parse_fasta(fasta_file)
            for contig_id, seq in contigs.items():
                write_fasta(fasta_out, contig_id, seq)
            n_written += 1

    return n_written, metagenome_path


# ============================================================================
# Generación de una simulación (una corrida de inserción) para un dataset
# ============================================================================

def generate_sim(ds_name: str, sim_idx: int, seed: int, base_genomes: list,
                  gene_pool: list, genes: dict, masked_dir: Path,
                  max_genes_per_genome: int, min_block_len: int, out_dir: Path):
    rng = random.Random(seed * 1000 + hash(ds_name) % 997)  # variación por dataset, reproducible

    sim_dir = out_dir / ds_name / "results" / f"sim{sim_idx}"
    sim_dir.mkdir(parents=True, exist_ok=True)

    metagenome_path = sim_dir / f"{ds_name}_sim{sim_idx}_metagenome.fasta"
    truthset_path = sim_dir / f"{ds_name}_sim{sim_idx}_truthset.tsv"

    n_total_inserted = 0
    n_genomes_with_genes = 0

    with open(metagenome_path, "w") as fasta_out, open(truthset_path, "w") as tsv_out:
        tsv_out.write("dataset\tsim\tgenome\tcontig\tgene\tcategory\tmode\tstrand\t"
                       "insert_start\tinsert_end\toriginal_block_len\tgene_len\n")

        for acc in base_genomes:
            fasta_file = masked_dir / f"{acc}_clean.fasta"
            if not fasta_file.exists():
                print(f"  [WARN] No se encontró {fasta_file}, se omite {acc}", file=sys.stderr)
                continue

            contigs = parse_fasta(fasta_file)  # siempre desde el genoma LIMPIO original
            n_genes = rng.randint(0, max_genes_per_genome)

            if n_genes > 0:
                targets = build_targets_for_genome(contigs, n_genes, min_block_len)
                targets_with_genes = []
                for t in targets:
                    gene_name = rng.choice(gene_pool)
                    gene = genes.get(gene_name)
                    if gene is None:
                        print(f"  [WARN] Gen '{gene_name}' del pool de {ds_name} no "
                              f"encontrado en genes_amr/, se omite esta inserción", file=sys.stderr)
                        continue
                    t2 = dict(t)
                    t2["gene"] = gene
                    t2["strand"] = rng.choice(["+", "-"])
                    targets_with_genes.append(t2)

                new_contigs, records = apply_insertions(contigs, targets_with_genes)
                n_total_inserted += len(records)
                if records:
                    n_genomes_with_genes += 1

                for r in records:
                    tsv_out.write(
                        f"{ds_name}\tsim{sim_idx}\t{acc}\t{r['contig']}\t{r['gene']}\t"
                        f"{r['category']}\t{r['mode']}\t{r['strand']}\t{r['insert_start']}\t"
                        f"{r['insert_end']}\t{r['original_block_len']}\t{r['gene_len']}\n"
                    )
            else:
                new_contigs = contigs

            for contig_id, seq in new_contigs.items():
                write_fasta(fasta_out, contig_id, seq)

    return n_genomes_with_genes, n_total_inserted


# ============================================================================
# MAIN
# ============================================================================

def main():
    ap = argparse.ArgumentParser(
        description="Inserta genes AMR (diseño por dataset, 10 sims pareadas)")
    ap.add_argument("--masked-dir", type=Path,
                     default=Path("genomas/control_no_amr/masked"))
    ap.add_argument("--species-map", type=Path,
                     default=Path("genomas/control_no_amr/stats/species_map.tsv"))
    ap.add_argument("--genes-dir", type=Path, default=Path("genes_amr"))
    ap.add_argument("--output-dir", type=Path, default=Path("datasets"))
    ap.add_argument("--base-seed", type=int, default=42,
                     help="Semilla para el reparto de los genomas base en 4 datasets (una sola vez)")
    ap.add_argument("--control-seed", type=int, default=999,
                     help="Semilla para elegir los genomas del dataset control (independiente, puede solapar con ds1-4)")
    ap.add_argument("--max-genes-per-genome", type=int, default=3)
    ap.add_argument("--min-block-len", type=int, default=20)
    args = ap.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Cargando genes AMR desde {args.genes_dir} ...")
    genes = load_genes(args.genes_dir)
    print(f"  {len(genes)} genes cargados")
    if len(genes) == 0:
        sys.exit("ERROR: no se cargó ningún gen, revisa --genes-dir")

    # Validar que todos los genes de cada pool existen realmente
    for ds_name, spec in DATASET_SPECS.items():
        missing = [g for g in spec["genes"] if g not in genes]
        if missing:
            print(f"  [WARN] {ds_name}: genes en el pool que no están en genes_amr/: "
                  f"{missing}", file=sys.stderr)

    print(f"Cargando mapa de especies desde {args.species_map} ...")
    species_of = load_species_map(args.species_map)
    print(f"  {len(species_of)} genomas con especie asignada")

    print(f"Listando genomas limpios en {args.masked_dir} ...")
    genome_files = sorted(args.masked_dir.glob("*_clean.fasta"))
    accessions = [f.name.replace("_clean.fasta", "") for f in genome_files]
    print(f"  {len(accessions)} genomas encontrados")

    print(f"Repartiendo 48 genomas en 4 datasets (base-seed={args.base_seed}) ...")
    base_assignment = assign_base_genomes(accessions, species_of, args.base_seed)

    with open(args.output_dir / "base_genome_assignment.tsv", "w") as f:
        f.write("dataset\tgenome\tgenus\n")
        for ds_name, accs in base_assignment.items():
            for acc in accs:
                f.write(f"{ds_name}\t{acc}\t{species_of.get(acc, 'Unknown')}\n")
    for ds_name, accs in base_assignment.items():
        print(f"  {ds_name}: {len(accs)} genomas base")

    print("\nGenerando simulaciones (10 por dataset) ...")
    for ds_name, spec in DATASET_SPECS.items():
        base_genomes = base_assignment[ds_name]
        gene_pool = spec["genes"]
        print(f"\n[{ds_name}] pool de {len(gene_pool)} genes, {len(base_genomes)} genomas base")

        for sim_idx, seed in enumerate(INSERTION_SEEDS, start=1):
            n_genomes_with_genes, n_inserted = generate_sim(
                ds_name=ds_name, sim_idx=sim_idx, seed=seed,
                base_genomes=base_genomes, gene_pool=gene_pool, genes=genes,
                masked_dir=args.masked_dir,
                max_genes_per_genome=args.max_genes_per_genome,
                min_block_len=args.min_block_len, out_dir=args.output_dir,
            )
            print(f"  sim{sim_idx:02d} (seed={seed}): {n_genomes_with_genes} genomas "
                  f"con AMR | {n_inserted} genes insertados")

    print(f"\n[control] seleccionando {CONTROL_SIZE} genomas "
          f"(control-seed={args.control_seed}, puede solapar con ds1-4) ...")
    control_genomes = select_control_genomes(accessions, CONTROL_SIZE, args.control_seed)
    n_written, control_path = generate_control(control_genomes, args.masked_dir, args.output_dir)
    print(f"  control: {n_written} genomas concatenados en {control_path} (sin inserción)")

    n_overlap = len(set(control_genomes) & set(
        acc for accs in base_assignment.values() for acc in accs
    ))
    if n_overlap > 0:
        print(f"  [INFO] {n_overlap}/{CONTROL_SIZE} genomas del control también "
              f"aparecen en ds1-4 (solape esperado y permitido, sin inserción en control)")

    with open(args.output_dir / "control_genome_list.tsv", "w") as f:
        f.write("genome\n")
        for acc in control_genomes:
            f.write(f"{acc}\n")

    print(f"\nListo. 40 metagenomas (4 datasets x 10 sims) + 1 metagenoma control "
          f"en {args.output_dir}/")
    print(f"IMPORTANTE: actualiza insilicoseq/nextflow.config -> control_fasta "
          f"para que apunte a {control_path} (12 genomas), no al "
          f"metagenome_clean_combined.fasta de 48/50 genomas del Paso 1.")


if __name__ == "__main__":
    main()
