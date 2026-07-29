#!/usr/bin/env python3
import csv, re, sys
from pathlib import Path

sample    = sys.argv[1]
truth_f   = sys.argv[2]
amrf_f    = sys.argv[3]
amrp_f    = sys.argv[4]
contigs_f = sys.argv[5] if len(sys.argv) > 5 else None
faa_f_arg = sys.argv[6] if len(sys.argv) > 6 else None  # .faa pasado directo

def normalize_gene(name):
    if not name:
        return ""
    return re.sub(r'[^a-z0-9]', '', name.lower().strip())

def extract_precise_gene(family_stem, product_name):
    """
    Reconstruye el alelo exacto de AMRProfiler usando la columna de producto,
    p.ej. blaTEM + 'beta-lactamase TEM-1' -> blaTEM-1

    AMRProfiler reporta el gen a nivel de familia en 'Resistance Gene'
    (p.ej. blaTEM, blaCTX-M, blaKPC) pero el alelo especifico suele
    aparecer en texto libre en 'Name of the product of the Resistance Gene'.
    Sin esto, comparar exacto contra el truthset (que usa alelos) genera
    FN + FP simultaneos para genes correctamente detectados.
    """
    if not product_name:
        return family_stem
    core = re.sub(r'^bla', '', family_stem, flags=re.IGNORECASE)
    core_esc = re.escape(core)
    pattern = rf'{core_esc}-?\d+[a-z]?\b'
    matches = re.findall(pattern, product_name, flags=re.IGNORECASE)
    if matches:
        prefix = family_stem[:len(family_stem) - len(core)] if core else family_stem
        return prefix + matches[-1]
    return family_stem

def match_genes(detected_set, expected_set):
    """Matching exacto. Ya no hace falta fuzzy matching (startswith) porque
    tanto AMRFinderPlus como AMRProfiler (via extract_precise_gene) reportan
    el alelo exacto igual que el truthset."""
    tp = len(detected_set & expected_set)
    fp = len(detected_set - expected_set)
    fn = len(expected_set - detected_set)
    return tp, fp, fn, detected_set - expected_set

def extract_from_fasta(fasta_path, seq_id, start=None, stop=None):
    if not fasta_path or not Path(fasta_path).exists():
        return "N/A"
    try:
        in_seq = False
        seq = []
        with open(fasta_path) as fh:
            for line in fh:
                line = line.strip()
                if line.startswith('>'):
                    header = line[1:].split()[0]
                    in_seq = (header == seq_id)
                elif in_seq:
                    seq.append(line)
        full_seq = ''.join(seq)
        if not full_seq:
            return "N/A"
        if start is not None and stop is not None:
            s = int(start) - 1
            e = int(stop)
            if s > e:
                s, e = e - 1, s + 1
            return full_seq[s:e]
        return full_seq
    except Exception as ex:
        return f"ERROR:{ex}"

# -- Extraer dataset Y sim del sample_id -------------------------------------
m = re.match(r'(ds\d+|control)(?:_sim(\d+))?', sample)
dataset = m.group(1) if m else sample
sim_idx = f"sim{m.group(2)}" if (m and m.group(2)) else None

# -- .faa de Prokka (para secuencias FP de AMRFinderPlus) --------------------
faa_f = faa_f_arg if (faa_f_arg and faa_f_arg != 'null' and Path(faa_f_arg).exists()) else None

# -- Truthset (control no tiene: sim_idx es None y no hay filas para 'control')
truth_genes = set()
if dataset != 'control':
    with open(truth_f) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            if row.get('dataset', '').strip() != dataset:
                continue
            if sim_idx is not None and row.get('sim', '').strip() != sim_idx:
                continue
            g = normalize_gene(row.get('gene', ''))
            if g:
                truth_genes.add(g)

# -- AMRFinderPlus - solo Type=AMR, excluir plus ------------------------------
amrf_genes = set()
amrf_rows  = []
if Path(amrf_f).exists() and Path(amrf_f).stat().st_size > 0:
    with open(amrf_f) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            if row.get('Type', '').strip().upper() != 'AMR':
                continue
            if row.get('Scope', '').strip().lower() == 'plus':
                continue
            g = normalize_gene(row.get('Element symbol', ''))
            if g:
                amrf_genes.add(g)
                amrf_rows.append(row)

# -- AMRProfiler ---------------------------------------------------------------
amrp_genes = set()
amrp_rows  = []
amrp_path  = Path(amrp_f) if amrp_f != 'null' else None
if amrp_path and amrp_path.exists() and amrp_path.stat().st_size > 0:
    with open(amrp_path) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            family  = row.get('Resistance Gene', '')
            product = row.get('Name of the product of the Resistance Gene', '')
            g_full  = extract_precise_gene(family, product)
            g = normalize_gene(g_full)
            if g:
                amrp_genes.add(g)
                amrp_rows.append(row)

# -- Metricas --------------------------------------------------------------------
def calc_metrics(tp, fp, fn):
    sens = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    prec = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    f1   = 2 * sens * prec / (sens + prec) if (sens + prec) > 0 else 0.0
    return sens, prec, f1

tp_f, fp_f, fn_f, fp_amrf = match_genes(amrf_genes, truth_genes)
tp_p, fp_p, fn_p, fp_amrp = match_genes(amrp_genes, truth_genes)
sens_f, prec_f, f1_f = calc_metrics(tp_f, fp_f, fn_f)
sens_p, prec_p, f1_p = calc_metrics(tp_p, fp_p, fn_p)

# -- FP AMRFinderPlus -> secuencia del FAA de Prokka ---------------------------
with open(f"{sample}_fp_amrfinder.tsv", 'w') as fh:
    fh.write("sample\tdataset\tgene_normalized\tElement symbol\t"
             "Protein id\tContig id\tStart\tStop\tClass\t%Identity\tsequence\n")
    for row in amrf_rows:
        g = normalize_gene(row.get('Element symbol', ''))
        if g in fp_amrf:
            protein_id = row.get('Protein id', '')
            seq = extract_from_fasta(faa_f, protein_id)
            fh.write(f"{sample}\t{dataset}\t{g}\t"
                     f"{row.get('Element symbol', '')}\t"
                     f"{protein_id}\t"
                     f"{row.get('Contig id', '')}\t"
                     f"{row.get('Start', '')}\t"
                     f"{row.get('Stop', '')}\t"
                     f"{row.get('Class', '')}\t"
                     f"{row.get('% Identity to reference', '')}\t"
                     f"{seq}\n")

# -- FP AMRProfiler -> secuencia del FASTA de contigs SPAdes -------------------
with open(f"{sample}_fp_amrprofiler.tsv", 'w') as fh:
    fh.write("sample\tdataset\tgene_normalized\tResistance Gene\t"
             "Contig ID\tQuery Start\tQuery Stop\t"
             "Antibiotic Class\tIdentity\tCoverage\tsequence\n")
    for row in amrp_rows:
        family  = row.get('Resistance Gene', '')
        product = row.get('Name of the product of the Resistance Gene', '')
        g = normalize_gene(extract_precise_gene(family, product))
        if g in fp_amrp:
            contig = row.get('Contig ID', '')
            start  = row.get('Query Start', '')
            stop   = row.get('Query Stop', '')
            seq = extract_from_fasta(contigs_f, contig, start, stop)
            fh.write(f"{sample}\t{dataset}\t{g}\t"
                     f"{row.get('Resistance Gene', '')}\t"
                     f"{contig}\t{start}\t{stop}\t"
                     f"{row.get('Antibiotic Class', '')}\t"
                     f"{row.get('Identity', '')}\t"
                     f"{row.get('Coverage', '')}\t"
                     f"{seq}\n")

# -- Metricas TSV ----------------------------------------------------------------
with open(f"{sample}_metrics.tsv", 'w') as fh:
    fh.write("sample\tdataset\tsim\ttool\tTP\tFP\tFN\tSensitivity\tPrecision\tF1\n")
    fh.write(f"{sample}\t{dataset}\t{sim_idx or 'NA'}\tAMRFinderPlus\t{tp_f}\t{fp_f}\t{fn_f}\t"
             f"{sens_f:.4f}\t{prec_f:.4f}\t{f1_f:.4f}\n")
    fh.write(f"{sample}\t{dataset}\t{sim_idx or 'NA'}\tAMRProfiler\t{tp_p}\t{fp_p}\t{fn_p}\t"
             f"{sens_p:.4f}\t{prec_p:.4f}\t{f1_p:.4f}\n")

# -- Summary -----------------------------------------------------------------------
with open(f"{sample}_summary.txt", 'w') as fh:
    fh.write(f"{'='*60}\n  Sample: {sample}\n  Dataset: {dataset}\n  Sim: {sim_idx or 'NA'}\n{'='*60}\n\n")
    fh.write(f"  Genes esperados ({len(truth_genes)}): {', '.join(sorted(truth_genes)) or 'ninguno'}\n\n")
    fh.write(f"  AMRFinderPlus (solo AMR adquirido):\n")
    fh.write(f"    TP={tp_f}  FP={fp_f}  FN={fn_f}  Sens={sens_f:.3f}  Prec={prec_f:.3f}  F1={f1_f:.3f}\n")
    if fp_amrf: fh.write(f"    FP: {', '.join(sorted(fp_amrf))}\n")
    fh.write(f"\n  AMRProfiler:\n")
    fh.write(f"    TP={tp_p}  FP={fp_p}  FN={fn_p}  Sens={sens_p:.3f}  Prec={prec_p:.3f}  F1={f1_p:.3f}\n")
    if fp_amrp: fh.write(f"    FP: {', '.join(sorted(fp_amrp))}\n")

print(f"[OK] {sample} ({dataset}/{sim_idx or 'NA'}): AMRFinder F1={f1_f:.3f} | AMRProfiler F1={f1_p:.3f}")
