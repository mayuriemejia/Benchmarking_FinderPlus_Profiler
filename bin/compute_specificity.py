#!/usr/bin/env python3
"""
Calcula Specificity para AMRFinderPlus y AMRProfiler, usando como universo
de "verdaderos negativos" el PANEL DE GENES POR DATASET: la union de todos
los genes que fueron insertados en CUALQUIER simulacion de ese dataset
(ej. panel(ds1) = todos los genes distintos que aparecen en filas ds1 de
truthset_master.tsv, sin importar el 'sim').

Para una muestra concreta (ej. ds1_sim3):
  - Clase positiva = genes insertados en ESA muestra especifica (ya usado
    para TP/FN en eval_benchmark.py).
  - Clase negativa = Panel(ds1) - genes insertados en esta muestra. Son
    genes que el diseno experimental si contempla para ds1, pero que en
    esta simulacion concreta no fueron insertados.
  - Un gen de la clase negativa es FP si alguna deteccion lo matchea
    EXACTO (mismo criterio que match_genes() en eval_benchmark.py tras el
    fix: sin fuzzy matching por prefijo). Si ninguna deteccion lo matchea
    exacto, es TN.
  - Specificity = TN / (TN + FP)

AMRProfiler reporta el gen a nivel de familia ('Resistance Gene', p.ej.
blaTEM) en vez del alelo especifico. Se reconstruye el alelo exacto via
extract_precise_gene() usando la columna 'Name of the product of the
Resistance Gene' (misma logica que eval_benchmark.py), para poder
comparar exacto contra el truthset sin perder ni inventar especificidad.

Detecciones que no matchean NINGUN gen del panel del dataset (ni en la
clase positiva ni en la negativa) quedan FUERA del calculo de Specificity
-- no pertenecen al espacio controlado del diseno experimental (podrian
ser resistoma nativo del genoma base, ruido, etc.) -- y se reportan aparte
como diagnostico, sin afectar el resultado.

Uso:
    python3 compute_specificity.py \
        --truthset truthset_master.tsv \
        --benchmark-dir benchmark_results \
        --out specificity_metrics.tsv
"""
import argparse
import csv
import re
import sys
from pathlib import Path
from collections import defaultdict


def normalize_gene(name):
    if not name:
        return ""
    return re.sub(r'[^a-z0-9]', '', name.lower().strip())


def extract_precise_gene(family_stem, product_name):
    """Reconstruye el alelo exacto de AMRProfiler usando la columna de producto,
    p.ej. blaTEM + 'beta-lactamase TEM-1' -> blaTEM-1. Misma logica que
    eval_benchmark.py -- mantenerlas sincronizadas."""
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


def gene_matches(detected, candidate):
    """Matching exacto (post-fix). Ya no hay fuzzy startswith -- tanto
    AMRFinderPlus como AMRProfiler (via extract_precise_gene) reportan el
    alelo exacto igual que el truthset."""
    return detected == candidate


def build_dataset_panels(truth_f):
    """
    panel[dataset] = set de genes normalizados insertados en CUALQUIER
    simulacion de ese dataset.
    """
    panels = defaultdict(set)
    with open(truth_f) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            dataset = row.get('dataset', '').strip()
            if not dataset or dataset == 'control':
                continue
            g = normalize_gene(row.get('gene', ''))
            if g:
                panels[dataset].add(g)
    return panels


def load_sample_truth(truth_f, dataset, sim_idx):
    genes = set()
    with open(truth_f) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            if row.get('dataset', '').strip() != dataset:
                continue
            if sim_idx is not None and row.get('sim', '').strip() != sim_idx:
                continue
            g = normalize_gene(row.get('gene', ''))
            if g:
                genes.add(g)
    return genes


def load_amrfinder_detected(amrf_f):
    detected = set()
    p = Path(amrf_f)
    if not p.exists() or p.stat().st_size == 0:
        return detected
    with open(p) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            if row.get('Type', '').strip().upper() != 'AMR':
                continue
            if row.get('Scope', '').strip().lower() == 'plus':
                continue
            g = normalize_gene(row.get('Element symbol', ''))
            if g:
                detected.add(g)
    return detected


def load_amrprofiler_detected(amrp_f):
    detected = set()
    if amrp_f is None:
        return detected
    p = Path(amrp_f)
    if not p.exists() or p.stat().st_size == 0:
        return detected
    with open(p) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            family = row.get('Resistance Gene', '')
            product = row.get('Name of the product of the Resistance Gene', '')
            g_full = extract_precise_gene(family, product)
            g = normalize_gene(g_full)
            if g:
                detected.add(g)
    return detected


def compute_specificity_for_sample(panel, truth_genes, detected_genes):
    negative_class = panel - truth_genes

    fp_panel_genes = set()
    for e in negative_class:
        if any(gene_matches(d, e) for d in detected_genes):
            fp_panel_genes.add(e)

    tn = len(negative_class) - len(fp_panel_genes)
    fp = len(fp_panel_genes)
    specificity = tn / (tn + fp) if (tn + fp) > 0 else float('nan')

    out_of_panel = [d for d in detected_genes
                    if not any(gene_matches(d, g) for g in panel)]

    return {
        'panel_size': len(panel),
        'truth_size': len(truth_genes),
        'neg_class_size': len(negative_class),
        'fp_panel': fp,
        'tn': tn,
        'specificity': specificity,
        'out_of_panel_count': len(out_of_panel),
        'out_of_panel_genes': ",".join(sorted(out_of_panel)) if out_of_panel else ""
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--truthset', required=True)
    ap.add_argument('--benchmark-dir', required=True)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    panels = build_dataset_panels(args.truthset)
    print(f"[info] Paneles por dataset:", file=sys.stderr)
    for ds, genes in sorted(panels.items()):
        print(f"    {ds}: {len(genes)} genes -> {sorted(genes)}", file=sys.stderr)

    bdir = Path(args.benchmark_dir)
    amrf_files = sorted(bdir.glob("*/05_amrfinder/*_amr_results.tsv"))
    if not amrf_files:
        print(f"[error] No se encontraron *_amr_results.tsv bajo {bdir}/*/05_amrfinder/", file=sys.stderr)
        sys.exit(1)

    rows_out = []
    for amrf_f in amrf_files:
        sample = amrf_f.parent.parent.name
        m = re.match(r'(ds\d+|control)(?:_sim(\d+))?', sample)
        dataset = m.group(1) if m else sample
        sim_idx = f"sim{m.group(2)}" if (m and m.group(2)) else None

        if dataset == 'control' or dataset not in panels:
            continue

        panel = panels[dataset]
        truth_genes = load_sample_truth(args.truthset, dataset, sim_idx)

        amrp_f = bdir / sample / "05_amrprofiler" / f"{sample}_final_results_tool1.csv"
        if not amrp_f.exists():
            amrp_f = None

        amrf_detected = load_amrfinder_detected(amrf_f)
        amrp_detected = load_amrprofiler_detected(amrp_f)

        for tool, detected in [('AMRFinderPlus', amrf_detected), ('AMRProfiler', amrp_detected)]:
            r = compute_specificity_for_sample(panel, truth_genes, detected)
            rows_out.append({
                'sample': sample, 'dataset': dataset, 'sim': sim_idx or 'NA', 'tool': tool,
                **r
            })

    with open(args.out, 'w') as fh:
        fh.write("sample\tdataset\tsim\ttool\tpanel_size\ttruth_size\tneg_class_size\t"
                  "fp_panel\ttn\tSpecificity\tout_of_panel_count\tout_of_panel_genes\n")
        for r in rows_out:
            fh.write(f"{r['sample']}\t{r['dataset']}\t{r['sim']}\t{r['tool']}\t"
                      f"{r['panel_size']}\t{r['truth_size']}\t{r['neg_class_size']}\t"
                      f"{r['fp_panel']}\t{r['tn']}\t{r['specificity']:.4f}\t"
                      f"{r['out_of_panel_count']}\t{r['out_of_panel_genes']}\n")

    print(f"\n[OK] {len(rows_out)} filas escritas en {args.out}", file=sys.stderr)

    from collections import defaultdict as dd
    agg = dd(lambda: {'tn': 0, 'fp': 0, 'n': 0, 'out_of_panel': 0})
    for r in rows_out:
        key = (r['dataset'], r['tool'])
        agg[key]['tn'] += r['tn']
        agg[key]['fp'] += r['fp_panel']
        agg[key]['n'] += 1
        agg[key]['out_of_panel'] += r['out_of_panel_count']

    print("\n=== Resumen agregado de Specificity por dataset x tool ===", file=sys.stderr)
    for (ds, tool), d in sorted(agg.items()):
        spec = d['tn'] / (d['tn'] + d['fp']) if (d['tn'] + d['fp']) > 0 else float('nan')
        print(f"  {ds:<8} {tool:<15} n={d['n']:<4} TN={d['tn']:<4} FP_panel={d['fp']:<4} "
              f"Specificity={spec:.3f}  (detecciones fuera de panel: {d['out_of_panel']})", file=sys.stderr)


if __name__ == "__main__":
    main()
