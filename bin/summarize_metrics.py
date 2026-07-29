#!/usr/bin/env python3
"""
Genera la tabla resumen agregada (promedios de TP/FP/FN/Sensitivity/
Precision/F1 por dataset x herramienta) a partir de los *_metrics.tsv
finales del pipeline. Pensada para pegar directo en la memoria del TFM.

Uso:
    python3 summarize_metrics.py <benchmark_results_dir> [--out summary.tsv]
"""
import sys
import csv
import argparse
from pathlib import Path
from collections import defaultdict


def calc_metrics(tp, fp, fn):
    sens = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    prec = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    f1 = 2 * sens * prec / (sens + prec) if (sens + prec) > 0 else 0.0
    return sens, prec, f1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('benchmark_dir')
    ap.add_argument('--out', default=None)
    args = ap.parse_args()

    base = Path(args.benchmark_dir)
    files = sorted(base.glob("*/06_benchmark/*_metrics.tsv"))

    if not files:
        print(f"No se encontraron *_metrics.tsv bajo {base}/*/06_benchmark/", file=sys.stderr)
        sys.exit(1)

    # Agregado "micro" (suma de TP/FP/FN across muestras, luego calcula métricas)
    agg_sum = defaultdict(lambda: {'tp': 0, 'fp': 0, 'fn': 0, 'n': 0})
    # Agregado "macro" (promedio de Sensitivity/Precision/F1 por-muestra)
    agg_avg = defaultdict(lambda: {'sens': [], 'prec': [], 'f1': [], 'n': 0})

    for f in files:
        with open(f) as fh:
            for row in csv.DictReader(fh, delimiter='\t'):
                dataset = row['dataset']
                tool = row['tool']
                key = (dataset, tool)

                tp, fp, fn = int(row['TP']), int(row['FP']), int(row['FN'])
                agg_sum[key]['tp'] += tp
                agg_sum[key]['fp'] += fp
                agg_sum[key]['fn'] += fn
                agg_sum[key]['n'] += 1

                agg_avg[key]['sens'].append(float(row['Sensitivity']))
                agg_avg[key]['prec'].append(float(row['Precision']))
                agg_avg[key]['f1'].append(float(row['F1']))
                agg_avg[key]['n'] += 1

    out_lines = []
    header = (f"{'Dataset':<10} {'Tool':<15} {'n':>3} | "
              f"{'TP':>4} {'FP':>4} {'FN':>4} | "
              f"{'Sens(micro)':>11} {'Prec(micro)':>11} {'F1(micro)':>9} | "
              f"{'Sens(macro)':>11} {'Prec(macro)':>11} {'F1(macro)':>9}")
    out_lines.append(header)
    out_lines.append("-" * len(header))

    for key in sorted(agg_sum.keys()):
        dataset, tool = key
        d = agg_sum[key]
        sens_micro, prec_micro, f1_micro = calc_metrics(d['tp'], d['fp'], d['fn'])

        a = agg_avg[key]
        n = a['n']
        sens_macro = sum(a['sens']) / n if n else 0.0
        prec_macro = sum(a['prec']) / n if n else 0.0
        f1_macro = sum(a['f1']) / n if n else 0.0

        out_lines.append(
            f"{dataset:<10} {tool:<15} {d['n']:>3} | "
            f"{d['tp']:>4} {d['fp']:>4} {d['fn']:>4} | "
            f"{sens_micro:>11.3f} {prec_micro:>11.3f} {f1_micro:>9.3f} | "
            f"{sens_macro:>11.3f} {prec_macro:>11.3f} {f1_macro:>9.3f}"
        )

    output = "\n".join(out_lines)
    print(output)

    if args.out:
        with open(args.out, 'w') as fh:
            fh.write(output + "\n")
        print(f"\n[OK] guardado en {args.out}", file=sys.stderr)

    # Tabla markdown lista para pegar en la memoria (macro, más estándar para reportar)
    print("\n\n=== Tabla Markdown (macro-promedio, lista para la memoria) ===\n")
    md_lines = [
        "| Dataset | Herramienta | Sensitivity | Precision | F1 |",
        "|---|---|---|---|---|",
    ]
    for key in sorted(agg_avg.keys()):
        dataset, tool = key
        a = agg_avg[key]
        n = a['n']
        sens_macro = sum(a['sens']) / n if n else 0.0
        prec_macro = sum(a['prec']) / n if n else 0.0
        f1_macro = sum(a['f1']) / n if n else 0.0
        md_lines.append(f"| {dataset} | {tool} | {sens_macro:.3f} | {prec_macro:.3f} | {f1_macro:.3f} |")
    print("\n".join(md_lines))


if __name__ == "__main__":
    main()
