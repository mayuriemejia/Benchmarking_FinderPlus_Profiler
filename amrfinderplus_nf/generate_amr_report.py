#!/usr/bin/env python3
import sys
from pathlib import Path

from openpyxl import Workbook


def main() -> int:
    if len(sys.argv) != 4:
        print(
            "Usage: generate_amr_report.py <amrfinder_results.tsv> <out.xlsx> <sample_id>",
            file=sys.stderr,
        )
        return 2

    tsv_path = Path(sys.argv[1])
    out_xlsx = Path(sys.argv[2])
    sample_id = sys.argv[3]

    if not tsv_path.exists():
        print(f"ERROR: TSV not found: {tsv_path}", file=sys.stderr)
        return 1

    wb = Workbook()
    ws = wb.active
    ws.title = "AMR"
    ws.append(["sample_id", sample_id])
    ws.append([])

    with tsv_path.open("r", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n")
        if header:
            ws.append(header.split("\t"))
        for line in fh:
            ws.append(line.rstrip("\n").split("\t"))

    out_xlsx.parent.mkdir(parents=True, exist_ok=True)
    wb.save(out_xlsx)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

