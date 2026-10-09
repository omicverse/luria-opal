#!/usr/bin/env python3
"""Recompute the Pfam enrichment screen from tables/family_enrichment.tsv.

Checks three claims made about the screen:
  1. how many families were tested, which is the Benjamini-Hochberg denominator
  2. the Fisher p value and the BH q value for DUF2924 (PF11149)
  3. how many families reach q < 0.05

The screen is one-sided (alternative: greater). It asks whether a family is
over-represented beside integrases, not whether it differs in either direction.
Run two-sided instead and 46 families pass rather than 32; the switch is the
 argument below.

Run:  python verify/verify_screen.py
"""
import csv, sys, pathlib
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

T = pathlib.Path(__file__).resolve().parent.parent / "tables" / "family_enrichment.tsv"
rows = list(csv.DictReader(open(T), delimiter="\t"))

N_CASE = N_BG = 10000           # proteins sampled from each side, seed 42
print(f"families in table           : {len(rows)}")

p = []
for r in rows:
    c, b = int(r["case"]), int(r["bg"])
    p.append(fisher_exact([[c, N_CASE - c], [b, N_BG - b]], alternative="greater")[1])
q = multipletests(p, method="fdr_bh")[1]

sig = sum(1 for x in q if x < 0.05)
print(f"families at q < 0.05        : {sig}        (reported: 32)")

for i, r in enumerate(rows):
    if r["acc"].startswith("PF11149"):
        print(f"\nDUF2924 / {r['acc']}")
        print(f"  case / background         : {r['case']} / {r['bg']}")
        print(f"  Fisher p  recomputed      : {p[i]:.3e}   (table: {r['p_all']})")
        print(f"  BH q      recomputed      : {q[i]:.3e}   (table: {r['q_fdr']})")
        print(f"  BH denominator used here  : {len(rows)}")
        break

print("\nNote. The campaign deliverable (deliverables/screen.zh.md) states the BH correction")
print("was taken over 30,134 tests, the number of Pfam models. The q values in the table")
print("match a correction over the families actually testable, which is the row count above.")
