# Independent recomputation

Each script reads only from `tables/` or `ledger/`, recomputes a statistic the manuscript reports,
and prints both values side by side. None of them reads a figure, a cached result, or anything the
agents wrote as prose.

| script | recomputes | runtime |
|---|---|---|
| `verify_ledger.py` | event counts per class, campaigns, span | instant |
| `verify_screen.py` | Fisher p and BH q for all 4,857 families; how many reach q < 0.05 | ~3 min |
| `verify_reverse.py` | serine and tyrosine denominators, the 58-fold odds ratio, the distance distribution | ~1 min, 2 GB RAM |

Requires `scipy` and `statsmodels`.

Expected output of `verify_reverse.py`:

```
serine   loci   524,414   with DUF2924  14,088    2.686%
tyrosine loci 1,403,550   with DUF2924     665   0.0474%
odds ratio    58.237
```

Expected output of `verify_screen.py`:

```
families in table           : 4857
families at q < 0.05        : 32        (reported: 32)
DUF2924 / PF11149.14
  Fisher p  recomputed      : 2.271e-10   (table: 2.27e-10)
  BH q      recomputed      : 9.193e-08   (table: 9.19e-08)
```
