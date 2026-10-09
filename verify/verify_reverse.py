#!/usr/bin/env python3
"""Recompute the reverse rate from the position tables, independently of
tables/reverse_rate.json.

Every serine recombinase locus and every tyrosine integrase locus in the
proteome search is taken in turn, and the presence of a DUF2924 protein within
five genes is recorded. The denominators are therefore self-contained: no
case/background sampling enters this number.

Run:  python verify/verify_reverse.py       (needs ~2 GB RAM, ~1 min)
"""
import json, pathlib, collections

T = pathlib.Path(__file__).resolve().parent.parent / "tables"
duf = json.load(open(T / "duf2924_pos.json"))
ser = json.load(open(T / "serine_pos.json"))
tyr = json.load(open(T / "tyrosine_pos.json"))
ref = json.load(open(T / "reverse_rate.json"))

def rate(target, window=5):
    """loci of `target` that carry a DUF2924 protein within `window` genes"""
    n = hit = 0
    hist = collections.Counter()
    for contig, idxs in target.items():
        d = set(duf.get(contig, []))
        for i in idxs:
            n += 1
            if not d:
                continue
            near = [abs(o) for o in range(-window, window + 1) if i + o in d]
            if near:
                hit += 1
                hist[min(near)] += 1
    return n, hit, hist

n_s, h_s, hist_s = rate(ser)
n_t, h_t, hist_t = rate(tyr)
odds = (h_s / (n_s - h_s)) / (h_t / (n_t - h_t))

print(f"serine   loci {n_s:>9,}   with DUF2924 {h_s:>7,}   {100*h_s/n_s:6.3f}%")
print(f"tyrosine loci {n_t:>9,}   with DUF2924 {h_t:>7,}   {100*h_t/n_t:6.4f}%")
print(f"odds ratio    {odds:.3f}")
print(f"\nreported in tables/reverse_rate.json:")
print(f"  n_ser {ref['n_ser']:,}  ser_with_duf {ref['ser_with_duf']:,}")
print(f"  n_tyr {ref['n_tyr']:,}  tyr_with_duf {ref['tyr_with_duf']:,}")
print(f"  odds_ratio {ref['odds_ratio']}")

print("\ndistance distribution, serine side (% of DUF2924-positive serine loci)")
for k in range(0, 6):
    v = hist_s.get(k, 0)
    print(f"  {k} gene{'s' if k!=1 else ' '}: {v:>6,}  {100*v/max(h_s,1):5.1f}%")
