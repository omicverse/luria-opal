#!/usr/bin/env python3
"""新 Panel D 的数据 → fig/panelD2_data.json
① 24 个预测结构的表面电荷  ② 183 条 vs 对照的净电荷分布  ③ foldseek 最佳命中"""
import json, glob, os, re, collections
import numpy as np, statistics as st
from Bio import SeqIO
import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx
from scipy.stats import mannwhitneyu
import warnings; warnings.filterwarnings("ignore")

D="../outputs/fe0970b9-86bd-9a46-2a9e-dde59b863117/discovery/dark"
SP="/tmp/claude-460513/-scratch-users-steorra-analysis/b345d1d6-1039-48d9-85a8-07d814549e6a/scratchpad"
POS={"ARG","LYS"}; NEG={"ASP","GLU"}

# ① 结构表面电荷
surf=[]
for p in sorted(glob.glob(f"{D}/fsq2/*.cif"))+sorted(glob.glob(f"{D}/fsq/*.cif")):
    try: a=pdbx.get_structure(pdbx.CIFFile.read(p),model=1)
    except Exception: continue
    a=a[struc.filter_amino_acids(a)]
    if len(a)<100: continue
    sa=struc.apply_residue_wise(a,struc.sasa(a,vdw_radii="Single"),np.nansum)
    ids=np.unique(a.res_id); names=[a.res_name[a.res_id==r][0] for r in ids]
    ex=sa>30
    surf.append({"q":os.path.basename(p)[:-4],
                 "n":len(ids),"nsurf":int(ex.sum()),
                 "pos":sum(1 for n_,e in zip(names,ex) if e and n_ in POS),
                 "neg":sum(1 for n_,e in zip(names,ex) if e and n_ in NEG)})
for s in surf: s["net"]=s["pos"]-s["neg"]

# ② 序列净电荷
def q(s): s=s.upper(); return (s.count("R")+s.count("K"))-(s.count("D")+s.count("E"))
def p100(s): return 100*q(s)/max(len(s),1)
duf=[str(r.seq).rstrip("*") for r in SeqIO.parse(f"{D}/duf2924_clean.faa","fasta")]
nb=[str(r.seq).rstrip("*") for r in SeqIO.parse(f"{D}/duf_neighbors.faa","fasta")]
rec=[str(r.seq).rstrip("*") for r in SeqIO.parse(f"{SP}/recomb.faa","fasta")]
nbm=[s for s in nb if 120<=len(s)<=200]
groups={"DUF2924":duf,"neighbours\n(120-200 aa)":nbm,"serine\nrecombinases":rec}
charge={k:{"n":len(v),"per100":[p100(s) for s in v],
           "median":st.median([p100(s) for s in v]),
           "frac_pos":sum(1 for s in v if q(s)>0)/len(v)} for k,v in groups.items()}
mwu=float(mannwhitneyu(charge["DUF2924"]["per100"],
                       charge["neighbours\n(120-200 aa)"]["per100"],alternative="greater")[1])

# ③ foldseek 最佳命中
best={}
for f in ["fs20b.tsv","all_fs.tsv","fs20.tsv","duf_fs.tsv"]:
    for l in open(f"{D}/{f}",errors="replace"):
        c=l.rstrip("\n").split("\t")
        if len(c)<8: continue
        try: ev=float(c[4])
        except ValueError: continue
        if c[0] not in best or ev<best[c[0]][0]: best[c[0]]=(ev,c[1],c[-1][:70])
folds=[{"q":k,"evalue":v[0],"target":v[1],"desc":v[2],
        "mpnd":bool(re.search(r"MPND",v[2],re.I)),
        "group":"RAMA-like" if k.startswith("R_") else ("non-RAMA" if k.startswith("N_") else "other")}
       for k,v in sorted(best.items(),key=lambda kv:kv[1][0])]
out={"surface":surf,"charge":charge,"mwu_p":mwu,"folds":folds,
     "surface_median_net":st.median([s["net"] for s in surf])}
json.dump(out,open("fig/panelD2_data.json","w"),indent=1)
print(f"  ① 结构 {len(surf)} 个，表面净电荷中位 {out['surface_median_net']:+.0f}，"
      f"全为正 {all(s['net']>0 for s in surf)}")
CK="neighbours\n(120-200 aa)"
print(f"  ② DUF2924 {charge['DUF2924']['median']:+.2f}/100aa  vs 对照 "
      f"{charge[CK]['median']:+.2f}  p={mwu:.1e}")
print(f"  ③ foldseek {len(folds)} 个，MPND 命中 {sum(1 for f in folds if f['mpnd'])} 个")
