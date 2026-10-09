#!/usr/bin/env python3
"""DUF2924 域的保守性分布，以及「碱性残基是否落在保守位点」→ fig/conserv_data.json"""
import json, numpy as np, collections
from Bio import AlignIO
from scipy.stats import mannwhitneyu
SP="/tmp/claude-460513/-scratch-users-steorra-analysis/b345d1d6-1039-48d9-85a8-07d814549e6a/scratchpad"
aln=AlignIO.read(f"{SP}/duf2924_dom.aln","fasta")
n=len(aln); L=aln.get_alignment_length()
cols=[]
for j in range(L):
    c=[r.seq[j] for r in aln]
    gap=c.count("-")/n
    aa=[x for x in c if x!="-"]
    if not aa: cols.append(None); continue
    cnt=collections.Counter(aa)
    p=np.array([v/len(aa) for v in cnt.values()])
    H=-(p*np.log2(p)).sum()
    cons=1-H/np.log2(20)                       # 0=全随机 1=完全保守
    basic=sum(cnt.get(x,0) for x in "RK")/len(aa)
    acid =sum(cnt.get(x,0) for x in "DE")/len(aa)
    cols.append({"j":j,"gap":gap,"cons":float(cons),"basic":float(basic),
                 "acid":float(acid),"n_aa":len(aa),
                 "top":cnt.most_common(1)[0][0],
                 "top_frac":cnt.most_common(1)[0][1]/len(aa)})
# 只看覆盖度 >50% 的列（真正属于这个域的位置）
core=[c for c in cols if c and c["gap"]<0.5]
cons=np.array([c["cons"] for c in core])
bas =np.array([c["basic"] for c in core])
# 碱性位点 = 该列 ≥40% 是 R/K
isb=bas>=0.40
u,p=mannwhitneyu(cons[isb],cons[~isb],alternative="greater")
out={"n_seq":n,"L":L,"n_core":len(core),
     "cols":core,
     "basic_cols":int(isb.sum()),
     "cons_basic_median":float(np.median(cons[isb])),
     "cons_other_median":float(np.median(cons[~isb])),
     "mwu_p":float(p),
     "top_conserved":[{"j":c["j"],"aa":c["top"],"frac":c["top_frac"],"cons":c["cons"]}
                      for c in sorted(core,key=lambda d:-d["cons"])[:12]]}
json.dump(out,open("fig/conserv_data.json","w"),indent=1)
print(f"  {n} 条 × {L} 列，覆盖>50% 的核心列 {len(core)}")
print(f"  碱性位点（≥40% R/K）{int(isb.sum())} 个")
print(f"  保守度中位: 碱性位点 {np.median(cons[isb]):.3f}  其余 {np.median(cons[~isb]):.3f}")
print(f"  Mann-Whitney（碱性更保守）p = {p:.2e}")
print("  最保守的 8 列:")
for c in out["top_conserved"][:8]:
    print(f"   列{c['j']:>4}  {c['aa']}  {100*c['frac']:>5.1f}%  保守度 {c['cons']:.3f}")
