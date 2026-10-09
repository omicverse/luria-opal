#!/usr/bin/env python3
"""Figure 2 的数据 → fig/fig2_data.json"""
import json, re, collections, glob
import numpy as np
D="../outputs/fe0970b9-86bd-9a46-2a9e-dde59b863117/discovery/dark"
R="../outputs/fe0970b9-86bd-9a46-2a9e-dde59b863117/discovery/reverse"

# ── 每条 DUF2924 的位置与 ±5 邻域域 ──
pos={}
tlen={}
for l in open(f"{D}/duf2924_all.domtbl",errors="replace"):
    if l.startswith("#"): continue
    p=l.split()
    if len(p)<5: continue
    sid=p[0].split("|")[-1]; tlen[sid]=int(p[2])
    m=re.match(r"^(.*)_(\d+)$",sid)
    if m: pos[sid]=(m.group(1),int(m.group(2)))
nb=collections.defaultdict(dict)
for l in open(f"{D}/duf_nb.domtbl",errors="replace"):
    if l.startswith("#"): continue
    p=l.split()
    if len(p)<8: continue
    m=re.match(r"^(.*)_(\d+)$",p[0])
    if not m: continue
    c,i=m.group(1),int(m.group(2)); ev=float(p[6])
    if i not in nb[c] or ev<nb[c][i][1]: nb[c][i]=(p[3],ev)

# ── A：±5 每个位置上出现什么域 ──
OFF=list(range(-5,6))
cnt={o:collections.Counter() for o in OFF}
nloc=0
for sid,(c,i) in pos.items():
    if c not in nb: continue
    nloc+=1
    for o in OFF:
        if o==0: cnt[o]["DUF2924"]+=1; continue
        d=nb[c].get(i+o)
        if d: cnt[o][d[0]]+=1
tot=collections.Counter()
for o in OFF:
    if o: tot.update(cnt[o])
TOP=[d for d,_ in tot.most_common(14)]
M=[[cnt[o][d] for o in OFF] for d in TOP]
neigh={"offsets":OFF,"domains":TOP,"matrix":M,"n_loci":nloc,
       "annotated":{str(o):sum(cnt[o].values()) for o in OFF}}

# ── B：长度分布 ──
lens=sorted(tlen.values())

# ── C：正反向距离分布 ──
rr=json.load(open(f"{R}/reverse_rate.json"))

# ── D：foldseek 最佳命中 ──
best={}
for f in ["fs20b.tsv","all_fs.tsv","fs20.tsv","duf_fs.tsv"]:
    for l in open(f"{D}/{f}",errors="replace"):
        p=l.rstrip("\n").split("\t")
        if len(p)<8: continue
        try: ev=float(p[4]); fid=float(p[2]); prob=float(p[7])
        except Exception: continue
        if p[0] not in best or ev<best[p[0]][0]: best[p[0]]=(ev,fid,prob,p[1],p[-1][:60])
folds=[{"q":q,"evalue":v[0],"fident":v[1],"prob":v[2],"target":v[3],"desc":v[4],
        "group":"RAMA-like" if q.startswith("R_") else ("non-RAMA" if q.startswith("N_") else "other")}
       for q,v in sorted(best.items(),key=lambda kv:kv[1][0])]

# ── E：甲基化酶出现在哪个位置 ──
MET={"N6_N4_Mtase","N6_Mtase","DNA_methylase","Methylase_S","MmeI_Mtase","Dam",
     "HTH_ParB_Mtase","ResIII","N6-adenine"}
mp=collections.Counter()
for sid,(c,i) in pos.items():
    if c not in nb: continue
    for o in range(-5,6):
        d=nb[c].get(i+o)
        if o and d and d[0] in MET: mp[o]+=1
out={"neigh":neigh,"lens":lens,"reverse":rr,"folds":folds,
     "methyl_pos":{str(k):v for k,v in sorted(mp.items())}}
json.dump(out,open("fig/fig2_data.json","w"),indent=1)
print(f"  邻域: {nloc} 个位点, top14 域, 各位置注释数 "
      f"{[neigh['annotated'][str(o)] for o in OFF]}")
print(f"  前 6 域: {TOP[:6]}")
print(f"  长度: 中位 {lens[len(lens)//2]}  融合体 {[x for x in lens if x>400]}")
print(f"  folds: {len(folds)} 个结构")
print(f"  甲基化位置: {dict(sorted(mp.items()))}")
