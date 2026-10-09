#!/usr/bin/env python3
"""183 条 DUF2924 域的系统发育树 + 每个叶的注释 → fig/tree_data.json

树: duf2924_clean.faa → 按 domtbl 的 ali 坐标切域 → mafft --auto → FastTree -lg → 中点定根
注释: RAMA-like（rama_split.json）、是否融合体（tlen>400）、±2 内有无重组酶 / 甲基化酶
"""
import json, re, collections
from Bio import Phylo
D="../outputs/fe0970b9-86bd-9a46-2a9e-dde59b863117/discovery/dark"
SP="/tmp/claude-460513/-scratch-users-steorra-analysis/b345d1d6-1039-48d9-85a8-07d814549e6a/scratchpad"

t=Phylo.read(f"{SP}/duf2924_dom.nwk","newick")
t.root_at_midpoint()
rama={x.split("|",1)[-1] for x in json.load(open(f"{D}/rama_split.json"))["rama_like"]}
tlen={}
pos={}
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
    if m: nb[m.group(1)].setdefault(int(m.group(2)),p[3])
REC={"Resolvase","Recombinase","Zn_ribbon_recom","Recombinase_2","Phage_integrase"}
MET={"N6_N4_Mtase","N6_Mtase","DNA_methylase","Methylase_S","MmeI_Mtase","Dam",
     "HTH_ParB_Mtase","ResIII","N6-adenine"}

def ann(sid):
    c,i=pos.get(sid,(None,None))
    near={nb.get(c,{}).get(i+k) for k in (-2,-1,1,2)} if c else set()
    return {"rama":sid in rama or f"case|{sid}" in rama,
            "fusion":tlen.get(sid,0)>400,
            "rec":bool(near&REC),"met":bool(near&MET),"len":tlen.get(sid,0)}

# 叶的绘图坐标：x 为到根的枝长，y 为叶序
leaves=t.get_terminals()
depth=t.depths()
tips=[]
for k,c in enumerate(leaves):
    sid=c.name
    tips.append({"id":sid,"x":depth[c],"y":k,**ann(sid)})
# 内部节点线段
segs=[]
def walk(cl):
    ys=[]
    for ch in cl.clades:
        ys.append(walk(ch))
    if not cl.clades:
        return leaves.index(cl)
    y0,y1=min(ys),max(ys); ym=sum(ys)/len(ys)
    segs.append({"type":"v","x":depth[cl],"y0":y0,"y1":y1})
    for ch,yy in zip(cl.clades,ys):
        segs.append({"type":"h","x0":depth[cl],"x1":depth[ch],"y":yy})
    return ym
walk(t.root)
out={"tips":tips,"segs":segs,"n":len(tips),
     "counts":{"rama":sum(1 for d in tips if d["rama"]),
               "fusion":sum(1 for d in tips if d["fusion"]),
               "rec":sum(1 for d in tips if d["rec"]),
               "met":sum(1 for d in tips if d["met"])}}
json.dump(out,open("fig/tree_data.json","w"),indent=1)
print(f"  {out['n']} 叶, {len(segs)} 条线段")
print(f"  注释: {out['counts']}")
