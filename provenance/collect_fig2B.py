#!/usr/bin/env python3
"""从 Panel A 的四个类别里各取代表基因座 → fig/fig2B_loci.json（先出代表名单与邻域序列）"""
import json, re, collections, glob, os
from Bio import SeqIO
SP="/tmp/claude-460513/-scratch-users-steorra-analysis/b345d1d6-1039-48d9-85a8-07d814549e6a/scratchpad"
D="../outputs/fe0970b9-86bd-9a46-2a9e-dde59b863117/discovery/dark"
R="../outputs/fe0970b9-86bd-9a46-2a9e-dde59b863117/discovery/reverse"
RT=json.load(open("fig/rtree_data.json"))
meta=json.load(open(f"{SP}/recomb_meta.json"))
dom=collections.defaultdict(set)
for l in open(f"{SP}/recomb_lsr.domtbl",errors="replace"):
    if l.startswith("#"): continue
    p=l.split()
    if len(p)>=8: dom[p[0]].add(p[3])
from Bio import Phylo
t=Phylo.read(f"{SP}/recomb_dom.nwk","newick"); t.root_at_midpoint()
leaves=[c.name for c in t.get_terminals()]
def grp(n):
    lsr=bool(dom.get(n)); duf=meta[n]["grp"]=="with DUF2924"
    return ("LSR + DUF2924" if (lsr and duf) else "LSR, no DUF2924" if lsr else
            "small serine + DUF2924" if duf else "small serine, no DUF2924")
# 每类取最大楔形里的 2 条（small serine + DUF2924 只有 3 条单支，取 1）
byg=collections.defaultdict(list)
for n in leaves: byg[grp(n)].append(n)
PICK={"LSR + DUF2924":2,"LSR, no DUF2924":2,"small serine + DUF2924":1,"small serine, no DUF2924":1}
sel=[]
for g,k in PICK.items():
    for n in byg[g][:k]: sel.append({"rec":n,"group":g})
# 取每个代表的 ±4 邻域（带坐标与链向）
HDR=re.compile(r"^>(\S+)\s+#\s+(\d+)\s+#\s+(\d+)\s+#\s+(-?1)\s+#")
want=collections.defaultdict(set)
for s_ in sel:
    c,i=s_["rec"].rsplit("_",1); i=int(i)
    for o in range(-4,5): want[c].add(i+o)
coords=collections.defaultdict(dict); seqs={}
for p in [f"{D}/dufctx/chunk_{{}}.faa"]+sorted(glob.glob(f"{R}/ctrlctx/chunk_*.faa")):
    if not os.path.exists(p): continue
    for r in SeqIO.parse(p,"fasta"):
        m=HDR.match(">"+r.description)
        if not m: continue
        sid=m.group(1); cc,ii=sid.rsplit("_",1); ii=int(ii)
        if cc in want and ii in want[cc] and ii not in coords[cc]:
            coords[cc][ii]=(int(m.group(2)),int(m.group(3)),int(m.group(4)))
            seqs[sid]=str(r.seq)
SeqIO.write([SeqIO.SeqRecord(__import__("Bio.Seq",fromlist=["Seq"]).Seq(v),id=k,description="")
             for k,v in seqs.items()],f"{SP}/fig2B_nb.faa","fasta")
json.dump({"sel":sel,"coords":{c:{str(k):v for k,v in d.items()} for c,d in coords.items()}},
          open(f"{SP}/fig2B_raw.json","w"))
print(f"  代表 {len(sel)} 条:")
for s_ in sel: print(f"   {s_['group']:<26} {s_['rec']}")
print(f"  邻域蛋白 {len(seqs)} 条 → {SP}/fig2B_nb.faa")
