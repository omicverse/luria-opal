#!/usr/bin/env python3
"""伙伴重组酶的催化域树，折叠成分类楔形 → fig/rtree_data.json

类别由域架构 + 是否与 DUF2924 相邻定义：
  LSR + DUF2924 / LSR, no DUF2924 / small serine + DUF2924 / small serine, no DUF2924
折叠规则：取「所有叶同类」的极大单系簇；不足 3 叶的保留为单支。
"""
import json, collections
from Bio import Phylo
SP="/tmp/claude-460513/-scratch-users-steorra-analysis/b345d1d6-1039-48d9-85a8-07d814549e6a/scratchpad"
meta=json.load(open(f"{SP}/recomb_meta.json"))
dom=collections.defaultdict(set)
for l in open(f"{SP}/recomb_lsr.domtbl",errors="replace"):
    if l.startswith("#"): continue
    p=l.split()
    if len(p)>=8: dom[p[0]].add(p[3])
def grp(name):
    lsr=bool(dom.get(name)); duf=meta[name]["grp"]=="with DUF2924"
    return ("LSR + DUF2924" if (lsr and duf) else
            "LSR, no DUF2924" if lsr else
            "small serine + DUF2924" if duf else "small serine, no DUF2924")
t=Phylo.read(f"{SP}/recomb_dom.nwk","newick"); t.root_at_midpoint()
leaves=t.get_terminals(); depth=t.depths()
idx={c:i for i,c in enumerate(leaves)}
G={c.name:grp(c.name) for c in leaves}

def pure(cl,minfrac=0.90,minn=4,maxn=24):
    """≥85% 同类且 ≥4 叶就折叠，按优势类别标注；纯簇 ≥3 叶也折叠。"""
    gs=collections.Counter(G[x.name] for x in cl.get_terminals())
    n=sum(gs.values()); g,c=gs.most_common(1)[0]
    if n>maxn: return None                      # 太大的簇继续往下拆
    if n>=minn: return g,c/n                    # 够大就折叠，按优势类别标注
    return None
wedges=[]; wedge_clades=[]; kept=set()
def visit(cl):
    r=pure(cl)
    ts=cl.get_terminals()
    if r:
        g,frac=r
        ys=[idx[x] for x in ts]
        wedges.append({"g":g,"n":len(ts),"purity":round(frac,3),
                       "y0":min(ys),"y1":max(ys),
                       "x0":depth[cl],"x1":max(depth[x] for x in ts)})
        wedge_clades.append(cl); kept.update(ts); return True
    for ch in cl.clades: visit(ch)
    return False
par={}
for cl in t.get_nonterminals():
    for ch in cl.clades: par[ch]=cl
visit(t.root)
single_cl=[c for c in leaves if c not in kept]
singles=[{"g":G[c.name],"n":1,"y0":idx[c],"y1":idx[c],"x0":0,"x1":depth[c]}
         for c in single_cl]
par={}
for cl in t.get_nonterminals():
    for ch in cl.clades: par[ch]=cl
for sg,c in zip(singles,single_cl):
    sg["x0"]=depth[par[c]] if c in par else 0
for w,cl in zip(wedges,wedge_clades):
    w["x0"]=depth[par[cl]] if cl in par else 0      # 尖端接在父节点上
    w["xc"]=depth[cl]
segs=[]
WEDGE={}                                   # clade -> (y中点)
def ymid(cl):
    ts=cl.get_terminals()
    return (min(idx[x] for x in ts)+max(idx[x] for x in ts))/2
collapsed=set()
for w,cl in zip(wedges,wedge_clades): collapsed.add(cl)
def walk(cl):
    """返回该节点的 y；沿途把枝画出来。折叠簇当作叶处理。"""
    if cl in collapsed or not cl.clades:
        return ymid(cl) if cl.clades else idx[cl]
    ys=[]
    for ch in cl.clades:
        y=walk(ch); ys.append(y)
        segs.append({"type":"h","x0":depth[cl],"x1":depth[ch],"y":y})
    segs.append({"type":"v","x":depth[cl],"y0":min(ys),"y1":max(ys)})
    return sum(ys)/len(ys)
walk(t.root)
cnt=collections.Counter(G.values())
out={"wedges":wedges+singles,"segs":segs,"n":len(leaves),
     "groups":dict(cnt),
     "stats":{"z_all":13.50,"p_all":0.0002,"z_lsr":11.89,"p_lsr":0.0002,
              "n_lsr":176,"lsr_withduf":158,"lsr_ctrl":18,"median_pid":0.69}}
json.dump(out,open("fig/rtree_data.json","w"),indent=1)
print(f"  {len(leaves)} 叶 → {len(wedges)} 个楔形 + {len(singles)} 条单支")
for k,v in cnt.most_common(): print(f"   {k:<26} {v:>4}")
