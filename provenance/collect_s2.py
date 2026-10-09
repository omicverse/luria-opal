#!/usr/bin/env python3
"""Supp Fig 2 的全部数字 → fig/supp2_data.json。五项证伪检验，含没通过的那一项。"""
import json, collections, re, numpy as np
from scipy.stats import fisher_exact
D="../outputs/fe0970b9-86bd-9a46-2a9e-dde59b863117/discovery/dark"
G="../outputs/fe0970b9-86bd-9a46-2a9e-dde59b863117/discovery"

# ① genome-matched 对照（分母见 dark/matched.py：NC_m=8491 / NB_m=7558）
m=json.load(open(f"{D}/matched_enrichment.json"))
duf=[r for r in m if "DUF2924" in (r.get("dom") or "")]
NC,NB=8491,7558
cm=sum(r["cm"] for r in duf); bm=sum(r["bm"] for r in duf)
A={"clusters":[{"n":r["n"],"cm":r["cm"],"bm":r["bm"],"dom":r["dom"]} for r in duf],
   "case_hit":cm,"case_n":NC,"bg_hit":bm,"bg_n":NB,
   "p":fisher_exact([[cm,NC-cm],[bm,NB-bm]],"two-sided")[1]}

# ② 分类学广度
g2c=dict(l.split()[:2] for l in open(f"{D}/duf2924_genomes.tsv") if len(l.split())>=2)
tax=json.load(open(f"{D}/tax2.json"))["reports"]
mp={}
for r in tax:
    o=(r.get("organism") or {}).get("organism_name","")
    for k in ("accession","current_accession","paired_accession"):
        if r.get(k): mp[re.sub(r"^GC[AF]_","",r[k])]=o
orgs=[mp.get(re.sub(r"^GC[AF]_","",a),"") for a in g2c.values()]
named=[o for o in orgs if o]
PLACE=("uncultured","Candidatus","metagenome")
proper=[o for o in named if not o.startswith(PLACE)]
gen=collections.Counter(o.split()[0] for o in proper)
B={"n_contig":len(g2c),"n_named":len(named),"n_placeholder":len(named)-len(proper),
   "n_genus":len(gen),"top":gen.most_common(10),
   "max_share_pct":100*gen.most_common(1)[0][1]/len(proper) if proper else 0}

# ③ 折叠归属：24 个预测结构的最佳 foldseek 命中
best={}
for f in ["fs20b.tsv","all_fs.tsv","fs20.tsv","duf_fs.tsv"]:
    for l in open(f"{D}/{f}",errors="replace"):
        p=l.rstrip("\n").split("\t")
        if len(p)<8: continue
        try: fid,ev,prob=float(p[2]),float(p[4]),float(p[7])
        except Exception: continue
        if p[0] not in best or ev<best[p[0]][1]: best[p[0]]=(p[1],ev,fid,prob,p[-1][:60])
def grp(q): return "RAMA-like" if q.startswith("R_") else ("non-RAMA" if q.startswith("N_") else "unassigned")
C={"structures":[{"q":q,"group":grp(q),"evalue":v[1],"fident":v[2],"prob":v[3],
                  "target":v[0],"desc":v[4]} for q,v in sorted(best.items(),key=lambda kv:kv[1][1])],
   "confident_threshold":1e-3}
cnt=collections.Counter(grp(q) for q in best)
hit=collections.Counter(grp(q) for q,v in best.items() if v[1]<0.1)
C["counts"]=dict(cnt); C["hits_lt_0.1"]=dict(hit)
r_n,n_n=cnt["RAMA-like"],cnt["non-RAMA"]
C["fisher_p"]=fisher_exact([[hit["RAMA-like"],r_n-hit["RAMA-like"]],
                            [hit["non-RAMA"], n_n-hit["non-RAMA"]]],"two-sided")[1]

# ④ RAMA-like 拆分，以及两组在重组酶邻域上有无差别
rs=json.load(open(f"{D}/rama_split.json"))
rama={x.split("|",1)[-1] for x in rs["rama_like"]}
REC={"Resolvase","Recombinase","Zn_ribbon_recom","Recombinase_2","Phage_integrase"}
nb=collections.defaultdict(dict)                       # contig -> {idx: pfam}
for l in open(f"{D}/duf_nb.domtbl",errors="replace"):
    if l.startswith("#"): continue
    p=l.split()
    if len(p)<8: continue
    mm=re.match(r"^(.*)_(\d+)$",p[0])
    if mm: nb[mm.group(1)].setdefault(int(mm.group(2)),p[3])
pos=collections.defaultdict(list)                      # 183 个 case 位点
for l in open(f"{D}/duf2924_all.domtbl",errors="replace"):
    if l.startswith("#"): continue
    p=l.split()
    if len(p)<5: continue
    sid=p[0].split("|")[-1]
    mm=re.match(r"^(.*)_(\d+)$",sid)
    if mm: pos[mm.group(1)].append((int(mm.group(2)),sid))
tab={"RAMA-like":[0,0],"non-RAMA":[0,0]}               # [有重组酶, 无]
for c,lst in pos.items():
    for idx,sid in lst:
        g="RAMA-like" if sid in rama else "non-RAMA"
        near={nb[c].get(idx+k) for k in (-2,-1,1,2)}
        tab[g][0 if (near & REC) else 1]+=1
pr=fisher_exact([tab["RAMA-like"],tab["non-RAMA"]],"two-sided")[1]
Dd={"n_rama":rs["n_rama"],"n_total":rs["n_total"],
    "recombinase_adjacent":tab,"fisher_p":pr}

# ⑤ GANTC（已在 Fig1 的 meta 里算过，这里复用同一来源）
M=json.load(open("fig/panelD_meta.json"))
E={"effect":M["ganc"]["effect"],"n":M["ganc"]}

out={"A":A,"B":B,"C":C,"D":Dd,"E":E}
json.dump(out,open("fig/supp2_data.json","w"),ensure_ascii=False,indent=1)
print(f"  ① matched: {cm}/{NC} vs {bm}/{NB}  p={A['p']:.2e}")
print(f"  ② 分类: {B['n_contig']} contig / {B['n_named']} 有名 / {B['n_genus']} 属, "
      f"最大属 {B['max_share_pct']:.1f}%")
print(f"  ③ 结构: {len(best)} 个 {dict(cnt)}  E<0.1: {dict(hit)}  p={C['fisher_p']:.3f}")
print(f"  ④ RAMA-like {Dd['n_rama']}/{Dd['n_total']}；重组酶邻域 "
      f"RAMA-like {Dd['recombinase_adjacent']['RAMA-like']} vs "
      f"non {Dd['recombinase_adjacent']['non-RAMA']}  p={Dd['fisher_p']:.3f}")
print(f"  ⑤ GANTC {E['effect']['GANTC']['d']:+.3f} SD, GATC {E['effect']['GATC']['d']:+.3f} SD")
