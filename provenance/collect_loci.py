#!/usr/bin/env python3
"""挑 4 个不同物种的 DUF2924 基因座 → fig/supp_loci.json。

来源：dark/duf_nb.domtbl（1,572 条邻域蛋白的域注释，hmmsearch 口径：
col1=序列 col3=长度 col4=Pfam名 col7=E-value）、dark/duf2924_genomes.tsv
（contig→基因组）、dark/tax2.json（基因组→物种名）。
"""
import json, collections, re
D="../outputs/fe0970b9-86bd-9a46-2a9e-dde59b863117/discovery/dark"

# 每条蛋白取 E-value 最小的域
best={}
for l in open(f"{D}/duf_nb.domtbl",errors="replace"):
    if l.startswith("#"): continue
    p=l.split()
    if len(p)<8: continue
    sid,tlen,pf,ev=p[0],int(p[2]),p[3],float(p[6])
    if sid not in best or ev<best[sid][2]: best[sid]=(pf,tlen,ev)

loci=collections.defaultdict(dict)        # contig -> {gene_idx: (pfam,len)}
for sid,(pf,tlen,ev) in best.items():
    m=re.match(r"^(.*)_(\d+)$",sid)
    if not m: continue
    loci[m.group(1)][int(m.group(2))]=(pf,tlen)

# DUF2924 自己的位置来自 duf2924_all.domtbl（183 条命中，ID 带 case| 前缀）
duf=collections.defaultdict(list)
for l in open(f"{D}/duf2924_all.domtbl",errors="replace"):
    if l.startswith("#"): continue
    p=l.split()
    if len(p)<5: continue
    sid=p[0].split("|")[-1]
    m=re.match(r"^(.*)_(\d+)$",sid)
    if m: duf[m.group(1)].append((int(m.group(2)),p[3],int(p[2])))

# contig → 基因组 → 物种
g2c={}
for l in open(f"{D}/duf2924_genomes.tsv"):
    p=l.split()
    if len(p)>=2: g2c[p[0]]=p[1]
tax=json.load(open(f"{D}/tax2.json"))["reports"]
acc2org={}
for r in tax:
    o=r.get("organism") or {}
    for k in ("accession","current_accession","paired_accession"):
        if r.get(k): acc2org[re.sub(r"^GC[AF]_","",r[k])]=o.get("organism_name","")
def org(contig):
    a=g2c.get(contig)
    return acc2org.get(re.sub(r"^GC[AF]_","",a),"") if a else ""

REC={"Resolvase","Recombinase","Zn_ribbon_recom","Phage_integrase","Recombinase_2"}
MET={"N6_N4_Mtase","N6_Mtase","DNA_methylase","Methylase_S","MmeI_Mtase","Dam",
     "HTH_ParB_Mtase","N6-adenine"}
cand=[]
for c,pos in duf.items():
    genes=loci.get(c,{})
    if len(genes)<3: continue
    idx=pos[0][0]
    near={g:v for g,v in genes.items() if abs(g-idx)<=4}
    near[idx]=("DUF2924",pos[0][2])
    pfs={v[0] for v in near.values()}
    cand.append({"contig":c,"organism":org(c),"duf_idx":idx,"duf_pfam":pos[0][1],
                 "genes":{str(g):{"pfam":v[0],"len":v[1]} for g,v in sorted(near.items())},
                 "has_rec":bool(pfs&REC),"has_met":bool(pfs&MET),"n":len(near)})
# 优先：正式双名种 > 有名字 > 无名；再按「同时带重组酶和甲基化酶」「基因数多」排
VAGUE=("uncultured","bacterium","Candidatus","sp.","metagenome")
def named(o):
    if not o: return 0
    w=o.split()
    if len(w)>=2 and w[0][:1].isupper() and not any(v in o for v in VAGUE): return 2
    return 1
cand.sort(key=lambda d:(-named(d["organism"]),-(d["has_rec"]+d["has_met"]),-d["n"]))
pick=[]; seen=set()
for d in cand:
    g=(d["organism"].split()[0] if d["organism"] else d["contig"][:6])
    if g in seen: continue
    seen.add(g); pick.append(d)
    if len(pick)==4: break
json.dump(pick,open("fig/supp_loci.json","w"),ensure_ascii=False,indent=1)
print(f"  候选 {len(cand)} 个位点，选出 {len(pick)}：")
for d in pick:
    print(f"   {d['organism'] or d['contig']:<46} genes={d['n']} "
          f"rec={d['has_rec']} met={d['has_met']}  {sorted(v['pfam'] for v in d['genes'].values())}")
