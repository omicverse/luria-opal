#!/usr/bin/env python3
"""带坐标与链方向的基因座 → fig/supp2_loci.json。

DUF2924 侧：dark/dufctx/chunk_{}.faa（prodigal 头，含 start/end/strand）
           + dark/duf_nb.domtbl（邻域域注释）+ dark/duf2924_all.domtbl（DUF2924 位置）
对照侧：  reverse/ctrlctx/*.faa + reverse/ctrl_recombinase_pos.tsv（重组酶位置与类型）
           + reverse/ctrl_methyl.domtbl（只扫了甲基化 HMM，故其余基因标为未注释）
"""
import json, re, glob, collections
D="../outputs/fe0970b9-86bd-9a46-2a9e-dde59b863117/discovery/dark"
R="../outputs/fe0970b9-86bd-9a46-2a9e-dde59b863117/discovery/reverse"
HDR=re.compile(r"^>(\S+)\s+#\s+(\d+)\s+#\s+(\d+)\s+#\s+(-?1)\s+#")

def coords(paths, want):
    """只读 header，取 want 里这些 contig 的基因坐标。"""
    out=collections.defaultdict(dict)
    for p in paths:
        for l in open(p,errors="replace"):
            if l[0]!=">": continue
            m=HDR.match(l)
            if not m: continue
            sid=m.group(1); mm=re.match(r"^(.*)_(\d+)$",sid)
            if not mm or mm.group(1) not in want: continue
            out[mm.group(1)][int(mm.group(2))]=(int(m.group(2)),int(m.group(3)),int(m.group(4)))
    return out

# ── DUF2924 侧 ──
pfam=collections.defaultdict(dict)
for l in open(f"{D}/duf_nb.domtbl",errors="replace"):
    if l.startswith("#"): continue
    p=l.split()
    if len(p)<8: continue
    mm=re.match(r"^(.*)_(\d+)$",p[0])
    if mm:
        c,i=mm.group(1),int(mm.group(2))
        if i not in pfam[c] or float(p[6])<pfam[c][i][1]: pfam[c][i]=(p[3],float(p[6]))
dufpos=collections.defaultdict(list)
for l in open(f"{D}/duf2924_all.domtbl",errors="replace"):
    if l.startswith("#"): continue
    p=l.split()
    if len(p)<5: continue
    mm=re.match(r"^(.*)_(\d+)$",p[0].split("|")[-1])
    if mm: dufpos[mm.group(1)].append(int(mm.group(2)))
sel=json.load(open("fig/supp_loci.json"))            # Supp1 已挑好的 4 个物种
want={d["contig"] for d in sel}
co=coords([f"{D}/dufctx/chunk_{{}}.faa"],want)
duf_loci=[]
for d in sel:
    c=d["contig"]; cc=co.get(c,{})
    if not cc: continue
    idx=d["duf_idx"]
    genes=[]
    for i in sorted(cc):
        if abs(i-idx)>4: continue
        st,en,sd=cc[i]
        pf=pfam.get(c,{}).get(i,(None,))[0]
        if i==idx: pf="DUF2924"
        genes.append({"i":i,"start":st,"end":en,"strand":sd,"pfam":pf})
    if genes: duf_loci.append({"organism":d["organism"],"contig":c,"duf_idx":idx,"genes":genes})

# ── 对照侧：挑丝氨酸重组酶、不带 DUF2924 的位点 ──
ctrl=[l.split() for l in open(f"{R}/ctrl_recombinase_pos.tsv") if len(l.split())>=3]
ser=[(c,int(i)) for c,i,t in ctrl if t=="serine"]
met=collections.defaultdict(dict)
for l in open(f"{R}/ctrl_methyl.domtbl",errors="replace"):
    if l.startswith("#"): continue
    p=l.split()
    mm=re.match(r"^(.*)_(\d+)$",p[0])
    if mm: met[mm.group(1)].setdefault(int(mm.group(2)),p[3])
cwant={c for c,_ in ser}
cco=coords(sorted(glob.glob(f"{R}/ctrlctx/chunk_*.faa")),cwant)
ctrl_loci=[]
for c,idx in ser:
    cc=cco.get(c,{})
    if len(cc)<5: continue
    genes=[]
    for i in sorted(cc):
        if abs(i-idx)>4: continue
        st,en,sd=cc[i]
        pf="serine recombinase" if i==idx else met.get(c,{}).get(i)
        genes.append({"i":i,"start":st,"end":en,"strand":sd,"pfam":pf})
    if len(genes)>=7:
        ctrl_loci.append({"organism":"control locus (no DUF2924)","contig":c,
                          "duf_idx":idx,"genes":genes})
    if len(ctrl_loci)>=2: break
json.dump({"duf":duf_loci,"ctrl":ctrl_loci},open("fig/supp2_loci.json","w"),
          ensure_ascii=False,indent=1)
print(f"  DUF2924 位点 {len(duf_loci)} 个（带坐标与链向）")
for d in duf_loci: print(f"   {d['organism']:<34} {len(d['genes'])} 基因")
print(f"  对照位点 {len(ctrl_loci)} 个")
for d in ctrl_loci: print(f"   {d['contig']:<34} {len(d['genes'])} 基因 "
                          f"recombinase@{d['duf_idx']}")
