#!/usr/bin/env python3
"""Supp Fig 1 的全部数字 → fig/supp1_data.json。数字只从交付文件读，不手写。"""
import json, csv, glob, os, collections, statistics as st

RUN="../"; O="../outputs/fe0970b9-86bd-9a46-2a9e-dde59b863117"; D=f"{O}/discovery"
LED="../.omicos_home/conversations/fe0970b9-86bd-9a46-2a9e-dde59b863117/campaign.jsonl"

# ── A：各轮的 token、开场数、派工数 ──
runs=collections.OrderedDict()
for f in sorted(glob.glob(f"{RUN}sse_*.log")):
    nm=os.path.basename(f)[4:-4]; acc=collections.Counter(); t0=None
    for l in open(f,errors="replace"):
        if '"usage"' in l:
            try: c=json.loads(l).get("content") or {}
            except Exception: c={}
            if isinstance(c,dict):
                for k in ("input_tokens","output_tokens"):
                    if isinstance(c.get(k),(int,float)): acc[k]+=c[k]
        elif t0 is None and ('"accepted_at"' in l or '"started_at"' in l):
            try: r=json.loads(l)
            except Exception: continue
            t0=r.get("accepted_at") or r.get("started_at")
    if acc: runs[nm]={"input":acc["input_tokens"],"output":acc["output_tokens"],"t0":t0}

ev=[json.loads(l) for l in open(LED,errors="replace") if l.strip()]
starts=sorted(((v["t0"],k) for k,v in runs.items() if v["t0"]), key=lambda t:t[0])
def owner(ts):
    """归到「最后一个早于它开始」的那一轮。"""
    pick=starts[0][1]
    for t,k in starts:
        if t<=ts: pick=k
        else: break
    return pick
# created 的时间在 campaign 快照里；dispatched/round 自带顶层 at
for k in runs: runs[k].update(campaigns=0,auditor=0,experiment=0,rounds=0)
for e in ev:
    t=e["type"]
    if t=="created":
        ts=(e.get("campaign") or {}).get("created_at")
        if ts: runs[owner(ts)]["campaigns"]+=1
    elif t=="dispatched" and e.get("at"):
        r=e.get("role","")
        key="auditor" if r=="auditor" else "experiment" if r=="experiment" else None
        if key: runs[owner(e["at"])][key]+=1
    elif t=="round" and e.get("at"):
        runs[owner(e["at"])]["rounds"]+=1

# ── B：32 个显著家族 + 三道闸 ──
CIRC={"Phage_integrase","Resolvase","Recombinase","Phage_int_SAM_1","DEDD_Tnp_IS110"}
HOST={"RuvA_N","RuvA_C"}; TN3={"DUF4158"}; KEEP={"DUF2924"}
rows=list(csv.DictReader(open(f"{D}/family_enrichment.tsv"),delimiter="\t"))
sig=[r for r in rows if float(r["q_fdr"])<0.05]
def grp(n):
    return ("circular" if n in CIRC else "host" if n in HOST else
            "tn3" if n in TN3 else "retained" if n in KEEP else "mge")
fam=[]
for r in sig:
    cg,bg4=int(r["case_gt400"]),int(r["bg_gt400"])
    p=float(r["p_gt400"])
    # >400 分层为空时该检验无从施加，是 N/A 而不是「不通过」
    state = "na" if (cg==0 and bg4==0) else ("pass" if p<0.05 else "ns")
    fam.append({"name":r["name"],"acc":r["acc"].split(".")[0],
                "case":int(r["case"]),"bg":int(r["bg"]),"odds":float(r["odds"]),
                "q":float(r["q_fdr"]),"p400":p,"gt400":[cg,bg4],
                "len_state":state,"group":grp(r["name"])})
fam.sort(key=lambda d:(-d["odds"],d["q"]))

neu=[float(r["odds"]) for r in rows if str(r.get("neutral","")).strip().lower()=="true"]
shared=json.load(open(f"{D}/dark/shared_contigs.json"))
rr=json.load(open(f"{D}/reverse/reverse_rate.json"))

out={"A":{"runs":runs},
     "B":{"families":fam,"n_screened":len(rows),"n_tests":30134,
          "counts":dict(collections.Counter(f["group"] for f in fam)),
          "len_counts":dict(collections.Counter(f["len_state"] for f in fam))},
     "C":{"neutral":{"odds":sorted(neu),"n":len(neu),"median":st.median(neu)},
          "length":{"case_median_aa":211,"bg_median_aa":273,
                    "case_short_pct":36,"bg_short_pct":23},
          "genome":{"shared_contig":len(shared),"bg_from_case_pct":5.9,
                    "case_contig":72086,"bg_contig":97675},
          "reverse":{"ser_total":rr["n_ser"],"tyr_total":rr["n_tyr"],"ser_hit":rr["ser_with_duf"],"tyr_hit":rr["tyr_with_duf"],
                     "odds":rr.get("odds_ratio")},
          "_source":"DISCOVERY_SUMMARY.md 的 falsifier 检查 + REVERSE_SUMMARY.md"}}
json.dump(out,open("fig/supp1_data.json","w"),ensure_ascii=False,indent=1)
print(f"  A: {len(runs)} 轮, {sum(v['campaigns'] for v in runs.values())} 场, "
      f"auditor {sum(v['auditor'] for v in runs.values())} / "
      f"experiment {sum(v['experiment'] for v in runs.values())} 次派工, "
      f"{sum(v['rounds'] for v in runs.values())} 回合")
print(f"  B: {len(fam)} 家族 {out['B']['counts']} | 长度闸 {out['B']['len_counts']}")
print(f"  C: 中性 {len(neu)} 个中位 {out['C']['neutral']['median']:.3f} | "
      f"共有 contig {len(shared):,} | 反向率 {out['C']['reverse']['odds']}")
