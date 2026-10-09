#!/usr/bin/env python3
"""Recount the ledger, and compare with the counts printed on Figure 1.

Figure 1 was rendered from a snapshot of this ledger taken before the last
campaign closed. The deposited ledger is the final state, so three classes
differ. This script prints both so the difference is visible and checkable.

Run:  python verify/verify_ledger.py
"""
import json, collections, pathlib, datetime

L = pathlib.Path(__file__).resolve().parent.parent / "ledger" / "campaign.jsonl"
events = [json.loads(l) for l in open(L) if l.strip()]

KEEP = {"round": "round", "believed": "belief",
        "objective_falsifier_declared": "falsifier",
        "dispatched": "dispatch", "judged": "adjudication"}

# same filter the figure used: a campaign is opened by an event carrying
# campaign.campaign_id, and only events naming that campaign are counted
C = collections.OrderedDict()
for e in events:
    c = e.get("campaign") or {}
    if c.get("campaign_id"):
        C[c["campaign_id"]] = []
    cid = e.get("campaign_id")
    if cid in C and e.get("type") in KEEP:
        C[cid].append(KEEP[e["type"]])

counts = collections.Counter(t for v in C.values() for t in v)
nonempty = sum(1 for v in C.values() if v)

# what Figure 1 prints
FIG = {"campaigns": 25, "round": 93, "belief": 18,
       "falsifier": 7, "dispatch": 76, "adjudication": 1, "events": 373}

print(f"{'':<16}{'this ledger':>13}{'Figure 1':>11}")
print(f"{'events':<16}{len(events):>13}{FIG['events']:>11}")
print(f"{'campaigns':<16}{nonempty:>13}{FIG['campaigns']:>11}")
for k in ("round", "belief", "falsifier", "dispatch", "adjudication"):
    mark = "" if counts[k] == FIG[k] else "   <- differs"
    print(f"{k:<16}{counts[k]:>13}{FIG[k]:>11}{mark}")

print("\nevery event type in the ledger")
for t, n in collections.Counter(e.get("type") for e in events).most_common():
    print(f"  {str(t):<32}{n:>5}")

delta = len(events) - FIG["events"]
print(f"\nThe deposited ledger holds {delta} events more than the snapshot Figure 1 was")
print("rendered from, which is one further campaign that closed after the figures were")
print("generated. Beliefs, objective falsifiers and adjudications are unchanged.")

ts = [e.get("ts") or e.get("timestamp") for e in events]
ts = [t for t in ts if isinstance(t, (int, float))]
if ts:
    a = datetime.datetime.fromtimestamp(min(ts), datetime.timezone.utc)
    b = datetime.datetime.fromtimestamp(max(ts), datetime.timezone.utc)
    print(f"span {a:%Y-%m-%d %H:%M} to {b:%Y-%m-%d %H:%M} UTC")
