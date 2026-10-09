# luria-opal

Audit material for an autonomous genome-mining campaign.

An autonomous multi-agent system, **Luria**, was given a brief to find a new programmable system.
It surveyed 618,638,921 proteins through five independent entry points and returned nothing on all
five. A by-product of the same run nominated one Pfam family of unknown function, DUF2924
(PF11149), which we name **OPAL**. Counted outward from every recombinase in the proteome, OPAL
lies beside serine recombinases 58 times more often than beside tyrosine integrases.

**This repository exists so that the analysis can be checked, not so that it can be read.** It
holds the ledger the agents wrote, the tables their claims rest on, and scripts that recompute the
headline statistics from those tables. The manuscript and figures are not here.

## Check it yourself

```bash
pip install scipy statsmodels
python verify/verify_ledger.py     # recount the ledger, compare with Figure 1
python verify/verify_screen.py     # recompute the enrichment screen, ~3 min
python verify/verify_reverse.py    # recompute the 58-fold rate, ~1 min, 2 GB RAM
```

Each script reads only from `tables/` or `ledger/` and prints its own result beside the value the
manuscript reports. Three discrepancies are expected and are printed by the scripts themselves:

| what | why |
|---|---|
| ledger holds 384 events / 26 campaigns, Figure 1 says 373 / 25 | the figures were rendered from a snapshot taken before the last campaign closed; beliefs, falsifiers and adjudications are unchanged |
| the screen deliverable says BH over 30,134 tests | the q values in `tables/family_enrichment.tsv` match a correction over the 4,857 families actually testable; the deliverable is wrong on this point and the table is right |
| the screen is one-sided | it asks whether a family is over-represented beside integrases, not whether it differs in either direction; two-sided gives 46 families at q < 0.05 rather than 32 |

## Layout

```
ledger/campaign.jsonl        the append-only ledger, 384 events over 26 campaigns
deliverables/                the agents' own campaign write-ups, verbatim, in the original Chinese
tables/                      the numbers every reported statistic rests on
  family_enrichment.tsv        4,857 families: case, background, odds, p, q, length-stratified gate
  reverse_rate.json            the serine-versus-tyrosine comparison
  duf2924_pos.json             gene index of every DUF2924 protein, by contig
  serine_pos.json              gene index of every serine recombinase, by contig
  tyrosine_pos.json            gene index of every tyrosine integrase, by contig
  panels/                      the per-panel tables behind each figure mark
provenance/                  the scripts that produced tables/panels from the run outputs
verify/                      independent recomputation of the headline statistics
```

## What is not here

- **The Luria harness.** The architecture is not published. What the manuscript claims about it is
  ledger behaviour, which `ledger/` and `verify/verify_ledger.py` let you check directly.
- **The raw agent transcript.** It contains API credentials and full agent prompts.
- **The manuscript and figures.**
- **The screen working files**, about 1.3 GB: per-chunk HMMER output, FASTA shards, provirus calls,
  predicted structures. These are deposited at Zenodo; see `docs/zenodo.md`.

## Key parameters

| | |
|---|---|
| Proteins surveyed | 618,638,921 (GTDB r232, viral RefSeq, phage genomes) |
| Gene calls | Prodigal v2.6.3, coordinates from FASTA headers |
| Domain assignment | HMMER 3.4, `hmmsearch --cut_ga` against Pfam-A.hmm (30,134 models) |
| Screen sampling | 10,000 case and 10,000 background proteins, seed 42 |
| Clustering | MMseqs2 `easy-cluster --min-seq-id 0.30 -c 0.5 --cov-mode 1` |
| Model behind Luria | DeepSeek-4-Pro |

Full parameters, including every tool version, are in the manuscript's STAR Methods.

## Licence

Code under MIT. Tables, ledger and deliverables under CC BY 4.0.
