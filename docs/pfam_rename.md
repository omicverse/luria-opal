# Proposed rename: PF11149, DUF2924 to OPAL

**Short name** `OPAL`
**Full name** Operon-coupled accessory of large serine recombinases
**Synonym to retain** DUF2924

## Current state of the entry

As of October 2026, PF11149 is typed as a family, belongs to **no clan**, and carries **no attached
literature**. A search of Pfam returns 6,558 of its 30,134 families under a label of unknown
function; this is one of them.

## Basis for the rename

The name records a measured genomic association, not a molecular function.

- 14,088 of 524,414 serine recombinase loci across the proteome carry the domain within five genes,
  against 665 of 1,403,550 tyrosine integrase loci: odds ratio 58.2, Fisher p below double-precision
  resolution.
- 12,790 of 13,494 pairable loci (95%) place the recombinase at exactly one gene.
- All 168 loci with gene coordinates are on the same strand as their recombinase; median intergenic
  distance −3 bp; 76% overlap; 95% within 20 bp.
- At 89 of those 168 loci the junction is ATGA, the stop codon TGA of the upstream gene sharing
  three bases with the start codon ATG of the downstream gene. The upstream cistron therefore
  terminates on an **opal** codon that is also its partner's start, which is what the acronym
  records.
- Overlap length takes only two of three values modulo 3 across the 144 overlapping loci
  (χ² = 123, p = 2 × 10⁻²⁷), so the overlap preserves a fixed reading-frame relationship.

## What the name does not claim

The fold is undetermined: none of 24 predicted structures reaches a confident foldseek match
(best E = 0.074). Co-folding has no discriminative power in this system, calibrated against
interactions that must exist. The GANTC methylation-reading mechanism implied by the RAMA name was
tested with a falsifier declared beforehand and is not supported.

## Relationship to RAMA (PF18755)

PF11149 partly matches RAMA, as noted in passing by Tan et al. 2024 (PMID 39106433). Splitting the
183 family members by whether they match the RAMA model gives 82 and 101 members that are equally
likely to carry a recombinase within two genes (77/82 against 98/101, Fisher p = 0.47), so the
association is a property of the whole family.
