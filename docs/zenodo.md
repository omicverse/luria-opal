# Working files deposited separately

About 1.3 GB of intermediate files are too large for this repository and are deposited at Zenodo.
They are not needed to check any reported statistic; `tables/` holds everything the numbers rest on.

| file set | size | what it is |
|---|---|---|
| `v2/*.domtbl`, `v2/*.log` | ~700 MB | per-chunk hmmsearch output for the 10,000 case and 10,000 background proteins |
| `reverse/*_contigs.faa` | ~155 MB | contig protein sets for the control recombinase loci |
| `reverse/recombinase_contigs.txt` | 22 MB | contig list for the whole-proteome recombinase scan |
| `dark/dufctx/chunk_*.faa` | ~90 MB | neighbourhood proteins around every DUF2924 hit, with Prodigal coordinates in the headers |
| `ctx/`, `trna/`, structure predictions | ~300 MB | gene-context extracts, tRNA windows, the 24 predicted structures and foldseek output |

DOI: [10.5281/zenodo.23251666](https://doi.org/10.5281/zenodo.23251666)

The concept DOI [10.5281/zenodo.23251665](https://doi.org/10.5281/zenodo.23251665) always resolves
to the latest version. Verify downloads against `SHA256SUMS` in the deposit; the checksums Zenodo
shows on its own pages are MD5.
