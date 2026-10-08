# Small RNA-Seq kits: adapters and trimming (NEXTFLEX v3 vs Revvity NEXTFLEX v4)

Compiled 2026-10-06.

## Summary

The 3' adapter core is the same in both kits, but the kits should not be analyzed identically. v3 has 4 random bases on each end of the insert that must be trimmed. v4 appears not to.

| | NEXTFLEX v3 | Revvity NEXTFLEX v4 (with UDIs) |
|---|---|---|
| 3' adapter (to clip) | `TGGAATTCTCGGGTGCCAAGG` | `TGGAATTCTCGGGTGCCAAGG` (same) |
| Random bases next to insert | Yes, 4 nt on each side (3' and 5' adapters) | Not in the published adapter sequences; a Biostars commenter says v4 doesn't use randomized bases |
| Indexing | Single index (Bioo barcodes) | Unique dual indexes (P5 and P7) |

## v4 sequences (from the manual)

- 3' adapter: `5' rApp /TGGAATTCTCGGGTGCCAAGG/ 3SpC3/`
- 5' adapter: `5' UCUUUCCCUACACGACGCUCUUCCGAUCU`
- RT primer: `5' CCTTGGCACCCGAGAATTCCA`
- P7 index primer: `5' CAAGCAGAAGACGGCATACGAGAT [8-nt index] GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA`
- P5 index primer: `5' AATGATACGGCGACCACCGAGATCTACAC [8-nt index] ACACTCTTTCCCTACACGACGCTCTTCCGATCT`

The manual gives no explicit trimming guidance.

## Practical consequences

- **v3**
  1. Clip `TGGAATTCTCGGGTGCCAAGG`.
  2. Trim 4 nt from both ends of the clipped read, e.g. `cutadapt -a TGGAATTCTCGGGTGCCAAGG -u 4 -u -4 -m 15`. Without this step, miRNA alignments and counts are distorted.
  3. Optionally use the 4-nt tags as UMIs. Ligation bias makes them imperfectly random (each miRNA has its own ligation preference).
- **v4**: Clip the same 3' adapter. Do not blindly trim 4 nt from each end, because that would cut real miRNA bases. Confirm by checking the first and last bases of your reads, or ask Revvity support.
- **Mixed cohorts**: Process each kit as its own batch with its own trimming parameters. Kit type is a likely batch effect, so include it as a covariate or avoid comparing across kits directly.

## Caveats

- Revvity's "Nextflex Small RNA trimming instructions" PDF returned a 404, so the v4 statements rest on the manual and a forum comment. Check current Revvity documentation before finalizing.
- The v4 manual text retrieved was a summary, so the original sequence table was not seen.

## Sources

- [UC Davis DNA Technologies Core: trimming FAQ](https://dnatech.ucdavis.edu/faqs/how-should-the-mirnasmall-rna-data-be-trimmed/)
- [NEXTFLEX Small RNA-Seq Kit v4 with UDIs manual](https://resources.revvity.com/pdfs/man-nova-5132-XX-mn-nextflex-small-rna-seq-kit-v4-with-udis.pdf)
- [Biostars: Extracting UMIs from NEXTFLEX Small RNA-Seq reads](https://biostars.org/p/9597331)
