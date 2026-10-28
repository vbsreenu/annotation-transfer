# transfer_annotation.py

Transfers genome annotation from an annotated **reference** (GenBank or GFF3) to an
unannotated **target** genome (FASTA) through orthologous (syntenic) regions, then reports
which genes transferred, which failed and how the target genes differ from the reference.

Written for viral, organelle, plasmid and bacterial genomes (tested up to 4 Mb).

## Requirements

Python ≥ 3.8 plus three packages:

```bash
pip install biopython numpy matplotlib
```

It has been tested with Python 3.9 / Biopython 1.79 and Python 3.13 / Biopython 1.88, and gives identical results on both. It needs no external aligners (BLAST, minimap2 and so on).

## Usage

```bash
# GenBank reference
python transfer_annotation.py -r reference.gb -t new_strain.fasta -o results

# GFF3 reference (sequence in a ##FASTA section, or given separately)
python transfer_annotation.py -r reference.gff3 --ref-fasta reference.fasta \
       -t new_strain.fasta -o results -p strainX

# worked example (Arabidopsis chloroplast -> simulated rearranged genome)
python transfer_annotation.py -r example/reference_NC_000932.gb \
       -t example/target_simulated.fasta -o example/output -p chloroplast_example
```

Useful options (`-h` lists all of them):

| option | default | meaning |
|---|---|---|
| `--min-identity` | 0.5 | minimum nucleotide identity for a feature to count as transferred |
| `--min-coverage` | 0.5 | minimum fraction of the reference feature that must align |
| `--transl-table` | 1 | genetic code for CDS without `/transl_table` (11 for bacteria) |
| `--topology` | auto | target topology; auto = circular for a complete genome of a circular reference |
| `--exclude-types` | – | e.g. `variation,misc_difference` |
| `-k` | auto | anchor k-mer size; try a smaller value (e.g. 11) for very divergent genomes |
| `--max-gap` | 30000 | largest anchor-free gap bridged inside one block (e.g. inverted repeats) |
| `--gff-fasta` | off | include the target sequence in the GFF3 |

## Outputs (`<outdir>/<prefix>…`)

| file | content |
|---|---|
| `.gbk`, `.gff3` | annotated target genome. Every feature carries its reference location, identity and coverage. CDS also carry a new `/translation` and any problems found |
| `_report.html` | self-contained report with summary cards, the genome map, failed genes and a sortable, filterable table of gene differences |
| `_genome_map.png` / `.svg` | genes drawn as arrows along the reference and target. Forward-strand genes sit above the line and point right, reverse-strand genes sit below and point left. Ribbons link the orthologous blocks |
| `_genes.tsv` | one row per reference gene: status, reason for failure, locations, length change, nt identity, coverage, mismatches, indels, aa identity, non-synonymous and synonymous codon counts, aa changes, CDS flags |
| `_features.tsv` | the same metrics for every feature (gene, mRNA, CDS, exon, mat_peptide, …) plus the reference and target proteins |
| `_variants.tsv` | SNPs, insertions, deletions and complex or unaligned regions in reference coordinates, with target position, gene and coding effect (e.g. `S:D614G`, `frameshift`) |
| `_blocks.tsv` | orthologous blocks (reference ↔ target, orientation, identity) |
| `_summary.txt` | plain-text summary (also printed to the terminal) |

**Success and failure.** A gene counts as transferred when its product (CDS or RNA) aligns with at
least `--min-coverage` and `--min-identity`. Genes that fail are listed with a reason: absent or too
divergent, low coverage, or low identity.

**CDS flags.** The possible flags are `frameshift`, `premature_stop:aaN`, `start_codon_lost`,
`stop_codon_lost`, and `5'/3'_truncated` (the gene runs off a contig end or a rearrangement breakpoint;
written as a partial `<`/`>` feature). If bases are deleted at the very end of a gene but the genome
continues, this is reported as a terminal deletion, not as truncation.

## How it works

1. **Anchors.** k-mers that occur exactly once in each genome (counting both strands) are matched,
   and matches on the same diagonal are merged.
2. **Orthologous blocks.** Anchors are chained by collinear dynamic programming, as in minimap2.
   Inversions and translocations become separate blocks. Blocks are stitched across the origin of
   circular genomes. Repeats such as inverted or terminal repeats have no unique anchors, so they
   are bridged and placed by their syntenic flanks. That way each copy maps to its orthologous copy.
3. **Base-level alignment.** Gaps between anchors are filled with Needleman–Wunsch alignment
   (Biopython), using recursive re-anchoring for long gaps. Block ends are extended until the
   alignment score drops. Stretches with no detectable homology stay unaligned instead of being
   forced into an alignment.
4. **Lift-over.** Each reference locus (a gene with its mRNA, CDS, exon and peptide children) is
   projected through the block where its product aligns best. All features of a gene therefore
   move together.
5. **Comparison.** For every feature it measures identity, coverage and indels. For each CDS it also
   translates the protein, checks the start and stop codons, internal stops and frameshifts, and
   records amino-acid changes and synonymous/non-synonymous codon changes.

## Tests

```bash
bash tests/run_tests.sh
```

The test script simulates targets from three real NCBI genomes (phiX174, the Arabidopsis chloroplast
and a Yersinia plasmid). The simulations add SNPs, indels, a gene deletion, a novel insertion, an
inversion, rotation, reverse-complementing, splitting into contigs and up to 30% divergence. Every
transferred feature location is then checked against the simulated truth. All locations are correct
except the expected cases: deleted genes, and genes cut by contig ends or inversion breakpoints,
which are reported as partial. In the chloroplast simulation, SNP calls matched the truth with 99.96% recall
and 100% precision.

## Limitations

- Orthology is detected at nucleotide level. It works down to about 65–70% nt identity. For more
  distant genomes, use a protein-based approach.
- CDS with frameshifts or a lost stop codon are reported, not corrected. The ORF is not
  re-predicted.
- Genes that exist only in the target are not predicted. Target regions with no counterpart in the
  reference (≥ 500 bp) are listed in the report.
- GFF3 cannot record the exon order of trans-spliced genes. Such genes are lifted in coordinate
  order.
- Memory use is about 1.3 GB for a 4 Mb bacterial genome. The tool is not intended for large
  eukaryotic chromosomes.
