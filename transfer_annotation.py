#!/usr/bin/env python3
"""
transfer_annotation.py
======================

Transfer genome annotation from an annotated reference genome (GenBank or
GFF3) to an unannotated target genome (FASTA) through orthologous
(syntenic) regions, then report what was transferred, what failed and how
the target genes differ from the reference genes.

Method
------
1. Anchors   - k-mers that are unique in both genomes (both strands) are
               matched; exact matches on the same diagonal are merged.
2. Synteny   - anchors are chained with a collinear dynamic-programming
               chainer (minimap2-like).  Each chain is an orthologous,
               syntenic block; inversions / rearrangements give separate
               blocks and circular genomes are stitched across the origin.
               Repeats (e.g. inverted repeats) are resolved by their
               syntenic context.
3. Alignment - every block becomes a base-level alignment: anchors are
               exact matches, gaps between anchors are filled with
               Needleman-Wunsch (Biopython PairwiseAligner, recursive
               re-anchoring for long gaps), block ends are extended with an
               X-drop style semi-global alignment.
4. Lift-over - each reference locus (a gene together with its mRNA / CDS /
               exon / mat_peptide children, or a stand-alone feature) is
               projected through the block that covers it best, so all
               features of a gene are moved consistently.
5. Checks    - nucleotide identity, coverage, substitutions and indels;
               for CDS also translation, start/stop codons, premature stop
               codons, frameshifts, amino-acid changes and synonymous /
               non-synonymous codon changes.

Outputs (written to --outdir, file names start with --prefix)
-------
  <prefix>.gbk                  annotated target genome (GenBank)
  <prefix>.gff3                 annotated target genome (GFF3)
  <prefix>_genes.tsv            one row per reference gene: status + differences
  <prefix>_features.tsv         one row per reference feature
  <prefix>_variants.tsv         nucleotide differences reference -> target
  <prefix>_blocks.tsv           syntenic (orthologous) blocks
  <prefix>_genome_map.png/.svg  genes drawn as arrows along both genomes
  <prefix>_report.html          self-contained HTML report
  <prefix>_summary.txt          plain-text summary

Requirements
------------
  Python >= 3.8, biopython >= 1.79, numpy, matplotlib
      pip install biopython numpy matplotlib

Usage
-----
  python transfer_annotation.py -r reference.gb -t target.fasta -o out_dir
  python transfer_annotation.py -r reference.gff3 --ref-fasta reference.fasta \\
         -t target.fasta -o out_dir -p sample1
  python transfer_annotation.py -h        (all options)
"""

import argparse
import base64
import datetime
import html
import io
import logging
import math
import os
import re
import sys
import time
import urllib.parse
import warnings
from collections import Counter, OrderedDict, defaultdict
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

try:
    from Bio import SeqIO, BiopythonWarning
    from Bio.Align import PairwiseAligner, substitution_matrices
    from Bio.Data import CodonTable
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import (SeqFeature, CompoundLocation, BeforePosition,
                                AfterPosition, ExactPosition)
    try:
        from Bio.SeqFeature import SimpleLocation
    except ImportError:                      # Biopython < 1.80
        from Bio.SeqFeature import FeatureLocation as SimpleLocation
except ImportError:                          # pragma: no cover
    sys.exit("ERROR: Biopython is required:  pip install biopython numpy matplotlib")

warnings.simplefilter("ignore", BiopythonWarning)

__version__ = "1.0.0"
LOG = logging.getLogger("transfer_annotation")

# ----------------------------------------------------------------------------
# Feature type vocabularies
# ----------------------------------------------------------------------------
SKIP_TYPES = {"source", "region", "chromosome", "contig", "supercontig",
              "databank_entry", "assembly_component"}
GENE_PARENT_TYPES = {"gene", "pseudogene", "ncRNA_gene", "protein_coding_gene"}
RNA_TYPES = {"mRNA", "tRNA", "rRNA", "ncRNA", "tmRNA", "misc_RNA", "precursor_RNA",
             "lnc_RNA", "lncRNA", "snRNA", "snoRNA", "scRNA", "miRNA", "transcript",
             "primary_transcript", "pseudogenic_transcript", "pseudogenic_tRNA",
             "pseudogenic_rRNA", "RNase_P_RNA", "RNase_MRP_RNA", "SRP_RNA",
             "antisense_RNA", "guide_RNA", "telomerase_RNA", "Y_RNA", "vault_RNA"}
CODING_TYPES = {"CDS"}
# a locus counts as a "gene" if it contains any of these
GENE_LOCUS_TYPES = GENE_PARENT_TYPES | RNA_TYPES | CODING_TYPES
# children that sit inside a CDS
PEPTIDE_TYPES = {"mat_peptide", "sig_peptide", "transit_peptide", "propeptide",
                 "mature_protein_region", "mature_protein_region_of_CDS",
                 "signal_peptide", "transit_peptide_region_of_CDS"}

# GenBank feature key -> Sequence Ontology term for GFF3 output
SO_TYPES = {
    "CDS": "CDS", "gene": "gene", "mRNA": "mRNA", "tRNA": "tRNA", "rRNA": "rRNA",
    "ncRNA": "ncRNA", "tmRNA": "tmRNA", "misc_RNA": "ncRNA", "precursor_RNA": "primary_transcript",
    "mat_peptide": "mature_protein_region_of_CDS", "sig_peptide": "signal_peptide",
    "transit_peptide": "transit_peptide", "propeptide": "propeptide",
    "5'UTR": "five_prime_UTR", "3'UTR": "three_prime_UTR", "exon": "exon", "intron": "intron",
    "misc_feature": "sequence_feature", "repeat_region": "repeat_region",
    "stem_loop": "stem_loop", "rep_origin": "origin_of_replication",
    "regulatory": "regulatory_region", "promoter": "promoter", "terminator": "terminator",
    "RBS": "ribosome_entry_site", "polyA_signal": "polyA_signal_sequence",
    "polyA_site": "polyA_site", "protein_bind": "protein_binding_site",
    "misc_binding": "binding_site", "primer_bind": "primer_binding_site",
    "mobile_element": "mobile_genetic_element", "LTR": "long_terminal_repeat",
    "variation": "sequence_variant", "misc_difference": "sequence_difference",
    "D-loop": "D_loop", "STS": "STS", "gap": "gap", "assembly_gap": "assembly_gap",
    "operon": "operon", "pseudogene": "pseudogene", "lnc_RNA": "lnc_RNA",
    "unsure": "sequence_uncertainty", "oriT": "origin_of_transfer",
}

# GFF3 attribute <-> GenBank qualifier
GFF2GB_ATTR = {"Note": "note", "Dbxref": "db_xref", "Alias": "gene_synonym",
               "Ontology_term": "note"}
GB2GFF_ATTR = {"note": "Note", "db_xref": "Dbxref"}


# ----------------------------------------------------------------------------
# Data classes
# ----------------------------------------------------------------------------
@dataclass
class SeqEntry:
    """A genome sequence (one record of a possibly multi-record genome)."""
    sid: str
    seq: str
    circular: bool = False
    description: str = ""
    molecule_type: str = "DNA"
    organism: str = ""

    @property
    def length(self):
        return len(self.seq)


@dataclass
class RefFeature:
    """A reference feature.  ``parts`` are (start, end, strand), 0-based,
    half-open, listed in biological (5'->3') order."""
    idx: int
    ftype: str
    seqid: str
    parts: List[Tuple[int, int, int]]
    qualifiers: "OrderedDict[str, List[str]]"
    fmt: str = "genbank"
    partial5: bool = False
    partial3: bool = False
    gff_id: Optional[str] = None
    gff_parents: List[str] = field(default_factory=list)
    gff_source: str = "."
    gff_score: str = "."
    locus_id: Optional[str] = None

    @property
    def strand(self):
        s = {p[2] for p in self.parts}
        return s.pop() if len(s) == 1 else 0

    @property
    def start(self):
        return min(p[0] for p in self.parts)

    @property
    def end(self):
        return max(p[1] for p in self.parts)

    @property
    def length(self):
        return sum(e - s for s, e, _ in self.parts)

    def q(self, key, default=None):
        v = self.qualifiers.get(key)
        return v[0] if v else default

    @property
    def name(self):
        for key in ("gene", "Name", "locus_tag", "standard_name", "product", "label"):
            v = self.q(key)
            if v:
                return v
        if self.gff_id:
            return self.gff_id
        return "%s_%d" % (self.ftype, self.idx)


@dataclass
class Locus:
    lid: str
    seqid: str
    features: List[RefFeature]
    primary: RefFeature = None
    is_gene: bool = False

    @property
    def start(self):
        return min(f.start for f in self.features)

    @property
    def end(self):
        return max(f.end for f in self.features)

    @property
    def name(self):
        for f in [self.primary] + self.features:
            for key in ("gene", "Name", "locus_tag", "standard_name", "product"):
                v = f.q(key)
                if v:
                    return v
        return self.primary.name

    @property
    def locus_tag(self):
        for f in self.features:
            v = f.q("locus_tag")
            if v:
                return v
        return ""

    @property
    def product(self):
        for f in self.core_features() + self.features:
            v = f.q("product")
            if v:
                return v
        return ""

    def core_features(self):
        """The functional products (CDS / RNA) or, if none, the primary feature."""
        core = [f for f in self.features if f.ftype in CODING_TYPES or f.ftype in RNA_TYPES]
        return core or [self.primary]

    @property
    def kind(self):
        types = [f.ftype for f in self.features]
        for t in ("CDS",):
            if t in types:
                return t
        for t in types:
            if t in RNA_TYPES and t != "mRNA":
                return t
        if "mRNA" in types:
            return "mRNA"
        return self.primary.ftype


# ----------------------------------------------------------------------------
# Small sequence helpers
# ----------------------------------------------------------------------------
_COMP = str.maketrans("ACGTRYKMSWBDHVNacgtrykmswbdhvn", "TGCAYRMKSWVHDBNtgcayrmkswvhdbn")


def revcomp(s):
    return s.translate(_COMP)[::-1]


def clean_seq(s):
    s = re.sub(r"[^A-Za-z]", "", s).upper()
    return s.replace("U", "T")


def extract_parts(seq, parts):
    """Concatenate (start, end, strand) parts of ``seq`` in the given order."""
    out = []
    for s, e, st in parts:
        sub = seq[s:e]
        out.append(revcomp(sub) if st < 0 else sub)
    return "".join(out)


def fmt_loc(seqid, parts, one_based=True):
    """Human readable location, e.g. NC_1:266..21555(+) or join(...)."""
    if not parts:
        return "-"
    segs = []
    for s, e, st in parts:
        segs.append("%d..%d" % (s + 1, e) if one_based else "%d..%d" % (s, e))
    strands = {st for _, _, st in parts}
    sym = {1: "+", -1: "-", 0: "."}
    if len(strands) == 1:
        body = segs[0] if len(segs) == 1 else "join(%s)" % ",".join(segs)
        return "%s:%s(%s)" % (seqid, body, sym[strands.pop()])
    body = ",".join("%s(%s)" % (sg, sym[st]) for sg, (_, _, st) in zip(segs, parts))
    return "%s:join(%s)" % (seqid, body)


# ----------------------------------------------------------------------------
# Input: target FASTA
# ----------------------------------------------------------------------------
def read_fasta(path, circular=False):
    entries = []
    seen = set()
    for rec in SeqIO.parse(path, "fasta"):
        sid = rec.id
        if sid in seen:
            raise SystemExit("ERROR: duplicate sequence id '%s' in %s" % (sid, path))
        seen.add(sid)
        seq = clean_seq(str(rec.seq))
        if not seq:
            LOG.warning("sequence %s in %s is empty - skipped", sid, path)
            continue
        desc = rec.description[len(rec.id):].strip()
        entries.append(SeqEntry(sid, seq, circular, desc))
    if not entries:
        raise SystemExit("ERROR: no sequences found in FASTA file %s" % path)
    return entries


# ----------------------------------------------------------------------------
# Input: reference GenBank
# ----------------------------------------------------------------------------
def _parts_from_location(loc):
    parts = []
    for p in loc.parts:
        st = p.strand if p.strand in (1, -1) else 0
        parts.append((int(p.start), int(p.end), st))
    return parts


def read_genbank(path):
    refs, feats = [], []
    n = 0
    for rec in SeqIO.parse(path, "genbank"):
        try:
            seq = clean_seq(str(rec.seq))
        except Exception:
            seq = ""
        if not seq:
            raise SystemExit("ERROR: GenBank record %s has no sequence (ORIGIN section). "
                             "A reference sequence is required." % rec.id)
        rid = rec.id if rec.id and rec.id != "<unknown id>" else rec.name
        organism = rec.annotations.get("organism", "")
        refs.append(SeqEntry(rid, seq, rec.annotations.get("topology", "linear") == "circular",
                             rec.description, rec.annotations.get("molecule_type", "DNA"), organism))
        for f in rec.features:
            if f.type in SKIP_TYPES:
                continue
            parts = _parts_from_location(f.location)
            if not parts:
                continue
            first, last = f.location.parts[0], f.location.parts[-1]
            if (first.strand or 1) >= 0:
                p5 = isinstance(first.start, BeforePosition)
            else:
                p5 = isinstance(first.end, AfterPosition)
            if (last.strand or 1) >= 0:
                p3 = isinstance(last.end, AfterPosition)
            else:
                p3 = isinstance(last.start, BeforePosition)
            quals = OrderedDict((k, list(v) if isinstance(v, (list, tuple)) else [str(v)])
                                for k, v in f.qualifiers.items())
            feats.append(RefFeature(n, f.type, rid, parts, quals, "genbank", p5, p3))
            n += 1
    if not refs:
        raise SystemExit("ERROR: no GenBank records found in %s" % path)
    return refs, feats


# ----------------------------------------------------------------------------
# Input: reference GFF3 (+ FASTA)
# ----------------------------------------------------------------------------
def _parse_gff_attrs(text):
    attrs = OrderedDict()
    for item in text.strip().strip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
        else:                                   # tolerate GTF-like "key value"
            bits = item.split(None, 1)
            k, v = bits[0], (bits[1].strip('"') if len(bits) > 1 else "")
        attrs[urllib.parse.unquote(k.strip())] = [urllib.parse.unquote(x) for x in v.split(",")]
    return attrs


def read_gff3(path, fasta_path=None):
    feats_by_key = OrderedDict()
    fasta_lines, in_fasta = [], False
    circular_ids, region_len = set(), {}
    anon = 0
    with open(path) as fh:
        for line in fh:
            if in_fasta:
                fasta_lines.append(line)
                continue
            line = line.rstrip("\r\n")
            if line.startswith("##FASTA"):
                in_fasta = True
                continue
            if line.startswith(">"):            # FASTA without ##FASTA directive
                in_fasta = True
                fasta_lines.append(line + "\n")
                continue
            if line.startswith("##sequence-region"):
                bits = line.split()
                if len(bits) >= 4:
                    region_len[bits[1]] = int(bits[3])
                continue
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) == 8:
                cols.append("")
            if len(cols) < 9:
                LOG.warning("GFF3: skipping malformed line: %s", line[:80])
                continue
            seqid, source, ftype, start, end, score, strand, phase, attr = cols[:9]
            attrs = _parse_gff_attrs(attr)
            if ftype in SKIP_TYPES:
                if attrs.get("Is_circular", ["false"])[0].lower() == "true":
                    circular_ids.add(seqid)
                continue
            s, e = int(start) - 1, int(end)
            st = {"+": 1, "-": -1}.get(strand, 0)
            fid = attrs.get("ID", [None])[0]
            key = (seqid, fid) if fid else ("__anon__", anon)
            if not fid:
                anon += 1
            if key in feats_by_key:
                f = feats_by_key[key]
                f.parts.append((s, e, st))
                f._phases.append(phase)
                continue
            quals = OrderedDict((k, v) for k, v in attrs.items() if k not in ("ID", "Parent"))
            f = RefFeature(len(feats_by_key), ftype, seqid, [(s, e, st)], quals, "gff3",
                           gff_id=fid, gff_parents=attrs.get("Parent", []),
                           gff_source=source, gff_score=score)
            f._phases = [phase]
            feats_by_key[key] = f
    feats = list(feats_by_key.values())

    # sequences
    seqs = OrderedDict()
    if fasta_path:
        for rec in SeqIO.parse(fasta_path, "fasta"):
            seqs[rec.id] = clean_seq(str(rec.seq))
    elif fasta_lines:
        for rec in SeqIO.parse(io.StringIO("".join(fasta_lines)), "fasta"):
            seqs[rec.id] = clean_seq(str(rec.seq))
    if not seqs:
        raise SystemExit("ERROR: the GFF3 file has no ##FASTA section; give the reference "
                         "sequence with --ref-fasta")
    used_ids = []
    for f in feats:
        if f.seqid not in used_ids:
            used_ids.append(f.seqid)
    missing = [s for s in used_ids if s not in seqs]
    if missing:
        raise SystemExit("ERROR: GFF3 sequence ids not found in reference FASTA: %s"
                         % ", ".join(missing[:5]))
    refs = [SeqEntry(sid, s, sid in circular_ids) for sid, s in seqs.items()]
    lens = {r.sid: r.length for r in refs}

    for f in feats:
        L = lens[f.seqid]
        fixed = []
        for s, e, st in f.parts:               # circular features may run past the end
            if e > L and s < L:
                fixed += [(s, L, st), (0, e - L, st)]
            elif s >= L:
                fixed.append((s - L, e - L, st))
            else:
                fixed.append((s, e, st))
        phases = list(f._phases) + ["."] * (len(fixed) - len(f._phases))
        # biological order: ascending (+) / descending (-), respecting origin wrap
        order = sorted(range(len(fixed)), key=lambda i: fixed[i][0])
        if len(order) > 1 and fixed[order[0]][0] == 0 and fixed[order[-1]][1] == L:
            gaps = [fixed[order[i + 1]][0] - fixed[order[i]][1] for i in range(len(order) - 1)]
            cut = int(np.argmax(gaps)) + 1           # rotate at the largest gap
            order = order[cut:] + order[:cut]
        if all(p[2] < 0 for p in fixed):
            order = order[::-1]
        f.parts = [fixed[i] for i in order]
        if len(f._phases) == len(fixed):
            phases = [phases[i] for i in order]
        if f.ftype == "CDS" and "codon_start" not in f.qualifiers and phases and phases[0] in "012":
            f.qualifiers["codon_start"] = [str(int(phases[0]) + 1)]
        strand = f.strand
        if "start_range" in f.qualifiers:
            if strand >= 0:
                f.partial5 = True
            else:
                f.partial3 = True
        if "end_range" in f.qualifiers:
            if strand >= 0:
                f.partial3 = True
            else:
                f.partial5 = True
        del f._phases
    return refs, feats


def detect_format(path):
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            if line.startswith("LOCUS"):
                return "genbank"
            if line.startswith("##gff-version") or line.count("\t") >= 8:
                return "gff3"
            if line.startswith(">"):
                raise SystemExit("ERROR: %s looks like FASTA; the reference must be GenBank or GFF3" % path)
            if line.startswith("#"):
                continue
            break
    ext = os.path.splitext(path)[1].lower()
    if ext in (".gb", ".gbk", ".genbank", ".gbff"):
        return "genbank"
    if ext in (".gff", ".gff3"):
        return "gff3"
    raise SystemExit("ERROR: cannot tell whether %s is GenBank or GFF3" % path)


# ----------------------------------------------------------------------------
# Group reference features into loci (gene + children)
# ----------------------------------------------------------------------------
def _choose_primary(features):
    for f in features:
        if f.ftype in GENE_PARENT_TYPES:
            return f
    cds = [f for f in features if f.ftype == "CDS"]
    if cds:
        return max(cds, key=lambda f: f.length)
    rna = [f for f in features if f.ftype in RNA_TYPES]
    if rna:
        return max(rna, key=lambda f: f.length)
    return features[0]


def build_loci(feats, fmt):
    loci = []
    if fmt == "gff3":
        by_id = {f.gff_id: f for f in feats if f.gff_id}

        def top(f):
            seen = set()
            while f.gff_parents and f.gff_parents[0] in by_id and f.gff_id not in seen:
                seen.add(f.gff_id)
                f = by_id[f.gff_parents[0]]
            return f
        groups = OrderedDict()
        for f in feats:
            t = top(f)
            key = t.gff_id or "feature_%d" % t.idx
            groups.setdefault((f.seqid, key), []).append(f)
        for (seqid, key), fs in groups.items():
            loci.append(Locus(key, seqid, fs))
    else:
        groups = OrderedDict()
        for f in feats:
            key = f.q("locus_tag") or f.q("gene")
            if key is None:
                key = "%s_%d" % (f.ftype, f.idx)
            groups.setdefault((f.seqid, key), []).append(f)
        for (seqid, key), fs in groups.items():
            # split identically named but separate genes into clusters
            fs_sorted = sorted(fs, key=lambda f: f.start)
            clusters, cur, cur_end = [], [], -1
            for f in fs_sorted:
                if cur and f.start >= cur_end:
                    clusters.append(cur)
                    cur, cur_end = [], -1
                cur.append(f)
                cur_end = max(cur_end, f.end)
            if cur:
                clusters.append(cur)
            for ci, cl in enumerate(clusters):
                lid = key if len(clusters) == 1 else "%s_%d" % (key, ci + 1)
                cl.sort(key=lambda f: f.idx)
                loci.append(Locus(lid, seqid, cl))
    seen = Counter()
    for loc in loci:
        seen[loc.lid] += 1
        if seen[loc.lid] > 1:                    # same id on two sequences
            loc.lid = "%s_%d" % (loc.lid, seen[loc.lid])
        loc.primary = _choose_primary(loc.features)
        loc.is_gene = any(f.ftype in GENE_LOCUS_TYPES for f in loc.features)
        for f in loc.features:
            f.locus_id = loc.lid
            f.locus_name = loc.name if loc.is_gene else None
    loci.sort(key=lambda l: (l.seqid, l.start))
    return loci


# ----------------------------------------------------------------------------
# k-mer anchors
# ----------------------------------------------------------------------------
_LUT = np.full(256, -1, dtype=np.int8)
for _c, _v in zip(b"ACGT", range(4)):
    _LUT[_c] = _v
    _LUT[ord(chr(_c).lower())] = _v
_EMPTY = (np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64))


def auto_kmer(total_len):
    k = int(math.ceil(math.log(max(total_len, 1000), 4))) + 4
    k = max(11, min(k, 25))
    return k + 1 if k % 2 == 0 else k        # odd k: no reverse-complement palindromes


def kmer_codes(seq, k):
    """2-bit encoded k-mers of ``seq`` (k <= 31) and their start positions,
    skipping k-mers that contain non-ACGT characters."""
    a = np.frombuffer(seq.encode("ascii", "replace"), dtype=np.uint8)
    n = a.size - k + 1
    if n <= 0:
        return np.zeros(0, dtype=np.uint64), np.zeros(0, dtype=np.int64)
    c = _LUT[a]
    bad = np.concatenate(([0], np.cumsum(c < 0)))
    ok = (bad[k:k + n] - bad[:n]) == 0
    cc = np.where(c < 0, 0, c).astype(np.uint64)
    h = np.zeros(n, dtype=np.uint64)
    two = np.uint64(2)
    for i in range(k):
        h = (h << two) | cc[i:i + n]
    pos = np.nonzero(ok)[0].astype(np.int64)
    return h[pos], pos


def _unique_codes(codes):
    """Return (sorted unique-once codes, index of their single occurrence)."""
    if codes.size == 0:
        return codes, np.zeros(0, dtype=np.int64)
    u, idx, cnt = np.unique(codes, return_index=True, return_counts=True)
    keep = cnt == 1
    return u[keep], idx[keep]


def merge_runs(group, I, J, k):
    """Merge k-mer matches that are consecutive on the same diagonal."""
    if I.size == 0:
        z = np.zeros(0, dtype=np.int64)
        return z, z, z, z
    D = J - I
    order = np.lexsort((I, D, group))
    g, i, j, d = group[order], I[order], J[order], D[order]
    brk = np.ones(i.size, dtype=bool)
    brk[1:] = (g[1:] != g[:-1]) | (d[1:] != d[:-1]) | (i[1:] != i[:-1] + 1)
    starts = np.nonzero(brk)[0]
    ends = np.r_[starts[1:], i.size]
    return g[starts], i[starts], j[starts], (ends - starts) + k - 1


def global_anchor_runs(refs, tgts, k):
    """Unique (in both genomes, both strands) shared k-mers -> merged runs.
    Returns arrays: group id (ref_idx, tgt_idx, orientation), ref pos,
    oriented target pos, run length."""
    t_codes, t_lab, t_pos = [], [], []
    for ti, t in enumerate(tgts):
        for flag, s in ((0, t.seq), (1, revcomp(t.seq))):
            h, p = kmer_codes(s, k)
            t_codes.append(h)
            t_pos.append(p)
            t_lab.append(np.full(h.size, ti * 2 + flag, dtype=np.int64))
    T = np.concatenate(t_codes)
    T_lab = np.concatenate(t_lab)
    T_pos = np.concatenate(t_pos)
    uT, iT = _unique_codes(T)

    rf_codes, rf_rec, rf_pos, rr_codes = [], [], [], []
    for ri, r in enumerate(refs):
        h, p = kmer_codes(r.seq, k)
        rf_codes.append(h)
        rf_pos.append(p)
        rf_rec.append(np.full(h.size, ri, dtype=np.int64))
        rr_codes.append(kmer_codes(revcomp(r.seq), k)[0])
    Rf = np.concatenate(rf_codes)
    Rf_rec = np.concatenate(rf_rec)
    Rf_pos = np.concatenate(rf_pos)
    R_all = np.concatenate([Rf] + rr_codes)
    uR, iR = _unique_codes(R_all)
    fwd = iR < Rf.size
    uR, iR = uR[fwd], iR[fwd]

    _, ia, ib = np.intersect1d(uR, uT, assume_unique=True, return_indices=True)
    ref_rec = Rf_rec[iR[ia]]
    ref_pos = Rf_pos[iR[ia]]
    lab = T_lab[iT[ib]]
    tgt_pos = T_pos[iT[ib]]
    group = ref_rec * (2 * len(tgts)) + lab          # lab = tgt_idx*2 + (1 if reverse)
    del T, T_lab, T_pos, R_all
    return merge_runs(group, ref_pos, tgt_pos, k), int(ia.size)


# ----------------------------------------------------------------------------
# Collinear chaining (minimap2-like dynamic programming)
# ----------------------------------------------------------------------------
def chain_runs(I, J, L, k, max_gap, min_score, lookback=64):
    """Chain exact-match runs (I = ref start, J = target start, L = length).
    Returns a list of (score, I, J, L) chains, best first."""
    n = I.size
    if n == 0:
        return []
    order = np.lexsort((J, I))
    I, J, L = I[order], J[order], L[order]
    Ie, Je = I + L, J + L
    f = L.astype(np.float64)
    pred = np.full(n, -1, dtype=np.int64)
    ck = 0.01 * k
    for a in range(1, n):
        lo = a - lookback if a > lookback else 0
        ia, ja = I[a], J[a]
        bI, bJ = I[lo:a], J[lo:a]
        di = ia - Ie[lo:a]
        dj = ja - Je[lo:a]
        ok = (bI < ia) & (bJ < ja) & (di <= max_gap) & (dj <= max_gap)
        if not ok.any():
            continue
        ov = np.maximum(np.maximum(-di, -dj), 0)
        gain = L[a] - ov
        dd = np.abs(di - dj)
        sc = f[lo:a] + gain - (ck * dd + 0.5 * np.log2(dd + 1.0))
        sc = np.where(ok & (gain > 0), sc, -np.inf)
        b = int(np.argmax(sc))
        if sc[b] > f[a]:
            f[a] = sc[b]
            pred[a] = lo + b
    chains = []
    used = np.zeros(n, dtype=bool)
    for a in np.argsort(-f, kind="stable"):
        if used[a]:
            continue
        path, x = [], a
        while x >= 0 and not used[x]:
            path.append(x)
            x = pred[x]
        used[path] = True
        score = f[a] - (f[x] if x >= 0 else 0.0)
        if score >= min_score:
            path.reverse()
            chains.append((float(score), I[path], J[path], L[path]))
    chains.sort(key=lambda c: -c[0])
    return chains


def runs_to_blocks(I, J, L):
    """Make chained runs non-overlapping -> list of (ref, tgt, length)."""
    blocks = []
    pe_i = pe_j = None
    for i, j, l in zip(I.tolist(), J.tolist(), L.tolist()):
        if blocks:
            sh = max(0, pe_i - i, pe_j - j)
            if sh:
                i, j, l = i + sh, j + sh, l - sh
            if l <= 0:
                continue
        blocks.append((i, j, l))
        pe_i, pe_j = i + l, j + l
    return blocks


# ----------------------------------------------------------------------------
# Base-level alignment between anchors
# ----------------------------------------------------------------------------
MATCH, MISMATCH, GAP_OPEN, GAP_EXTEND = 2, -3, -5, -2
MAX_CELLS = 12_000_000
EXT_LEN = 2000


def _nt_aligner(extend=False):
    a = PairwiseAligner()
    a.mode = "global"
    a.match_score = MATCH
    a.mismatch_score = MISMATCH
    a.open_gap_score = GAP_OPEN
    a.extend_gap_score = GAP_EXTEND
    if extend:                       # free overhang at the right end of either sequence
        a.target_right_open_gap_score = 0
        a.target_right_extend_gap_score = 0
        a.query_right_open_gap_score = 0
        a.query_right_extend_gap_score = 0
    return a


_GLOBAL = _nt_aligner()
_EXTEND = _nt_aligner(extend=True)


def _aligned_blocks(aln):
    """Gap-free aligned blocks [(a0, a1, b0, b1), ...] of a pairwise alignment
    (works with old and new Biopython alignment objects)."""
    coords = getattr(aln, "coordinates", None)
    if coords is not None:
        c = np.asarray(coords)
        out = []
        for x in range(c.shape[1] - 1):
            a0, b0 = int(c[0, x]), int(c[1, x])
            a1, b1 = int(c[0, x + 1]), int(c[1, x + 1])
            if a1 > a0 and b1 > b0:
                out.append((a0, a1, b0, b1))
        return out
    return [(int(a0), int(a1), int(b0), int(b1))
            for (a0, a1), (b0, b1) in zip(aln.aligned[0], aln.aligned[1])]


def _aligned_pairs(aln):
    rs, ts = [], []
    for a0, a1, b0, b1 in _aligned_blocks(aln):
        rs.append(np.arange(a0, a1, dtype=np.int64))
        ts.append(np.arange(b0, b1, dtype=np.int64))
    if not rs:
        return _EMPTY
    return np.concatenate(rs), np.concatenate(ts)


def _nw(r, t):
    try:
        return _aligned_pairs(_GLOBAL.align(r, t)[0])
    except (MemoryError, OverflowError):
        return _EMPTY


def local_chain(r, t, k):
    """Best collinear chain of k-mers unique within the two segments."""
    hr, pr = kmer_codes(r, k)
    ht, pt = kmer_codes(t, k)
    ur, ir = _unique_codes(hr)
    ut, it = _unique_codes(ht)
    _, a, b = np.intersect1d(ur, ut, assume_unique=True, return_indices=True)
    if a.size == 0:
        return None
    I, J = pr[ir[a]], pt[it[b]]
    _, I, J, L = merge_runs(np.zeros(I.size, dtype=np.int64), I, J, k)
    ch = chain_runs(I, J, L, k, max_gap=max(len(r), len(t)), min_score=3 * k)
    return ch[0] if ch else None


def _gap_cost(n):
    return 0 if n <= 0 else -(GAP_OPEN + GAP_EXTEND * (n - 1))


def _flush_align(r, t):
    """Short vs much longer segment: align the short one contiguously to the
    left or right end of the long one (the rest is a single indel), which
    keeps a long insertion/deletion in one piece instead of splitting it."""
    short_is_r = len(r) <= len(t)
    sh, lg = (r, t) if short_is_r else (t, r)
    m, M = len(sh), len(lg)
    w = min(M, m + 10)
    best = None
    for side in ("left", "right"):
        win = lg[:w] if side == "left" else lg[M - w:]
        sc = _GLOBAL.score(sh, win) if short_is_r else _GLOBAL.score(win, sh)
        if best is None or sc > best[0]:
            best = (sc, side, win)
    sc, side, win = best
    a, b = _nw(sh, win) if short_is_r else _nw(win, sh)
    if side == "right":
        if short_is_r:
            b = b + (M - w)
        else:
            a = a + (M - w)
    return a, b


def align_global(r, t, k, depth=0):
    """Alignment of two segments that are anchored at both ends (between two
    anchors).  Returns arrays of aligned (ref_offset, tgt_offset) pairs.
    Long segments are re-anchored recursively with smaller k-mers; segments
    without detectable homology are left unaligned."""
    nr, nt = len(r), len(t)
    if nr == 0 or nt == 0:
        return _EMPTY
    if min(nr, nt) <= 50 and max(nr, nt) >= 4 * min(nr, nt) + 20:
        return _flush_align(r, t)
    if nr * nt <= MAX_CELLS:
        try:
            aln = _GLOBAL.align(r, t)[0]
        except (MemoryError, OverflowError):
            return _EMPTY
        n = min(nr, nt)
        if n >= 100 and (aln.score + _gap_cost(abs(nr - nt))) / n <= -0.25:
            return _EMPTY                          # unrelated sequence (see _homologous)
        return _aligned_pairs(aln)
    if depth < 4:
        for kk in sorted({max(9, k - 2), max(9, k - 4)}, reverse=True):
            ch = local_chain(r, t, kk)
            if ch is None or ch[0] < max(3 * kk, 0.02 * min(nr, nt)):
                continue
            blocks = runs_to_blocks(ch[1], ch[2], ch[3])
            R, T = [], []
            pi = pj = 0
            for (i, j, l) in blocks:
                a, b = align_global(r[pi:i], t[pj:j], kk, depth + 1)
                R.append(pi + a)
                T.append(pj + b)
                R.append(np.arange(i, i + l, dtype=np.int64))
                T.append(np.arange(j, j + l, dtype=np.int64))
                pi, pj = i + l, j + l
            a, b = align_global(r[pi:], t[pj:], kk, depth + 1)
            R.append(pi + a)
            T.append(pj + b)
            return np.concatenate(R), np.concatenate(T)
    return _EMPTY                                  # no detectable homology: leave unaligned


def _trim_xdrop(r, t, pr, pt):
    """Keep the prefix (from the anchor outwards) with the best cumulative score."""
    if pr.size == 0:
        return _EMPTY
    rb = np.frombuffer(r.encode(), dtype=np.uint8)[pr]
    tb = np.frombuffer(t.encode(), dtype=np.uint8)[pt]
    col = np.where(rb == tb, MATCH, MISMATCH).astype(np.float64)
    gr = np.diff(pr, prepend=-1) - 1
    gt = np.diff(pt, prepend=-1) - 1
    pen = np.where(gr > 0, -GAP_OPEN - GAP_EXTEND * (gr - 1), 0) + \
        np.where(gt > 0, -GAP_OPEN - GAP_EXTEND * (gt - 1), 0)
    cum = np.cumsum(col - pen)
    best = int(np.argmax(cum))
    if cum[best] <= 0:
        return _EMPTY
    return pr[:best + 1], pt[:best + 1]


def extend_alignment(r, t):
    """Extend outwards from an anchor (sequences start next to the anchor).
    Semi-global alignment, then trimmed at the best cumulative score."""
    if not r or not t:
        return _EMPTY
    try:
        pr, pt = _aligned_pairs(_EXTEND.align(r, t)[0])
    except (MemoryError, OverflowError):
        return _EMPTY
    return _trim_xdrop(r, t, pr, pt)


def flank_alignment(r, t, k):
    """Extend from a block end over a flank of up to max-gap length.  If the
    homology reaches the end of the first EXT_LEN window (e.g. a long
    terminal repeat), the whole flank is aligned and then X-drop trimmed."""
    a, b = extend_alignment(r[:EXT_LEN], t[:EXT_LEN])
    if len(r) > EXT_LEN and len(t) > EXT_LEN and a.size and \
            a[-1] > 0.9 * EXT_LEN and b[-1] > 0.9 * EXT_LEN:
        ratio = len(r) / len(t)
        if 0.5 <= ratio <= 2:
            fa, fb = align_global(r, t, k)
            fa, fb = _trim_xdrop(r, t, fa, fb)
            if fa.size > a.size:
                return fa, fb
    return a, b


# ----------------------------------------------------------------------------
# Syntenic blocks (chains) with base-level maps
# ----------------------------------------------------------------------------
class Chain:
    """An orthologous block.  ``rmap[p - R0]`` is the oriented target
    coordinate aligned to reference position p (-1 = deleted in target).
    Coordinates can be 'virtual' (>= sequence length) for blocks that were
    stitched across the origin of a circular sequence."""

    def __init__(self, cid, ref, tgt, orient, R, T, score):
        self.cid = cid
        self.ref = ref
        self.tgt = tgt
        self.orient = orient
        self.score = score
        self.R0 = int(R[0])
        self.R1 = int(R[-1]) + 1
        self.rmap = np.full(self.R1 - self.R0, -1, dtype=np.int64)
        self.rmap[R - self.R0] = T
        m = np.nonzero(self.rmap >= 0)[0]
        self.ins_after = np.zeros(self.rmap.size, dtype=np.int64)
        if m.size > 1:
            self.ins_after[m[:-1]] = np.diff(self.rmap[m]) - 1
        self.n_aligned = int(m.size)
        self.T0, self.T1 = int(T[0]), int(T[-1]) + 1
        self.identity = None

    def lookup(self, s, e):
        out = np.full(e - s, -1, dtype=np.int64)
        lo, hi = max(s, self.R0), min(e, self.R1)
        if lo < hi:
            out[lo - s:hi - s] = self.rmap[lo - self.R0:hi - self.R0]
        return out

    def to_target(self, q0, q1):
        """Oriented (possibly virtual) interval [q0, q1) -> list of forward
        target intervals, in oriented order."""
        L = self.tgt.length
        pieces = []
        if q0 < L < q1:
            segs = [(q0, L), (0, q1 - L)]
        elif q0 >= L:
            segs = [(q0 - L, q1 - L)]
        else:
            segs = [(q0, q1)]
        for a, b in segs:
            pieces.append((a, b) if self.orient == 1 else (L - b, L - a))
        return pieces

    def tpos(self, q):
        L = self.tgt.length
        q = q % L
        return q if self.orient == 1 else L - 1 - q

    def segments(self):
        """(ref_start, ref_end, tgt_start, tgt_end) forward-coordinate pieces
        of the block, split at sequence origins (for drawing / tables)."""
        Lr, Lt = self.ref.length, self.tgt.length
        m = np.nonzero(self.rmap >= 0)[0]
        if m.size == 0:
            return []
        p = (m + self.R0)
        q = self.rmap[m]
        # split where reference or target wraps
        brk = np.nonzero((np.diff(p // Lr) != 0) | (np.diff(q // Lt) != 0))[0] + 1
        out = []
        for seg in np.split(np.arange(m.size), brk):
            if seg.size == 0:
                continue
            ps, pe = int(p[seg[0]] % Lr), int(p[seg[-1]] % Lr) + 1
            qs, qe = int(q[seg[0]] % Lt), int(q[seg[-1]] % Lt) + 1
            if self.orient == 1:
                ts, te = qs, qe
            else:
                ts, te = Lt - qe, Lt - qs
            out.append((ps, pe, ts, te))
        return out


def _stitch_circular(chains, refs, tgts, k, max_gap):
    """Merge blocks that continue across the origin of a circular sequence.
    chains: list of dicts with keys ri, ti, orient, score, I, J, L.
    Merges that keep reference coordinates contiguous (shifting only the
    target into 'virtual' coordinates >= its length) are preferred."""

    def try_merge(a, b, allow_ref_shift):
        if a is b or (a["ri"], a["ti"], a["orient"]) != (b["ri"], b["ti"], b["orient"]):
            return False
        r, t = refs[a["ri"]], tgts[a["ti"]]
        a_ie = int(a["I"][-1] + a["L"][-1])
        a_je = int(a["J"][-1] + a["L"][-1])
        b_i0, b_j0 = int(b["I"][0]), int(b["J"][0])
        if (b["I"][-1] + b["L"][-1]) > r.length or (b["J"][-1] + b["L"][-1]) > t.length:
            return False                                    # b is already virtual
        near_r = min(max_gap, r.length // 4)
        near_t = min(max_gap, t.length // 4)
        sr = r.length if (r.circular and r.length - near_r < a_ie <= r.length + k
                          and b_i0 < near_r) else 0
        st = t.length if (t.circular and t.length - near_t < a_je <= t.length + k
                          and b_j0 < near_t) else 0
        if (not sr and not st) or (sr and not allow_ref_shift):
            return False
        bi0, bj0 = b_i0 + sr, b_j0 + st
        if not (bi0 >= a_ie - k and bj0 >= a_je - k and bi0 - a_ie <= max_gap
                and bj0 - a_je <= max_gap):
            return False
        b_ie = int(b["I"][-1] + b["L"][-1]) + sr
        b_je = int(b["J"][-1] + b["L"][-1]) + st
        if b_ie - int(a["I"][0]) > r.length or b_je - int(a["J"][0]) > t.length:
            return False                                    # would wrap more than once
        a["I"] = np.concatenate([a["I"], b["I"] + sr])
        a["J"] = np.concatenate([a["J"], b["J"] + st])
        a["L"] = np.concatenate([a["L"], b["L"]])
        a["score"] += b["score"]
        return True

    changed = True
    while changed:
        changed = False
        for allow_ref_shift in (False, True):
            for a in sorted(chains, key=lambda c: -c["score"]):
                for b in list(chains):
                    if try_merge(a, b, allow_ref_shift):
                        chains.remove(b)
                        changed = True
                        break
                if changed:
                    break
            if changed:
                break
    return chains


def build_blocks(refs, tgts, k, max_gap, min_score):
    """Find orthologous syntenic blocks between reference and target."""
    t0 = time.time()
    (G, I, J, L), n_anchor = global_anchor_runs(refs, tgts, k)
    LOG.info("  %d unique shared %d-mers (%d runs)  [%.1fs]", n_anchor, k, G.size, time.time() - t0)
    nt2 = 2 * len(tgts)
    raw = []
    for g in np.unique(G):
        sel = G == g
        ri, lab = int(g // nt2), int(g % nt2)
        ti, orient = lab // 2, (-1 if lab % 2 else 1)
        for score, ci, cj, cl in chain_runs(I[sel], J[sel], L[sel], k, max_gap, min_score):
            raw.append(dict(ri=ri, ti=ti, orient=orient, score=score, I=ci, J=cj, L=cl))
    raw = _stitch_circular(raw, refs, tgts, k, max_gap)
    raw.sort(key=lambda c: -c["score"])
    LOG.info("  %d syntenic chains  [%.1fs]", len(raw), time.time() - t0)

    ref_ext = [r.seq + r.seq if r.circular else r.seq for r in refs]
    ref_bytes = [np.frombuffer(s.encode(), dtype=np.uint8) for s in ref_ext]
    tgt_ext, tgt_bytes = {}, {}
    for ti, t in enumerate(tgts):
        for orient in (1, -1):
            s = t.seq if orient == 1 else revcomp(t.seq)
            tgt_ext[(ti, orient)] = s + s if t.circular else s
            tgt_bytes[(ti, orient)] = np.frombuffer(tgt_ext[(ti, orient)].encode(), dtype=np.uint8)

    chains = []
    for c in raw:
        rs, ts = ref_ext[c["ri"]], tgt_ext[(c["ti"], c["orient"])]
        blocks = runs_to_blocks(c["I"], c["J"], c["L"])
        if not blocks:
            continue
        R, T = [], []
        i0, j0, _ = blocks[0]
        ie_last = blocks[-1][0] + blocks[-1][2]
        je_last = blocks[-1][1] + blocks[-1][2]
        rL, tL = refs[c["ri"]].length, tgts[c["ti"]].length
        wrap = (refs[c["ri"]].circular and tgts[c["ti"]].circular
                and 0 <= i0 + rL - ie_last <= max_gap and 0 <= j0 + tL - je_last <= max_gap
                and ie_last <= 2 * rL and je_last <= 2 * tL)
        lr, lt = min(max_gap, i0), min(max_gap, j0)
        lr, lt = min(lr, int(1.3 * lt) + 200), min(lt, int(1.3 * lr) + 200)
        if lr and lt and not wrap:
            a, b = flank_alignment(rs[i0 - lr:i0][::-1], ts[j0 - lt:j0][::-1], k)
            R.append(i0 - 1 - a[::-1])
            T.append(j0 - 1 - b[::-1])
        prev = None
        for (i, j, l) in blocks:
            if prev is not None:
                ai, aj = prev
                if i > ai and j > aj:
                    a, b = align_global(rs[ai:i], ts[aj:j], k)
                    R.append(ai + a)
                    T.append(aj + b)
            R.append(np.arange(i, i + l, dtype=np.int64))
            T.append(np.arange(j, j + l, dtype=np.int64))
            prev = (i + l, j + l)
        ie, je = prev
        if wrap:                       # close the circle: fill up to the block's own start
            a, b = align_global(rs[ie:i0 + rL], ts[je:j0 + tL], k)
            R.append(ie + a)
            T.append(je + b)
        else:
            lr, lt = min(max_gap, len(rs) - ie), min(max_gap, len(ts) - je)
            r_first = int(R[0][0]) if R and R[0].size else i0
            t_first = int(T[0][0]) if T and T[0].size else j0
            if refs[c["ri"]].circular:
                lr = min(lr, r_first + rL - ie)
            if tgts[c["ti"]].circular:
                lt = min(lt, t_first + tL - je)
            lr, lt = min(lr, int(1.3 * lt) + 200), min(lt, int(1.3 * lr) + 200)
            if lr > 0 and lt > 0:
                a, b = flank_alignment(rs[ie:ie + lr], ts[je:je + lt], k)
                R.append(ie + a)
                T.append(je + b)
        R = np.concatenate(R)
        T = np.concatenate(T)
        good = np.ones(R.size, dtype=bool)        # safety: keep strictly collinear pairs
        if R.size > 1:
            bad = (np.diff(R) <= 0) | (np.diff(T) <= 0)
            if bad.any():
                LOG.debug("chain %d: %d non-collinear pairs removed", len(chains), int(bad.sum()))
                keep = [0]
                for x in range(1, R.size):
                    if R[x] > R[keep[-1]] and T[x] > T[keep[-1]]:
                        keep.append(x)
                good[:] = False
                good[keep] = True
        R, T = R[good], T[good]
        ch = Chain(len(chains) + 1, refs[c["ri"]], tgts[c["ti"]], c["orient"], R, T, c["score"])
        ch.ref_idx, ch.tgt_idx = c["ri"], c["ti"]
        rb, tb = ref_bytes[c["ri"]], tgt_bytes[(c["ti"], c["orient"])]
        ch.identity = float((rb[R] == tb[T]).mean()) if R.size else 0.0
        ch.ref_bytes, ch.tgt_bytes = rb, tb
        chains.append(ch)
    LOG.info("  base-level alignment of %d blocks done  [%.1fs]", len(chains), time.time() - t0)
    return chains


# ----------------------------------------------------------------------------
# Lift-over of features and comparison with the reference
# ----------------------------------------------------------------------------
@dataclass
class FeatureResult:
    feat: RefFeature
    status: str = "failed"                 # transferred | failed
    reason: str = ""
    chain: Optional[Chain] = None
    tgt_id: str = ""
    tgt_parts: List[Tuple[int, int, int]] = field(default_factory=list)
    ref_len: int = 0
    tgt_len: int = 0
    aligned: int = 0
    matches: int = 0
    mismatches: int = 0
    ins_bases: int = 0
    ins_events: int = 0
    del_bases: int = 0
    del_events: int = 0
    fs_indels: int = 0                     # indels whose length is not a multiple of 3
    miss5: int = 0
    miss3: int = 0
    parts_lost: int = 0
    identity: float = 0.0
    coverage: float = 0.0
    partial5: bool = False
    partial3: bool = False
    # CDS-specific
    ref_prot: str = ""
    tgt_prot: str = ""
    codon_start: int = 1
    aa_identity: Optional[float] = None
    aa_changes: List[str] = field(default_factory=list)
    syn: int = 0
    nonsyn: int = 0
    cds_flags: List[str] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)

    @property
    def tgt_start(self):
        return min(p[0] for p in self.tgt_parts) if self.tgt_parts else -1

    @property
    def tgt_end(self):
        return max(p[1] for p in self.tgt_parts) if self.tgt_parts else -1

    @property
    def tgt_strand(self):
        s = {p[2] for p in self.tgt_parts}
        return s.pop() if len(s) == 1 else 0

    @property
    def is_identical(self):
        return (self.status == "transferred" and self.matches == self.ref_len == self.tgt_len
                and self.ins_bases == 0 and self.del_bases == 0)

    @property
    def disrupted(self):
        bad = ("frameshift", "premature_stop", "start_codon_lost", "stop_codon_lost",
               "5'_truncated", "3'_truncated")
        return any(fl.split(":")[0] in bad for fl in self.cds_flags)


@dataclass
class LocusResult:
    locus: Locus
    status: str = "failed"
    reason: str = ""
    chain: Optional[Chain] = None
    results: List[FeatureResult] = field(default_factory=list)
    rep: Optional[FeatureResult] = None      # representative feature (CDS/RNA/primary)
    span: Optional[FeatureResult] = None     # gene-span feature (gene or primary)


_PROT_ALIGNER = None


def _prot_aligner():
    global _PROT_ALIGNER
    if _PROT_ALIGNER is None:
        a = PairwiseAligner()
        a.mode = "global"
        a.substitution_matrix = substitution_matrices.load("BLOSUM62")
        a.open_gap_score = -11
        a.extend_gap_score = -1
        _PROT_ALIGNER = a
    return _PROT_ALIGNER


_AA_OK = set("ARNDCQEGHILKMFPSTWYVBZX*")


def compare_proteins(rp, tp):
    """Return (aa identity relative to the reference length, list of changes)
    using reference numbering, e.g. D614G, Q57*, H69_V70del, K12_L13insAG."""
    if not rp:
        return None, []
    if rp == tp:
        return 1.0, []
    if len(rp) == len(tp):
        diffs = [i for i in range(len(rp)) if rp[i] != tp[i]]
        if len(diffs) <= 0.3 * len(rp):
            return 1 - len(diffs) / len(rp), ["%s%d%s" % (rp[i], i + 1, tp[i]) for i in diffs]
    r2 = "".join(c if c in _AA_OK else "X" for c in rp)
    t2 = "".join(c if c in _AA_OK else "X" for c in tp)
    if not t2:
        return 0.0, ["%s1_%s%ddel" % (rp[0], rp[-1], len(rp))]
    if len(r2) * len(t2) > 60_000_000:
        n = min(len(rp), len(tp))
        same = sum(1 for i in range(n) if rp[i] == tp[i])
        return same / len(rp), ["(proteins too long to align; %d positional differences)" % (len(rp) - same)]
    aln = _prot_aligner().align(r2, t2)[0]
    changes, ident = [], 0
    pr = pt = 0
    for a0, a1, b0, b1 in _aligned_blocks(aln):
        if a0 > pr:                                   # deletion in target
            changes.append(("%s%ddel" % (rp[pr], pr + 1)) if a0 - pr == 1 else
                           "%s%d_%s%ddel" % (rp[pr], pr + 1, rp[a0 - 1], a0))
        if b0 > pt:                                   # insertion in target
            left = "%s%d" % (rp[a0 - 1], a0) if a0 > 0 else "0"
            right = "%s%d" % (rp[a0], a0 + 1) if a0 < len(rp) else "end"
            changes.append("%s_%sins%s" % (left, right, tp[pt:b0]))
        for x in range(a1 - a0):
            if rp[a0 + x] == tp[b0 + x]:
                ident += 1
            else:
                changes.append("%s%d%s" % (rp[a0 + x], a0 + x + 1, tp[b0 + x]))
        pr, pt = a1, b1
    if pr < len(rp):
        changes.append(("%s%ddel" % (rp[pr], pr + 1)) if len(rp) - pr == 1 else
                       "%s%d_%s%ddel" % (rp[pr], pr + 1, rp[-1], len(rp)))
    if pt < len(tp):
        changes.append("%s%d_endins%s" % (rp[-1], len(rp), tp[pt:]))
    return ident / len(rp), changes


def translate(nt, table_id):
    nt = nt[: len(nt) - len(nt) % 3]
    if not nt:
        return ""
    return str(Seq(nt).translate(table=table_id))


def _trim_terminal_jumps(q):
    """Unmap a few terminal bases that are separated from the rest of the
    feature by a long insertion in the target (ambiguous gap placement at a
    feature boundary)."""
    idx = np.nonzero(q >= 0)[0]
    if idx.size < 10:
        return q
    big = max(20, int(0.1 * q.size))
    lim = max(3, min(30, int(0.1 * idx.size)))
    steps = np.diff(q[idx])
    q = q.copy()
    head = np.nonzero(steps[:lim] > big)[0]
    if head.size:
        q[idx[:head[-1] + 1]] = -1
    tail = np.nonzero(steps[-lim:] > big)[0]
    if tail.size:
        cut = steps.size - lim + tail[0] + 1
        q[idx[cut:]] = -1
    return q


def _merge_adjacent(pieces):
    """Merge touching forward-coordinate intervals (given in oriented order)."""
    out = []
    for a, b in pieces:
        if out and (out[-1][1] == a or out[-1][0] == b):
            out[-1] = (min(out[-1][0], a), max(out[-1][1], b))
        else:
            out.append((a, b))
    return out


class Lifter:
    def __init__(self, refs, tgts, chains, args):
        self.refs = {r.sid: r for r in refs}
        self.tgts = {t.sid: t for t in tgts}
        self.args = args
        self.chains_by_ref = defaultdict(list)
        for c in chains:
            self.chains_by_ref[c.ref.sid].append(c)
        self.effects = defaultdict(list)          # (ref_id, pos) -> ["gene:D614G", ...]
        self.cds_index = []                        # (ref_id, parts, name, codon_start) for variant effects

    # -- mapping helpers ---------------------------------------------------
    def _part_q(self, chain, s, e):
        """Oriented target coordinate for each reference position in [s, e)
        and the (virtual) reference coordinate it was looked up at."""
        q = chain.lookup(s, e)
        pos = np.arange(s, e, dtype=np.int64)
        ref = self.refs[chain.ref.sid]
        if ref.circular:
            L = ref.length
            q2 = chain.lookup(s + L, e + L)
            use2 = (q < 0) & (q2 >= 0)
            if use2.any():
                qc = np.where(use2, q2, q)
                max_step = max(1000, 3 * (e - s))
                tries = [qc]
                if chain.tgt.circular:          # undo a jump across the target origin
                    Lt = chain.tgt.length
                    tries.append(np.where(use2 & (q2 >= Lt), q2 - Lt, qc))
                    tries.append(np.where(~use2 & (q >= Lt), q - Lt, qc))
                for qt in tries:
                    m = qt[qt >= 0]
                    steps = np.diff(m)
                    if m.size and m.min() >= 0 and (steps.size == 0 or
                                                    (steps.min() > 0 and steps.max() <= max_step)):
                        return qt, np.where(use2, pos + L, pos)
                if (q2 >= 0).sum() > (q >= 0).sum():
                    return q2, pos + L
        return q, pos

    def _coverage(self, chain, feats):
        """Number of reference bases of the locus with an identical aligned
        target base in this block (used to choose the orthologous block),
        and the number of aligned bases."""
        n = n_al = 0
        seen = set()
        for f in feats:
            for s, e, _ in f.parts:
                if (s, e) in seen:
                    continue
                seen.add((s, e))
                if s == e:                       # zero-length site: use flanking bases
                    s, e = max(0, s - 1), s + 1
                q, pos = self._part_q(chain, s, e)
                m = q >= 0
                if m.any():
                    n += int((chain.ref_bytes[pos[m]] == chain.tgt_bytes[q[m]]).sum())
                    n_al += int(m.sum())
        return n, n_al

    # -- per locus -----------------------------------------------------------
    def lift_locus(self, locus):
        args = self.args
        lr = LocusResult(locus)
        cands = [c for c in self.chains_by_ref.get(locus.seqid, [])]
        best, best_key = None, (0, 0, 0)
        core_feats = locus.core_features()
        for c in cands:
            m_core, _ = self._coverage(c, core_feats)
            m_all, al_all = self._coverage(c, locus.features)
            key = (m_core, m_all, al_all)
            if key > best_key or (key == best_key and best is not None and c.score > best.score):
                best, best_key = c, key
        best_cov = best_key[2]
        lr.chain = best
        for f in locus.features:
            fr = self.lift_feature(f, best) if best is not None and best_cov > 0 else \
                FeatureResult(f, reason="no orthologous sequence found (absent or too divergent)",
                              ref_len=f.length)
            lr.results.append(fr)
        core = set(id(f) for f in locus.core_features())
        core_res = [r for r in lr.results if id(r.feat) in core]
        ok_core = [r for r in core_res if r.status == "transferred"]
        if ok_core:
            lr.status = "transferred"
            # keep mapped parent features (gene/mRNA) so the hierarchy stays intact
            for r in lr.results:
                if r.status != "transferred" and r.tgt_parts and \
                        (r.feat.ftype in GENE_PARENT_TYPES or r.feat.ftype == "mRNA"):
                    r.status = "transferred"
                    r.notes.append("kept as parent of transferred features (coverage %.1f%%, "
                                   "identity %.1f%%)" % (100 * r.coverage, 100 * r.identity))
        else:
            lr.status = "failed"
            best_r = max(core_res, key=lambda r: (r.coverage, r.identity))
            lr.reason = best_r.reason or "below thresholds"
            for r in lr.results:
                if r.status == "transferred":
                    r.status = "failed"
                    r.reason = "gene product not transferred"
        pool = ok_core or core_res
        lr.rep = max(pool, key=lambda r: (r.feat.ftype == "CDS", r.feat.length))
        spans = [r for r in lr.results if r.feat is locus.primary]
        lr.span = spans[0] if spans else lr.rep
        if lr.status == "transferred" and lr.span.status != "transferred":
            lr.span = lr.rep
        return lr

    # -- per feature -----------------------------------------------------------
    def _terminal_deleted(self, f, chain, q_bio, five):
        """True if the reference base just outside the feature end aligns
        directly next to the feature's last aligned target base."""
        mapped = np.nonzero(q_bio >= 0)[0]
        q_edge = int(q_bio[mapped[0]] if five else q_bio[mapped[-1]])
        s, e, st = f.parts[0] if five else f.parts[-1]
        fwd = (st or 1) > 0
        if five:
            beyond = s - 1 if fwd else e
        else:
            beyond = e if fwd else s - 1
        L = self.refs[f.seqid].length
        circ = self.refs[f.seqid].circular
        direction = 1 if beyond == e else -1      # walk away from the feature
        qb = -1
        for k in range(60):                     # first aligned base outside the feature
            p = beyond + direction * k
            if p < 0 or p >= L:
                if not circ:
                    return False
                p %= L
            qq, _ = self._part_q(chain, p, p + 1)
            if qq[0] >= 0:
                qb = int(qq[0])
                break
        if qb < 0:
            return False
        step = qb - q_edge
        if chain.tgt.circular:
            Lt = chain.tgt.length
            step = (step + Lt) % Lt
            step = step - Lt if step > Lt // 2 else step
        # in reference-forward order the outside base precedes (5' of + strand) or follows
        expect = -1 if (five == fwd) else 1
        return step == expect

    def _lift_point(self, f, chain, fr):
        """Zero-length feature (a site between two bases, e.g. 5933^5934)."""
        s, _, st = f.parts[0]
        if s == 0:
            fr.reason = "site at sequence start"
            return fr
        q, _ = self._part_q(chain, s - 1, s + 1)
        if (q >= 0).all() and q[1] - q[0] == 1:
            t0, t1 = chain.tpos(int(q[0])), chain.tpos(int(q[1]))
            p = max(t0, t1)
            fr.tgt_parts = [(p, p, st * chain.orient if st else 0)]
            fr.status, fr.identity, fr.coverage = "transferred", 1.0, 1.0
        else:
            fr.reason = "site not conserved (flanking bases not aligned)"
        return fr

    def lift_feature(self, f, chain):
        args = self.args
        fr = FeatureResult(f, chain=chain, tgt_id=chain.tgt.sid, ref_len=f.length)
        if f.length == 0:
            return self._lift_point(f, chain, fr)
        rb, tb = chain.ref_bytes, chain.tgt_bytes
        bio_q, bio_rpos, bio_part, bio_strand = [], [], [], []
        pieces_all = []
        for pi, (s, e, st) in enumerate(f.parts):
            q, pos = self._part_q(chain, s, e)
            q = _trim_terminal_jumps(q)
            ss = st if st else 1
            if ss < 0:
                bio_q.append(q[::-1])
                bio_rpos.append(pos[::-1])
            else:
                bio_q.append(q)
                bio_rpos.append(pos)
            bio_part.append(np.full(e - s, pi, dtype=np.int64))
            bio_strand.append(np.full(e - s, ss, dtype=np.int64))
            m = q >= 0
            if not m.any():
                fr.parts_lost += 1
                continue
            idx = np.nonzero(m)[0]
            f0, f1 = int(idx[0]), int(idx[-1])
            fr.aligned += idx.size
            fr.matches += int((rb[pos[idx]] == tb[q[idx]]).sum())
            mm = m[f0:f1 + 1]
            fr.del_bases += int((~mm).sum())
            # deletion events and their lengths
            if (~mm).any():
                d = np.diff(np.concatenate(([0], (~mm).astype(np.int8), [0])))
                starts, ends = np.nonzero(d == 1)[0], np.nonzero(d == -1)[0]
                fr.del_events += starts.size
                fr.fs_indels += int(((ends - starts) % 3 != 0).sum())
            qi = q[idx]
            steps = np.diff(qi)
            jumps = steps <= 0                     # seam of a block stitched across an origin
            gaps = np.where(jumps, 0, steps - 1)
            fr.ins_bases += int(gaps.sum())
            fr.ins_events += int((gaps > 0).sum())
            fr.fs_indels += int(((gaps > 0) & (gaps % 3 != 0)).sum())
            tst = st * chain.orient if st else 0
            pieces = []
            for run in np.split(qi, np.nonzero(jumps)[0] + 1):
                pieces += chain.to_target(int(run[0]), int(run[-1]) + 1)
            pieces = _merge_adjacent(pieces)
            if ss < 0:
                pieces = pieces[::-1]
            pieces_all += [(a, b, tst) for a, b in pieces]
        fr.tgt_parts = pieces_all
        fr.tgt_len = sum(b - a for a, b, _ in pieces_all)
        fr.mismatches = fr.aligned - fr.matches
        q_bio = np.concatenate(bio_q)
        mapped_bio = q_bio >= 0
        if mapped_bio.any():
            nz = np.nonzero(mapped_bio)[0]
            fr.miss5 = int(nz[0])
            fr.miss3 = int(q_bio.size - 1 - nz[-1])
        else:
            fr.miss5 = fr.miss3 = int(q_bio.size)
        # bases missing at a feature end while the target continues in register
        # are a terminal deletion, not a truncation of the feature
        if 0 < fr.miss5 < q_bio.size and self._terminal_deleted(f, chain, q_bio, five=True):
            fr.del_bases += fr.miss5
            fr.del_events += 1
            fr.fs_indels += int(fr.miss5 % 3 != 0)
            fr.notes.append("5' terminal deletion of %d nt" % fr.miss5)
            fr.miss5 = 0
        if 0 < fr.miss3 < q_bio.size and self._terminal_deleted(f, chain, q_bio, five=False):
            fr.del_bases += fr.miss3
            fr.del_events += 1
            fr.fs_indels += int(fr.miss3 % 3 != 0)
            fr.notes.append("3' terminal deletion of %d nt" % fr.miss3)
            fr.miss3 = 0
        fr.coverage = fr.aligned / fr.ref_len if fr.ref_len else 0.0
        cols = fr.aligned + fr.del_bases + fr.ins_bases
        fr.identity = fr.matches / cols if cols else 0.0
        fr.partial5 = f.partial5 or fr.miss5 > 0
        fr.partial3 = f.partial3 or fr.miss3 > 0
        if fr.parts_lost:
            fr.notes.append("%d of %d parts not aligned" % (fr.parts_lost, len(f.parts)))

        if not pieces_all:
            fr.reason = "not aligned within the orthologous block"
        elif fr.coverage < args.min_coverage:
            fr.reason = "low coverage (%.1f%% < %.0f%%)" % (100 * fr.coverage, 100 * args.min_coverage)
        elif fr.identity < args.min_identity and fr.ref_len >= 20:   # tiny sites: position only
            fr.reason = "low identity (%.1f%% < %.0f%%)" % (100 * fr.identity, 100 * args.min_identity)
        else:
            fr.status = "transferred"
        if fr.status == "transferred" and f.ftype == "CDS":
            self.analyse_cds(fr, chain, q_bio, np.concatenate(bio_rpos),
                             np.concatenate(bio_part), np.concatenate(bio_strand))
        return fr

    # -- CDS checks ------------------------------------------------------------
    def analyse_cds(self, fr, chain, q_bio, rpos_bio, part_bio, strand_bio):
        f = fr.feat
        try:
            table_id = int(f.q("transl_table", self.args.transl_table))
            table = CodonTable.unambiguous_dna_by_id[table_id]
        except (ValueError, KeyError):
            table_id, table = 1, CodonTable.unambiguous_dna_by_id[1]
        try:
            cs = int(f.q("codon_start", 1))
        except ValueError:
            cs = 1
        off = cs - 1
        ref_seq = self.refs[f.seqid].seq
        tgt_seq = self.tgts[fr.tgt_id].seq
        ref_cds = extract_parts(ref_seq, f.parts)
        tgt_cds = extract_parts(tgt_seq, fr.tgt_parts)
        t_off = (off - fr.miss5) % 3
        fr.codon_start = t_off + 1
        rp = translate(ref_cds[off:], table_id)
        tp = translate(tgt_cds[t_off:], table_id)
        fr.ref_prot, fr.tgt_prot = rp, tp
        starts, stops = set(table.start_codons), set(table.stop_codons)
        flags = []
        ref_start_ok = (not f.partial5) and off == 0 and ref_cds[:3] in starts
        ref_stop_ok = (not f.partial3) and (len(ref_cds) - off) % 3 == 0 and ref_cds[-3:] in stops
        tgt_start_ok = fr.miss5 == 0 and t_off == 0 and tgt_cds[:3] in starts
        tgt_stop_ok = fr.miss3 == 0 and (len(tgt_cds) - t_off) % 3 == 0 and tgt_cds[-3:] in stops
        if fr.miss5 > 0 and not f.partial5:
            flags.append("5'_truncated:%dnt" % fr.miss5)
        elif ref_start_ok and not tgt_start_ok:
            flags.append("start_codon_lost:%s>%s" % (ref_cds[:3], tgt_cds[:3]))
        if fr.miss3 > 0 and not f.partial3:
            flags.append("3'_truncated:%dnt" % fr.miss3)
        elif ref_stop_ok and not tgt_stop_ok:
            flags.append("stop_codon_lost")
        if (len(tgt_cds) - t_off) % 3 != (len(ref_cds) - off) % 3 or fr.fs_indels:
            flags.append("frameshift:%d_indel(s)" % fr.fs_indels if fr.fs_indels else "frameshift")
        r_core = rp[:-1] if (ref_stop_ok and rp.endswith("*")) else rp
        t_core = tp[:-1] if (tgt_stop_ok and tp.endswith("*")) else tp
        if t_core.count("*") > r_core.count("*"):
            r_stops = {i for i, c in enumerate(r_core) if c == "*"}
            first = next(i for i, c in enumerate(t_core) if c == "*" and i not in r_stops)
            flags.append("premature_stop:aa%d" % (first + 1))
        fr.cds_flags = flags
        fr.aa_identity, fr.aa_changes = compare_proteins(r_core, t_core)

        # codon-level synonymous / non-synonymous changes from the alignment
        rb, tb = chain.ref_bytes, chain.tgt_bytes
        n = q_bio.size
        ncod = (n - off) // 3
        if ncod <= 0:
            return
        comp = np.arange(256, dtype=np.uint8)
        for a, b in zip(b"ACGTN", b"TGCAN"):
            comp[a] = b
        rbase = rb[rpos_bio]
        rbase = np.where(strand_bio < 0, comp[rbase], rbase)
        qq = np.where(q_bio >= 0, q_bio, 0)
        tbase = tb[qq]
        tbase = np.where(strand_bio < 0, comp[tbase], tbase)
        ci = off + 3 * np.arange(ncod)
        q3 = np.stack([q_bio[ci], q_bio[ci + 1], q_bio[ci + 2]], 1)
        p3 = np.stack([part_bio[ci], part_bio[ci + 1], part_bio[ci + 2]], 1)
        ok = (q3 >= 0).all(1)
        ok &= (p3[:, 0] != p3[:, 1]) | (np.abs(q3[:, 1] - q3[:, 0]) == 1)
        ok &= (p3[:, 1] != p3[:, 2]) | (np.abs(q3[:, 2] - q3[:, 1]) == 1)
        r3 = np.stack([rbase[ci], rbase[ci + 1], rbase[ci + 2]], 1)
        t3 = np.stack([tbase[ci], tbase[ci + 1], tbase[ci + 2]], 1)
        diff = ok & (r3 != t3).any(1)
        fwd = table.forward_table
        gname = f.name
        L = self.refs[f.seqid].length
        for c in np.nonzero(diff)[0]:
            rc = bytes(r3[c]).decode()
            tc = bytes(t3[c]).decode()
            ra = "*" if rc in stops else fwd.get(rc, "X")
            ta = "*" if tc in stops else fwd.get(tc, "X")
            if ra == ta:
                fr.syn += 1
                eff = "%s:synonymous(%s%d%s)" % (gname, ra, c + 1, ta)
            else:
                fr.nonsyn += 1
                eff = "%s:%s%d%s" % (gname, ra, c + 1, ta)
            for x in range(3):
                if r3[c, x] != t3[c, x]:
                    self.effects[(f.seqid, int(rpos_bio[ci[c] + x]) % L)].append(eff)


# ----------------------------------------------------------------------------
# Genome-wide nucleotide differences
# ----------------------------------------------------------------------------
def call_variants(chains, refs):
    """Nucleotide differences reference -> target.  Each reference position
    is taken from the best-scoring block that aligns it; positions deleted in
    the target are taken from the best block spanning them.
    Types: SNP, INS, DEL, COMPLEX (small unaligned replacement) and
    UNALIGNED (larger region without detectable homology)."""
    out = []
    by_ref = defaultdict(list)
    for c in chains:
        by_ref[c.ref.sid].append(c)
    for r in refs:
        L = r.length
        cs = sorted(by_ref.get(r.sid, []), key=lambda c: -c.score)
        owner = np.full(L, -1, dtype=np.int64)
        owner_v = np.full(L, -1, dtype=np.int64)
        info = []
        for c in cs:                                  # pass 1: aligned positions
            pos = np.arange(c.R0, c.R1) % L
            mapped = c.rmap >= 0
            idx = np.nonzero(mapped)[0]
            inside = np.zeros(c.rmap.size, dtype=bool)
            if idx.size:
                inside[idx[0]:idx[-1] + 1] = True
            info.append((c, pos, mapped, inside))
            free = np.nonzero(mapped & (owner[pos] == -1))[0]
            _, first = np.unique(pos[free], return_index=True)
            free = free[first]
            owner[pos[free]] = c.cid
            owner_v[pos[free]] = free
        for c, pos, mapped, inside in info:           # pass 2: deleted positions
            free = np.nonzero(inside & ~mapped & (owner[pos] == -1))[0]
            _, first = np.unique(pos[free], return_index=True)
            free = free[first]
            owner[pos[free]] = c.cid
            owner_v[pos[free]] = free
        for c, pos, mapped, inside in info:
            vidx = np.arange(c.rmap.size)
            own = inside & (owner[pos] == c.cid) & (owner_v[pos] == vidx)
            rb, tb = c.ref_bytes, c.tgt_bytes
            pv = np.arange(c.R0, c.R1)
            qv = np.where(mapped, c.rmap, 0)
            for x in np.nonzero(own & mapped & (rb[pv] != tb[qv]))[0]:
                out.append(dict(ref=r.sid, pos=int(pos[x]) + 1, ref_allele=chr(rb[pv[x]]),
                                alt=chr(tb[c.rmap[x]]), type="SNP", tgt=c.tgt.sid,
                                tpos=c.tpos(int(c.rmap[x])) + 1, chain=c.cid, length=1))
            gap_after_del = set()
            dl = inside & ~mapped
            if dl.any():
                d = np.diff(np.concatenate(([0], dl.astype(np.int8), [0])))
                for a, b in zip(np.nonzero(d == 1)[0], np.nonzero(d == -1)[0]):
                    if own[a:b].mean() < 0.5:
                        continue                         # mostly explained by another block
                    q_prev = int(c.rmap[a - 1])
                    n_t = int(c.ins_after[a - 1])        # target bases opposite the gap
                    gap_after_del.add(a - 1)
                    ref_al = bytes(rb[pv[a]:pv[b - 1] + 1]).decode()
                    alt = bytes(tb[q_prev + 1:q_prev + 1 + n_t]).decode() if n_t else "-"
                    if n_t == 0:
                        vtype = "DEL"
                    elif max(b - a, n_t) <= 50:
                        vtype = "COMPLEX"
                    else:
                        vtype = "UNALIGNED"
                        ref_al = "%d bp" % (b - a)
                        alt = "%d bp" % n_t
                    out.append(dict(ref=r.sid, pos=int(pos[a]) + 1, ref_allele=ref_al, alt=alt,
                                    type=vtype, tgt=c.tgt.sid, tpos=c.tpos(q_prev + 1) + 1 if n_t
                                    else c.tpos(q_prev) + 1, chain=c.cid, length=int(b - a),
                                    tlen=n_t))
            for x in np.nonzero(own & mapped & (c.ins_after > 0))[0]:
                if x in gap_after_del:
                    continue
                q0 = int(c.rmap[x])
                n = int(c.ins_after[x])
                out.append(dict(ref=r.sid, pos=int(pos[x]) + 1, ref_allele="-",
                                alt=bytes(tb[q0 + 1:q0 + 1 + n]).decode() if n <= 1000 else "%d bp" % n,
                                type="INS", tgt=c.tgt.sid, tpos=c.tpos(q0 + 1) + 1, chain=c.cid,
                                length=n))
    out.sort(key=lambda v: (v["ref"], v["pos"]))
    return out


def annotate_variants(variants, lifter, loci_results):
    """Add overlapping genes and coding effects to variants."""
    if not variants:
        return
    by_ref = defaultdict(list)
    for i, v in enumerate(variants):
        by_ref[v["ref"]].append(i)
    genes = defaultdict(list)
    effects = defaultdict(list)
    for ref_id, idxs in by_ref.items():
        pos0 = np.array([variants[i]["pos"] - 1 for i in idxs])
        order = np.argsort(pos0, kind="stable")
        sp = pos0[order]
        for lr in loci_results:
            loc = lr.locus
            if loc.seqid != ref_id:
                continue
            name = loc.name
            for f in loc.features:
                if f is not loc.primary and f.ftype != "CDS":
                    continue
                for s, e, _ in f.parts:
                    a, b = np.searchsorted(sp, s, "left"), np.searchsorted(sp, e, "left")
                    for o in order[a:b]:
                        vi = idxs[o]
                        if f is loc.primary and name not in genes[vi]:
                            genes[vi].append(name)
                        if f.ftype == "CDS" and lr.status == "transferred":
                            v = variants[vi]
                            if v["type"] == "SNP":
                                eff = lifter.effects.get((ref_id, v["pos"] - 1))
                                for e_ in (eff or ["%s:complex" % name]):
                                    if e_ not in effects[vi]:
                                        effects[vi].append(e_)
                            else:
                                net = v["length"] if v["type"] != "COMPLEX" else abs(v["length"] - v.get("tlen", 0))
                                if v["type"] == "UNALIGNED":
                                    kind = "unaligned_region"
                                elif net % 3:
                                    kind = "frameshift"
                                elif v["type"] == "COMPLEX":
                                    kind = "inframe_complex"
                                else:
                                    kind = "inframe_deletion" if v["type"] == "DEL" else "inframe_insertion"
                                e_ = "%s:%s" % (name, kind)
                                if e_ not in effects[vi]:
                                    effects[vi].append(e_)
    for i, v in enumerate(variants):
        v["genes"] = ",".join(genes.get(i, [])) or "-"
        if effects.get(i):
            v["effect"] = ";".join(effects[i])
        else:
            v["effect"] = "non-coding" if genes.get(i) else "intergenic"


# ----------------------------------------------------------------------------
# Output: GenBank
# ----------------------------------------------------------------------------
INSDC_KEYS = set(SO_TYPES) | {"misc_signal", "misc_structure", "modified_base", "old_sequence",
                              "oriT", "telomere", "centromere", "iDNA", "C_region", "J_segment",
                              "N_region", "S_region", "V_region", "V_segment", "D_segment"}
SO2GB = {}
for _gb, _so in SO_TYPES.items():
    SO2GB.setdefault(_so, _gb)
SO2GB.update({"ncRNA_gene": "gene", "protein_coding_gene": "gene", "pseudogene": "gene",
              "transcript": "misc_RNA", "lnc_RNA": "ncRNA", "lncRNA": "ncRNA",
              "snRNA": "ncRNA", "snoRNA": "ncRNA", "miRNA": "ncRNA", "ncRNA": "ncRNA",
              "mature_protein_region": "mat_peptide", "signal_peptide": "sig_peptide",
              "pseudogenic_transcript": "misc_RNA", "pseudogenic_tRNA": "tRNA",
              "pseudogenic_rRNA": "rRNA", "primary_transcript": "precursor_RNA"})
GFF_SKIP_ATTRS = {"Is_circular", "Target", "Gap", "Derives_from", "start_range", "end_range",
                  "partial", "codon_start"}


def _make_location(parts, partial5, partial3):
    locs = []
    n = len(parts)
    for i, (s, e, st) in enumerate(parts):
        strand = st if st else None
        sp, ep = ExactPosition(s), ExactPosition(e)
        fwd = (strand or 1) > 0
        if i == 0 and partial5:
            if fwd:
                sp = BeforePosition(s)
            else:
                ep = AfterPosition(e)
        if i == n - 1 and partial3:
            if fwd:
                ep = AfterPosition(e)
            else:
                sp = BeforePosition(s)
        locs.append(SimpleLocation(sp, ep, strand=strand))
    return locs[0] if n == 1 else CompoundLocation(locs, "join")


def _transfer_note(fr, ref_id):
    txt = "annotation transferred from %s (nt identity %.1f%%, coverage %.1f%%)" % (
        fmt_loc(ref_id, fr.feat.parts), 100 * fr.identity, 100 * fr.coverage)
    if fr.feat.ftype == "CDS" and fr.aa_identity is not None:
        txt += "; aa identity %.1f%%" % (100 * fr.aa_identity)
    return txt


def _gb_qualifiers(fr, ref_id):
    f = fr.feat
    q = OrderedDict()
    if f.fmt == "genbank":
        for k, v in f.qualifiers.items():
            if k in ("translation", "protein_id", "codon_start"):
                continue
            q[k] = list(v)
    else:
        has_gene = "gene" in f.qualifiers
        for k, v in f.qualifiers.items():
            if k in GFF_SKIP_ATTRS:
                continue
            if k == "Name":
                if has_gene and v == f.qualifiers.get("gene"):
                    continue
                key = "gene" if (not has_gene and (f.ftype in GENE_LOCUS_TYPES)) else "standard_name"
            else:
                key = GFF2GB_ATTR.get(k, k)
            q.setdefault(key, [])
            q[key].extend(v)
    gb_type = _gb_type(f)
    if f.fmt == "gff3" and "gene" not in q and getattr(f, "locus_name", None) and \
            gb_type not in ("repeat_region", "misc_feature"):
        q["gene"] = [f.locus_name]                 # GenBank convention: children carry /gene
        q.move_to_end("gene", last=False)
    if gb_type == "misc_feature" and f.ftype not in ("misc_feature", "sequence_feature"):
        q.setdefault("note", []).insert(0, "feature type: %s" % f.ftype)
    notes = [_transfer_note(fr, ref_id)]
    if f.ftype == "CDS":
        q["codon_start"] = [str(fr.codon_start)]
        pid = f.q("protein_id")
        if pid:
            notes.append("reference protein_id %s" % pid)
        if fr.cds_flags:
            notes.append("differences from reference: %s" % ", ".join(fr.cds_flags))
        internal_stop = any(fl.startswith("premature_stop") for fl in fr.cds_flags)
        if fr.tgt_prot and not internal_stop and not fr.disrupted:
            q["translation"] = [fr.tgt_prot.rstrip("*")]
        elif internal_stop:
            notes.append("possible pseudogene (internal stop codon)")
    notes += fr.notes
    q.setdefault("note", []).extend(notes)
    return q


def _gb_type(f):
    if f.fmt == "genbank":
        return f.ftype
    if f.ftype in INSDC_KEYS:
        return f.ftype
    return SO2GB.get(f.ftype, "misc_feature")


def write_genbank(path, tgts, loci_results, ref_by_id, args):
    by_tgt = defaultdict(list)
    for lr in loci_results:
        for fr in lr.results:
            if fr.status == "transferred":
                by_tgt[fr.tgt_id].append(fr)
    type_rank = {"gene": 0, "mRNA": 1, "CDS": 3}
    mol = "DNA"
    org = ""
    for r in ref_by_id.values():
        mol = r.molecule_type or mol
        org = org or r.organism
    today = datetime.date.today().strftime("%d-%b-%Y").upper()
    records = []
    for t in tgts:
        rec = SeqRecord(Seq(t.seq), id=t.sid, name=re.sub(r"[^\w.-]", "_", t.sid)[:16],
                        description=t.description or "annotated by transfer_annotation.py")
        rec.annotations["molecule_type"] = mol
        rec.annotations["topology"] = "circular" if t.circular else "linear"
        rec.annotations["data_file_division"] = "UNA"
        rec.annotations["date"] = today
        src_q = OrderedDict()
        if org and args.copy_organism:
            src_q["organism"] = [org]
        src_q["mol_type"] = ["genomic RNA" if "RNA" in mol else "genomic DNA"]
        src_q["note"] = ["annotation transferred from %s with transfer_annotation.py v%s"
                         % (",".join(ref_by_id), __version__)]
        rec.features.append(SeqFeature(SimpleLocation(0, t.length, strand=1), type="source",
                                       qualifiers=src_q))
        frs = sorted(by_tgt.get(t.sid, []), key=lambda fr: (fr.tgt_start, -fr.tgt_end,
                                                              type_rank.get(fr.feat.ftype, 2)))
        for fr in frs:
            loc = _make_location(fr.tgt_parts, fr.partial5, fr.partial3)
            rec.features.append(SeqFeature(loc, type=_gb_type(fr.feat),
                                           qualifiers=_gb_qualifiers(fr, fr.feat.seqid)))
        records.append(rec)
    with open(path, "w") as fh:
        SeqIO.write(records, fh, "genbank")


# ----------------------------------------------------------------------------
# Output: GFF3
# ----------------------------------------------------------------------------
def _gff_esc(s):
    s = str(s)
    for a, b in (("%", "%25"), (";", "%3B"), ("=", "%3D"), ("&", "%26"), (",", "%2C"),
                 ("\t", "%09"), ("\n", "%0A"), ("\r", "%0D")):
        s = s.replace(a, b)
    return s


def _cds_phases(parts, codon_start):
    phases, done = [], 0
    p0 = codon_start - 1
    for i, (s, e, _) in enumerate(parts):
        if i == 0:
            phases.append(p0)
        else:
            phases.append((3 - ((done - p0) % 3)) % 3)
        done += e - s
    return phases


def write_gff3(path, tgts, loci_results, include_fasta=False):
    lines = ["##gff-version 3"]
    for t in tgts:
        lines.append("##sequence-region %s 1 %d" % (t.sid, t.length))
    used_ids = set()

    def uniq(x):
        base, i = x, 1
        while x in used_ids:
            i += 1
            x = "%s-%d" % (base, i)
        used_ids.add(x)
        return x

    rows = []
    tgt_order = {t.sid: i for i, t in enumerate(tgts)}
    for t in tgts:
        if t.circular:
            rows.append((tgt_order[t.sid], 0, -t.length, -1, [
                "%s\ttransfer_annotation\tregion\t1\t%d\t.\t+\t.\tID=%s;Is_circular=true"
                % (t.sid, t.length, _gff_esc(uniq("region-" + t.sid)))]))
    for lr in loci_results:
        if lr.status != "transferred":
            continue
        frs = [fr for fr in lr.results if fr.status == "transferred"]
        if not frs:
            continue
        loc = lr.locus
        base = re.sub(r"[^\w.:-]", "_", loc.lid)
        ids = {}
        # --- IDs and parents ---
        if frs[0].feat.fmt == "gff3":
            transferred_ids = {fr.feat.gff_id for fr in frs if fr.feat.gff_id}
            by_id = {fr.feat.gff_id: fr.feat for fr in lr.results if fr.feat.gff_id}
            for fr in frs:
                ids[id(fr)] = uniq(fr.feat.gff_id) if fr.feat.gff_id else None
            parents = {}
            for fr in frs:
                ps = []
                for p in fr.feat.gff_parents:       # climb to nearest transferred ancestor
                    seen = set()
                    while p and p not in transferred_ids and p in by_id and p not in seen:
                        seen.add(p)
                        nxt = by_id[p].gff_parents
                        p = nxt[0] if nxt else None
                    if p in transferred_ids and p not in ps:
                        ps.append(p)
                parents[id(fr)] = ps
        else:
            genes = [fr for fr in frs if fr.feat.ftype in GENE_PARENT_TYPES]
            rnas = [fr for fr in frs if fr.feat.ftype in RNA_TYPES]
            mrnas = [fr for fr in rnas if fr.feat.ftype == "mRNA"]
            cdss = [fr for fr in frs if fr.feat.ftype == "CDS"]
            for fr in genes:
                ids[id(fr)] = uniq("gene-" + base)
            for fr in rnas:
                ids[id(fr)] = uniq("rna-" + base)
            for fr in cdss:
                ids[id(fr)] = uniq("cds-" + base)
            for fr in frs:
                if id(fr) not in ids:
                    ids[id(fr)] = uniq("%s-%s" % (re.sub(r"[^\w]", "_", fr.feat.ftype), base))
            gene_id = ids[id(genes[0])] if genes else None
            parents = {}
            for fr in frs:
                t = fr.feat.ftype
                if t in GENE_PARENT_TYPES:
                    p = None
                elif t in RNA_TYPES:
                    p = gene_id
                elif t == "CDS":
                    p = ids[id(mrnas[0])] if len(mrnas) == 1 else gene_id
                elif t in ("exon", "intron", "5'UTR", "3'UTR"):
                    p = ids[id(rnas[0])] if len(rnas) == 1 else gene_id
                elif t in PEPTIDE_TYPES:
                    p = ids[id(cdss[0])] if len(cdss) == 1 else gene_id
                else:
                    p = gene_id
                parents[id(fr)] = [p] if p else []
        # --- lines ---
        rank = {"gene": 0, "pseudogene": 0, "mRNA": 1, "exon": 2, "CDS": 3}
        for fr in frs:
            f = fr.feat
            attrs = OrderedDict()                  # key -> list of raw values
            if f.fmt == "gff3":
                ftype = f.ftype
                if ids[id(fr)]:
                    attrs["ID"] = [ids[id(fr)]]
                if parents[id(fr)]:
                    attrs["Parent"] = list(parents[id(fr)])
                for k, v in f.qualifiers.items():
                    if k in ("codon_start", "start_range", "end_range", "partial"):
                        continue
                    attrs[k] = list(v)
            else:
                ftype = SO_TYPES.get(f.ftype, f.ftype)
                attrs["ID"] = [ids[id(fr)]]
                if parents[id(fr)]:
                    attrs["Parent"] = [parents[id(fr)][0]]
                nm = f.q("gene") or f.q("locus_tag") or f.q("product") or f.q("label")
                if nm:
                    attrs["Name"] = [nm]
                for k, v in f.qualifiers.items():
                    if k in ("translation", "codon_start"):
                        continue
                    key = GB2GFF_ATTR.get(k, k)
                    attrs.setdefault(key, []).extend([x for x in v if x != ""] or ["true"])
            attrs["transfer_ref_location"] = [fmt_loc(f.seqid, f.parts)]
            attrs["transfer_identity"] = ["%.4f" % fr.identity]
            attrs["transfer_coverage"] = ["%.4f" % fr.coverage]
            if f.ftype == "CDS":
                if fr.aa_identity is not None:
                    attrs["transfer_aa_identity"] = ["%.4f" % fr.aa_identity]
                if fr.cds_flags:
                    attrs["transfer_flags"] = list(fr.cds_flags)
            if fr.partial5 or fr.partial3:
                attrs["partial"] = ["true"]
            attr_txt = ";".join("%s=%s" % (_gff_esc(k), ",".join(_gff_esc(x) for x in v))
                                for k, v in attrs.items())
            phases = _cds_phases(fr.tgt_parts, fr.codon_start) if f.ftype == "CDS" else \
                ["."] * len(fr.tgt_parts)
            out = []
            for (s, e, st), ph in zip(fr.tgt_parts, phases):
                strand = {1: "+", -1: "-"}.get(st, ".")
                if s == e:                         # zero-length site: start = end = base to the left
                    s, e = max(0, s - 1), max(1, s)
                out.append("\t".join([fr.tgt_id, "transfer_annotation", ftype, str(s + 1), str(e),
                                      f.gff_score if f.fmt == "gff3" else ".", strand, str(ph),
                                      attr_txt]))
            rows.append((tgt_order[fr.tgt_id], lr.span.tgt_start if lr.span.tgt_parts else fr.tgt_start,
                         rank.get(f.ftype, 4), fr.tgt_start, out))
    rows.sort(key=lambda r: r[:4])
    for r in rows:
        lines.extend(r[4])
    if include_fasta:
        lines.append("##FASTA")
        for t in tgts:
            lines.append(">" + t.sid)
            lines.extend(t.seq[i:i + 60] for i in range(0, t.length, 60))
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


# ----------------------------------------------------------------------------
# Output: tables
# ----------------------------------------------------------------------------
def _pct(x):
    return "" if x is None else "%.2f" % (100 * x)


def write_tsv(path, header, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join("" if v is None else str(v).replace("\t", " ") for v in r) + "\n")


GENE_HEADER = ["locus_id", "gene", "locus_tag", "type", "product", "status", "reason",
               "ref_location", "target_location", "orientation", "ref_length", "target_length",
               "length_diff", "nt_identity_pct", "coverage_pct", "mismatches", "insertions",
               "inserted_bases", "deletions", "deleted_bases", "missing_5prime", "missing_3prime",
               "aa_identity_pct", "synonymous", "nonsynonymous", "aa_changes", "cds_flags",
               "features_transferred", "block"]


def gene_row(lr):
    loc, rep, span = lr.locus, lr.rep, lr.span
    ok = lr.status == "transferred"
    tloc = fmt_loc(span.tgt_id, span.tgt_parts) if ok and span.tgt_parts else "-"
    orient = "-"
    if ok and lr.chain is not None:
        orient = "same" if lr.chain.orient == 1 else "inverted"
    n_ok = sum(1 for r in lr.results if r.status == "transferred")
    return [loc.lid, loc.name, loc.locus_tag, loc.kind, loc.product, lr.status, lr.reason,
            fmt_loc(loc.seqid, span.feat.parts), tloc, orient, rep.ref_len,
            rep.tgt_len if rep.tgt_parts else "", (rep.tgt_len - rep.ref_len) if rep.tgt_parts else "",
            _pct(rep.identity) if rep.tgt_parts else "", _pct(rep.coverage) if rep.tgt_parts else "",
            rep.mismatches, rep.ins_events, rep.ins_bases, rep.del_events, rep.del_bases,
            rep.miss5, rep.miss3, _pct(rep.aa_identity), rep.syn if rep.feat.ftype == "CDS" else "",
            rep.nonsyn if rep.feat.ftype == "CDS" else "", ",".join(rep.aa_changes),
            ",".join(rep.cds_flags), "%d/%d" % (n_ok, len(lr.results)),
            lr.chain.cid if lr.chain is not None else ""]


FEATURE_HEADER = ["locus_id", "feature_type", "name", "status", "reason", "ref_location",
                  "target_location", "ref_length", "target_length", "nt_identity_pct",
                  "coverage_pct", "mismatches", "insertions", "inserted_bases", "deletions",
                  "deleted_bases", "missing_5prime", "missing_3prime", "aa_identity_pct",
                  "synonymous", "nonsynonymous", "aa_changes", "cds_flags", "notes",
                  "ref_protein", "target_protein"]


def feature_row(fr):
    f = fr.feat
    return [f.locus_id, f.ftype, f.name, fr.status, fr.reason, fmt_loc(f.seqid, f.parts),
            fmt_loc(fr.tgt_id, fr.tgt_parts) if fr.tgt_parts else "-", fr.ref_len,
            fr.tgt_len if fr.tgt_parts else "", _pct(fr.identity) if fr.tgt_parts else "",
            _pct(fr.coverage) if fr.tgt_parts else "", fr.mismatches, fr.ins_events, fr.ins_bases,
            fr.del_events, fr.del_bases, fr.miss5, fr.miss3, _pct(fr.aa_identity),
            fr.syn if f.ftype == "CDS" else "", fr.nonsyn if f.ftype == "CDS" else "",
            ",".join(fr.aa_changes), ",".join(fr.cds_flags), "; ".join(fr.notes),
            fr.ref_prot, fr.tgt_prot]


def target_coverage(chains, tgts):
    cov = {t.sid: np.zeros(t.length, dtype=bool) for t in tgts}
    for c in chains:
        m = c.rmap[c.rmap >= 0]
        if m.size:
            L = c.tgt.length
            q = m % L
            if c.orient == -1:
                q = L - 1 - q
            cov[c.tgt.sid][q] = True
    return cov


def ref_coverage(chains, refs):
    cov = {r.sid: np.zeros(r.length, dtype=bool) for r in refs}
    for c in chains:
        idx = np.nonzero(c.rmap >= 0)[0]
        if idx.size:
            cov[c.ref.sid][(idx + c.R0) % c.ref.length] = True
    return cov


def uncovered_regions(cov, min_len=500):
    out = []
    for sid, arr in cov.items():
        d = np.diff(np.concatenate(([0], (~arr).astype(np.int8), [0])))
        for a, b in zip(np.nonzero(d == 1)[0], np.nonzero(d == -1)[0]):
            if b - a >= min_len:
                out.append((sid, int(a), int(b)))
    return out


# ----------------------------------------------------------------------------
# Figure: genes as arrows along the reference and target genomes
# ----------------------------------------------------------------------------
C_CDS, C_DISRUPT, C_RNA, C_OTHER = "#2f6db5", "#e8871e", "#8452a8", "#3d9a6a"
C_REF_OK, C_FAIL = "#9fb0c4", "#d62728"
C_SYN_SAME, C_SYN_INV = "#9db6d3", "#f2a36b"


def _item_color(lr):
    rep = lr.rep
    if rep.feat.ftype == "CDS":
        return C_DISRUPT if rep.disrupted else C_CDS
    if rep.feat.ftype in RNA_TYPES:
        return C_RNA
    return C_OTHER


def draw_genome_map(png_path, svg_path, refs, tgts, loci_results, chains, ref_name, tgt_name):
    """Two tracks (reference above, target below) with genes drawn as arrows:
    forward-strand genes above the genome line pointing right, reverse-strand
    genes below pointing left; overlapping genes are stacked in lanes.
    Ribbons between the tracks show the orthologous blocks."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import PatchCollection, LineCollection
    from matplotlib.patches import Polygon, Patch
    from matplotlib.ticker import FuncFormatter, MaxNLocator
    from bisect import bisect_left

    FIG_W, LEFT, RIGHT = 16.0, 0.035, 0.99
    UNIT = 0.2                                     # inches per lane unit
    FS = 6.5                                       # label font size
    MAXLAB = 24                                    # max characters of outside labels

    def layout(entries):
        total = sum(e.length for e in entries)
        gap = int(total * min(0.015, 0.15 / len(entries))) if len(entries) > 1 else 0
        offs, x = {}, 0
        for e in entries:
            offs[e.sid] = x
            x += e.length + gap
        return offs, x - gap

    r_off, r_tot = layout(refs)
    t_off, t_tot = layout(tgts)
    xmax = max(r_tot, t_tot, 1)
    rec_len = {e.sid: e.length for e in refs + tgts}
    ax_w_pt = FIG_W * (RIGHT - LEFT) * 72
    dpp = xmax * 1.02 / ax_w_pt                    # data units per point
    char_w = 0.58 * FS * dpp
    pad = 0.5 * dpp                                # stack only genes that (nearly) overlap
    head = xmax * 0.007

    ref_items, tgt_items = [], []
    for lr in loci_results:
        if not lr.locus.is_gene:
            continue
        f = lr.rep.feat
        ref_items.append(dict(
            parts=[(r_off[f.seqid] + s, r_off[f.seqid] + e, st if st else 1) for s, e, st in f.parts],
            label=lr.locus.name, L=rec_len[f.seqid],
            color=C_REF_OK if lr.status == "transferred" else C_FAIL))
        if lr.status == "transferred" and lr.rep.tgt_parts:
            fr = lr.rep
            tgt_items.append(dict(
                parts=[(t_off[fr.tgt_id] + s, t_off[fr.tgt_id] + e, st if st else 1)
                       for s, e, st in fr.tgt_parts],
                label=lr.locus.name, color=_item_color(lr), L=rec_len[fr.tgt_id]))

    def prepare(items):
        """Assign lanes and decide which labels fit inside the arrows."""
        lanes = {1: [], -1: []}
        for it in sorted(items, key=lambda it: min(p[0] for p in it["parts"])):
            strand = it["parts"][0][2]
            it["strand"] = strand
            ivs = [(p[0], p[1] + pad) for p in it["parts"]]
            for li, (starts, ends) in enumerate(lanes[strand]):
                if all(not (bisect_left(starts, e) > 0 and ends[bisect_left(starts, e) - 1] > s)
                       for s, e in ivs):
                    break
            else:
                lanes[strand].append(([], []))
                li = len(lanes[strand]) - 1
            starts, ends = lanes[strand][li]
            for s, e in ivs:
                i = bisect_left(starts, s)
                starts.insert(i, s)
                ends.insert(i, e)
            it["lane"] = li
            longest = max(it["parts"], key=lambda p: p[1] - p[0])
            it["inside"] = len(it["label"]) * char_w < 0.85 * (longest[1] - longest[0])
        g = dict(up=len(lanes[1]), dn=len(lanes[-1]))
        show_out = len(items) <= 250
        for st, key in ((1, "band_up"), (-1, "band_dn")):
            outs = [it for it in items if it["strand"] == st and not it["inside"]] if show_out else []
            longest = max((min(len(it["label"]), MAXLAB) for it in outs), default=0)
            g[key] = (longest * 0.58 * (FS - 0.5) / (UNIT * 72) + 1.0) if outs else 0.0
        g["show_out"] = show_out
        g["top"] = (g["up"] + 0.5 if g["up"] else 0.4) + g["band_up"]
        g["bot"] = -((g["dn"] + 0.5 if g["dn"] else 0.4) + g["band_dn"])
        g["h"] = (g["top"] - g["bot"]) * UNIT
        return g

    g_r, g_t = prepare(ref_items), prepare(tgt_items)
    TITLE, TICKS, RIB, LEG, TOP = 0.32, 0.3, 1.0, 0.45, 0.12
    fig_h = TOP + TITLE + g_r["h"] + TICKS + RIB + TITLE + g_t["h"] + TICKS + LEG
    fig = plt.figure(figsize=(FIG_W, fig_h))

    def frac(y_in):
        return y_in / fig_h

    y = fig_h - TOP - TITLE - g_r["h"]
    ax_r = fig.add_axes([LEFT, frac(y), RIGHT - LEFT, frac(g_r["h"])])
    y_rib_top = y - TICKS
    y = y_rib_top - RIB
    ax_m = fig.add_axes([LEFT, frac(y), RIGHT - LEFT, frac(RIB)])
    y = y - TITLE - g_t["h"]
    ax_t = fig.add_axes([LEFT, frac(y), RIGHT - LEFT, frac(g_t["h"])])

    def draw_track(ax, items, entries, offs, g, title):
        h = 0.72
        patches, colors, lines = [], [], []
        out_labels = {1: [], -1: []}
        for it in items:
            st = it["strand"]
            yc = st * (0.85 + it["lane"])
            parts = sorted(it["parts"], key=lambda p: p[0])
            three_prime = it["parts"][-1]           # biological order: last part is 3'
            for (s, e, _) in parts:
                ln = e - s
                if (s, e) == three_prime[:2]:
                    hl = min(head, ln * 0.6)
                    if st > 0:
                        pts = [(s, yc - h / 2), (e - hl, yc - h / 2), (e, yc), (e - hl, yc + h / 2), (s, yc + h / 2)]
                    else:
                        pts = [(e, yc - h / 2), (s + hl, yc - h / 2), (s, yc), (s + hl, yc + h / 2), (e, yc + h / 2)]
                else:
                    pts = [(s, yc - h * 0.4), (e, yc - h * 0.4), (e, yc + h * 0.4), (s, yc + h * 0.4)]
                patches.append(Polygon(pts, closed=True))
                colors.append(it["color"])
            for (s1, e1, _), (s2, e2, _) in zip(parts[:-1], parts[1:]):
                if s2 - e1 < 0.3 * it["L"]:          # intron line (not across the origin)
                    lines.append([(e1, yc), (s2, yc)])
            longest = max(parts, key=lambda p: p[1] - p[0])
            xc = (longest[0] + longest[1]) / 2
            if it["inside"]:
                ax.text(xc, yc, it["label"], ha="center", va="center", fontsize=FS,
                        color="white", clip_on=True)
            elif g["show_out"]:
                lab = it["label"] if len(it["label"]) <= MAXLAB else it["label"][:MAXLAB - 1] + "…"
                out_labels[st].append((xc, yc, lab))
        ax.add_collection(LineCollection(lines, colors="#555555", linewidths=0.6, zorder=1))
        ax.add_collection(PatchCollection(patches, facecolors=colors, edgecolors="#333333",
                                          linewidths=0.3 if len(items) < 400 else 0, zorder=2))
        # outside labels: vertical text in a band, spread apart, with leader lines
        min_sep = FS * 1.2 * dpp
        for st, labs in out_labels.items():
            if not labs:
                continue
            labs.sort()
            xs = []
            for x, _, _ in labs:
                xs.append(max(x, xs[-1] + min_sep) if xs else x)
            over = xs[-1] - xmax * 1.01
            if over > 0:
                xs = [x - over for x in xs]
                for i in range(len(xs) - 2, -1, -1):
                    xs[i] = min(xs[i], xs[i + 1] - min_sep)
            edge = (g["up"] + 0.5) if st > 0 else -(g["dn"] + 0.5)
            ly = edge + st * 0.6
            leader = []
            for (x, yc, lab), xl in zip(labs, xs):
                leader.append([(x, yc + st * 0.36), (x, edge), (xl, ly - st * 0.1)])
                ax.text(xl, ly, lab, rotation=90, fontsize=FS - 0.5, ha="center",
                        va="bottom" if st > 0 else "top", color="#222222")
            ax.add_collection(LineCollection(leader, colors="#999999", linewidths=0.35, zorder=0))
        for e in entries:
            x0 = offs[e.sid]
            ax.plot([x0, x0 + e.length], [0, 0], color="#444444", lw=1.6,
                    solid_capstyle="butt", zorder=1)
        ax.set_xlim(-xmax * 0.005, xmax * 1.015)
        ax.set_ylim(g["bot"], g["top"])
        ax.set_yticks([])
        for sp in ("left", "right", "top"):
            ax.spines[sp].set_visible(False)
        ax.spines["bottom"].set_color("#aaaaaa")
        ax.set_title(title, loc="left", fontsize=9.5, fontweight="bold", color="#222222", pad=6)
        ax.text(-xmax * 0.004, 0.85 if g["up"] else 0.3, "+", fontsize=9, ha="right",
                va="center", color="#555555")
        ax.text(-xmax * 0.004, -0.85 if g["dn"] else -0.3, "−", fontsize=9, ha="right",
                va="center", color="#555555")
        if len(entries) == 1:
            def kb(x, _pos):
                return ("%.1f Mb" % (x / 1e6)) if xmax >= 2e6 else ("%g kb" % (x / 1e3))
            ax.xaxis.set_major_locator(MaxNLocator(12))
            ax.xaxis.set_major_formatter(FuncFormatter(kb))
        elif len(entries) <= 30:                 # one tick (named) per sequence
            ax.set_xticks([offs[e.sid] + e.length / 2 for e in entries])
            ax.set_xticklabels([e.sid for e in entries], rotation=30 if len(entries) > 6 else 0,
                               ha="right" if len(entries) > 6 else "center")
        else:
            ax.set_xticks([])
            ax.set_xlabel("%d sequences, shown end to end" % len(entries), fontsize=7, color="#444444")
        ax.tick_params(axis="x", labelsize=7, colors="#444444")

    n_ok = sum(1 for lr in loci_results if lr.locus.is_gene and lr.status == "transferred")
    n_all = sum(1 for lr in loci_results if lr.locus.is_gene)
    draw_track(ax_r, ref_items, refs, r_off, g_r,
               "Reference %s — %d genes (grey = transferred, red = failed)" % (ref_name, n_all))
    draw_track(ax_t, tgt_items, tgts, t_off, g_t,
               "Target %s — %d genes transferred" % (tgt_name, n_ok))

    polys, pcol = [], []
    for c in chains:
        for (rs, re_, ts, te) in c.segments():
            if re_ - rs < xmax * 0.001:
                continue
            a, b = r_off[c.ref.sid] + rs, r_off[c.ref.sid] + re_
            x0, x1 = t_off[c.tgt.sid] + ts, t_off[c.tgt.sid] + te
            pts = [(a, 1), (b, 1), (x1, 0), (x0, 0)] if c.orient == 1 else \
                [(a, 1), (b, 1), (x0, 0), (x1, 0)]
            polys.append(Polygon(pts, closed=True))
            pcol.append(C_SYN_SAME if c.orient == 1 else C_SYN_INV)
    ax_m.add_collection(PatchCollection(polys, facecolors=pcol, edgecolors="none", alpha=0.55))
    ax_m.set_xlim(ax_r.get_xlim())
    ax_m.set_ylim(0, 1)
    ax_m.axis("off")

    handles = [Patch(color=C_CDS, label="CDS gene"),
               Patch(color=C_DISRUPT, label="CDS with disruptive change"),
               Patch(color=C_RNA, label="RNA gene"), Patch(color=C_OTHER, label="other gene"),
               Patch(color=C_REF_OK, label="reference: transferred"),
               Patch(color=C_FAIL, label="reference: failed"),
               Patch(color=C_SYN_SAME, alpha=0.6, label="orthologous block"),
               Patch(color=C_SYN_INV, alpha=0.7, label="inverted block")]
    fig.legend(handles=handles, loc="lower center", ncol=8, fontsize=7, frameon=False,
               bbox_to_anchor=(0.5, 0.0))
    fig.savefig(png_path, dpi=200)
    fig.savefig(svg_path)
    plt.close(fig)


# ----------------------------------------------------------------------------
# Summary statistics and reports
# ----------------------------------------------------------------------------
def summarise(loci_results, variants, chains, refs, tgts):
    S = OrderedDict()
    genes = [lr for lr in loci_results if lr.locus.is_gene]
    other = [lr for lr in loci_results if not lr.locus.is_gene]
    ok = [lr for lr in genes if lr.status == "transferred"]
    S["genes_total"] = len(genes)
    S["genes_transferred"] = len(ok)
    S["genes_failed"] = len(genes) - len(ok)
    S["success_rate"] = len(ok) / len(genes) if genes else 0.0
    S["genes_identical"] = sum(1 for lr in ok if lr.rep.is_identical)
    S["genes_with_differences"] = len(ok) - S["genes_identical"]
    cds_ok = [lr for lr in ok if lr.rep.feat.ftype == "CDS"]
    S["cds_transferred"] = len(cds_ok)
    S["cds_disrupted"] = sum(1 for lr in cds_ok if lr.rep.disrupted)
    S["genes_partial"] = sum(1 for lr in ok if lr.rep.miss5 or lr.rep.miss3 or lr.rep.parts_lost)
    S["genes_inverted"] = sum(1 for lr in ok if lr.chain is not None and lr.chain.orient == -1)
    S["by_type"] = Counter((lr.locus.kind, lr.status) for lr in genes)
    S["fail_reasons"] = Counter(re.sub(r"\s*\(.*\)", "", lr.reason) for lr in genes
                                if lr.status != "transferred")
    S["other_total"] = len(other)
    S["other_transferred"] = sum(1 for lr in other if lr.status == "transferred")
    S["features_total"] = sum(len(lr.results) for lr in loci_results)
    S["features_transferred"] = sum(1 for lr in loci_results for r in lr.results
                                    if r.status == "transferred")
    S["ref_length"] = sum(r.length for r in refs)
    S["tgt_length"] = sum(t.length for t in tgts)
    rc, tc = ref_coverage(chains, refs), target_coverage(chains, tgts)
    S["ref_aligned"] = sum(int(a.sum()) for a in rc.values())
    S["tgt_aligned"] = sum(int(a.sum()) for a in tc.values())
    n_al = sum(c.n_aligned for c in chains)
    S["block_identity"] = (sum(c.identity * c.n_aligned for c in chains) / n_al) if n_al else 0.0
    S["blocks"] = len(chains)
    S["blocks_inverted"] = sum(1 for c in chains if c.orient == -1)
    S["uncovered_target"] = uncovered_regions(tc)
    vt = Counter(v["type"] for v in variants)
    S["snps"], S["insertions"], S["deletions"] = vt.get("SNP", 0), vt.get("INS", 0), vt.get("DEL", 0)
    S["complex"], S["unaligned"] = vt.get("COMPLEX", 0), vt.get("UNALIGNED", 0)
    eff = Counter()
    for v in variants:
        e = v.get("effect", "")
        if e in ("intergenic", "non-coding"):
            eff[e] += 1
            continue
        kinds = set()
        for x in e.split(";"):
            x = x.split(":", 1)[-1]
            if x.startswith("synonymous"):
                kinds.add("synonymous")
            elif x in ("frameshift", "inframe_deletion", "inframe_insertion", "complex",
                       "inframe_complex", "unaligned_region"):
                kinds.add(x)
            elif x.endswith("*"):
                kinds.add("nonsense")
            elif re.match(r"^\*\d+", x):
                kinds.add("stop_lost")
            else:
                kinds.add("missense")
        for k in kinds:
            eff[k] += 1
    S["variant_effects"] = eff
    return S


def write_summary_text(path, S, loci_results, args, files):
    L = []
    L.append("Annotation transfer summary  (transfer_annotation.py v%s)" % __version__)
    L.append("=" * 72)
    L.append("Reference : %s" % args.reference)
    L.append("Target    : %s" % args.target)
    L.append("")
    L.append("GENES")
    L.append("  reference genes            : %d" % S["genes_total"])
    L.append("  transferred successfully   : %d (%.1f%%)" % (S["genes_transferred"], 100 * S["success_rate"]))
    L.append("  failed                     : %d" % S["genes_failed"])
    L.append("    identical to reference   : %d" % S["genes_identical"])
    L.append("    with differences         : %d" % S["genes_with_differences"])
    L.append("    partial / truncated      : %d" % S["genes_partial"])
    L.append("    CDS with disruptive change: %d of %d CDS" % (S["cds_disrupted"], S["cds_transferred"]))
    if S["fail_reasons"]:
        L.append("  failure reasons:")
        for r, n in S["fail_reasons"].most_common():
            L.append("    %-40s %d" % (r, n))
    L.append("  other features (non-gene)  : %d of %d transferred" % (S["other_transferred"], S["other_total"]))
    L.append("  all features               : %d of %d transferred" % (S["features_transferred"], S["features_total"]))
    L.append("")
    L.append("GENOME ALIGNMENT")
    L.append("  reference length           : %d bp (%.1f%% in orthologous blocks)" % (
        S["ref_length"], 100 * S["ref_aligned"] / max(1, S["ref_length"])))
    L.append("  target length              : %d bp (%.1f%% in orthologous blocks)" % (
        S["tgt_length"], 100 * S["tgt_aligned"] / max(1, S["tgt_length"])))
    L.append("  orthologous blocks         : %d (%d inverted)" % (S["blocks"], S["blocks_inverted"]))
    L.append("  identity in blocks         : %.2f%%" % (100 * S["block_identity"]))
    L.append("  SNPs / insertions / deletions: %d / %d / %d" % (S["snps"], S["insertions"], S["deletions"]))
    if S["complex"] or S["unaligned"]:
        L.append("  complex / unaligned regions: %d / %d" % (S["complex"], S["unaligned"]))
    if S["uncovered_target"]:
        L.append("  target regions >=500 bp without reference counterpart: %d" % len(S["uncovered_target"]))
    L.append("")
    failed = [lr for lr in loci_results if lr.locus.is_gene and lr.status != "transferred"]
    if failed:
        L.append("FAILED GENES")
        for lr in failed:
            L.append("  %-18s %-6s %-36s %s" % (lr.locus.name[:18], lr.locus.kind,
                                                fmt_loc(lr.locus.seqid, lr.locus.primary.parts)[:36], lr.reason))
        L.append("")
    diff = [lr for lr in loci_results if lr.locus.is_gene and lr.status == "transferred"
            and not lr.rep.is_identical]
    if diff:
        L.append("GENES WITH DIFFERENCES (representative CDS/RNA)")
        L.append("  %-18s %8s %8s %7s %6s %5s %5s  %s" % ("gene", "ref_len", "tgt_len", "nt_id%",
                                                          "aa_id%", "nsyn", "syn", "flags / aa changes"))
        for lr in diff:
            r = lr.rep
            extra = ",".join(r.cds_flags)
            ch = r.aa_changes[:8]
            if ch:
                extra += (" | " if extra else "") + ",".join(ch) + (" ..." if len(r.aa_changes) > 8 else "")
            L.append("  %-18s %8d %8d %7.2f %6s %5s %5s  %s" % (
                lr.locus.name[:18], r.ref_len, r.tgt_len, 100 * r.identity,
                ("%.2f" % (100 * r.aa_identity)) if r.aa_identity is not None else "-",
                r.nonsyn if r.feat.ftype == "CDS" else "-", r.syn if r.feat.ftype == "CDS" else "-", extra))
        L.append("")
    L.append("OUTPUT FILES")
    for desc, fn in files:
        L.append("  %-28s %s" % (desc, fn))
    txt = "\n".join(L) + "\n"
    with open(path, "w") as fh:
        fh.write(txt)
    return txt


_CSS = """
:root{--fg:#1d2430;--muted:#5d6878;--line:#dfe4ea;--bg:#ffffff;--soft:#f5f7fa;--ok:#2e7d4f;
--bad:#c62828;--warn:#d9822b;--acc:#2f6db5}
*{box-sizing:border-box}body{margin:0;font:14px/1.45 -apple-system,BlinkMacSystemFont,"Segoe UI",
Roboto,Helvetica,Arial,sans-serif;color:var(--fg);background:var(--bg)}
main{max-width:1280px;margin:0 auto;padding:28px 20px 60px}
h1{font-size:24px;margin:0 0 4px}h2{font-size:18px;margin:34px 0 10px;padding-bottom:6px;
border-bottom:1px solid var(--line)}.sub{color:var(--muted);font-size:13px}
.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:12px;margin:18px 0}
.card{background:var(--soft);border:1px solid var(--line);border-radius:10px;padding:12px 14px}
.card .v{font-size:26px;font-weight:650}.card .l{color:var(--muted);font-size:12.5px}
.ok{color:var(--ok)}.bad{color:var(--bad)}.warn{color:var(--warn)}
table{border-collapse:collapse;width:100%;font-size:12.5px}th,td{padding:5px 8px;
border-bottom:1px solid var(--line);text-align:left;vertical-align:top}
th{background:var(--soft);position:sticky;top:0;cursor:pointer;white-space:nowrap}
td.n{text-align:right;font-variant-numeric:tabular-nums}.wrap{max-height:620px;overflow:auto;
border:1px solid var(--line);border-radius:8px}
.kv td:first-child{color:var(--muted);width:290px}.fig{width:100%;border:1px solid var(--line);
border-radius:8px}.tag{display:inline-block;padding:1px 6px;border-radius:4px;font-size:11.5px;
margin:1px 2px 1px 0;background:#fdeee0;color:#8a4b0f}.tag.r{background:#fde7e7;color:#9b1c1c}
.ctrl{margin:8px 0;display:flex;gap:14px;align-items:center;flex-wrap:wrap}
input[type=search]{padding:6px 9px;border:1px solid var(--line);border-radius:6px;min-width:260px}
.mono{font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;font-size:12px}
details summary{cursor:pointer;color:var(--acc)}
td.loc{font-size:11px;word-break:break-all;min-width:190px;max-width:260px}
td.aa{font-size:11px;word-break:break-word;min-width:200px;max-width:380px}
#gt{min-width:1100px}
"""

_JS = """
function filt(id,q,onlyDiff){q=q.toLowerCase();document.querySelectorAll('#'+id+' tbody tr').forEach(r=>{
let s=r.textContent.toLowerCase().includes(q);if(onlyDiff&&r.dataset.diff==='0')s=false;
r.style.display=s?'':'none'})}
function sortT(th){const t=th.closest('table'),i=[...th.parentNode.children].indexOf(th),
b=t.tBodies[0],rows=[...b.rows],asc=th.dataset.asc!=='1';th.dataset.asc=asc?'1':'0';
rows.sort((a,c)=>{let x=a.cells[i].dataset.v??a.cells[i].textContent,y=c.cells[i].dataset.v??c.cells[i].textContent;
let nx=parseFloat(x),ny=parseFloat(y);if(!isNaN(nx)&&!isNaN(ny))return asc?nx-ny:ny-nx;
return asc?x.localeCompare(y):y.localeCompare(x)});rows.forEach(r=>b.appendChild(r))}
document.querySelectorAll('th').forEach(th=>th.addEventListener('click',()=>sortT(th)));
"""


def _e(x):
    return html.escape(str(x))


def write_html_report(path, S, loci_results, variants, chains, args, png_path, files, refs, tgts):
    with open(png_path, "rb") as fh:
        img64 = base64.b64encode(fh.read()).decode()
    genes = [lr for lr in loci_results if lr.locus.is_gene]
    failed = [lr for lr in genes if lr.status != "transferred"]
    ok = [lr for lr in genes if lr.status == "transferred"]
    H = []
    H.append("<!DOCTYPE html><html lang='en'><head><meta charset='utf-8'>"
             "<meta name='viewport' content='width=device-width,initial-scale=1'>"
             "<title>Annotation transfer report</title><style>%s</style></head><body><main>" % _CSS)
    H.append("<h1>Annotation transfer report</h1>")
    H.append("<div class='sub'>Target <b>%s</b> annotated from reference <b>%s</b> &middot; %s &middot; "
             "transfer_annotation.py v%s</div>" % (
                 _e(os.path.basename(args.target)), _e(os.path.basename(args.reference)),
                 datetime.datetime.now().strftime("%Y-%m-%d %H:%M"), __version__))
    rate_cls = "ok" if S["success_rate"] >= 0.9 else ("warn" if S["success_rate"] >= 0.5 else "bad")
    cards = [
        (S["genes_total"], "genes in reference", ""),
        (S["genes_transferred"], "genes transferred", "ok"),
        (S["genes_failed"], "genes failed", "bad" if S["genes_failed"] else ""),
        ("%.1f%%" % (100 * S["success_rate"]), "success rate", rate_cls),
        (S["genes_identical"], "identical to reference", ""),
        (S["genes_with_differences"], "with differences", ""),
        ("%d / %d" % (S["cds_disrupted"], S["cds_transferred"]), "CDS with disruptive change",
         "warn" if S["cds_disrupted"] else ""),
        ("%d / %d" % (S["other_transferred"], S["other_total"]), "other features transferred", ""),
    ]
    H.append("<div class='cards'>")
    for v, l, c in cards:
        H.append("<div class='card'><div class='v %s'>%s</div><div class='l'>%s</div></div>" % (c, _e(v), _e(l)))
    H.append("</div>")

    H.append("<h2>Genome map</h2>")
    H.append("<div class='sub'>Genes are drawn as arrows: forward-strand genes above each genome "
             "line pointing right, reverse-strand genes below pointing left. Ribbons link orthologous "
             "blocks (orange = inverted).</div>")
    H.append("<img class='fig' alt='genome map' src='data:image/png;base64,%s'>" % img64)

    H.append("<h2>Genome alignment</h2><table class='kv'>")
    kv = [
        ("Reference", "%s (%d sequence%s, %s bp)" % (_e(", ".join(r.sid for r in refs[:6])) +
                                                     (" ..." if len(refs) > 6 else ""), len(refs),
                                                     "s" if len(refs) > 1 else "", "{:,}".format(S["ref_length"]))),
        ("Target", "%s (%d sequence%s, %s bp)" % (_e(", ".join(t.sid for t in tgts[:6])) +
                                                  (" ..." if len(tgts) > 6 else ""), len(tgts),
                                                  "s" if len(tgts) > 1 else "", "{:,}".format(S["tgt_length"]))),
        ("Reference in orthologous blocks", "%.2f%%" % (100 * S["ref_aligned"] / max(1, S["ref_length"]))),
        ("Target in orthologous blocks", "%.2f%%" % (100 * S["tgt_aligned"] / max(1, S["tgt_length"]))),
        ("Orthologous blocks", "%d (%d inverted)" % (S["blocks"], S["blocks_inverted"])),
        ("Nucleotide identity in blocks", "%.2f%%" % (100 * S["block_identity"])),
        ("SNPs / insertions / deletions", "{:,} / {:,} / {:,}".format(S["snps"], S["insertions"], S["deletions"])),
        ("Complex substitutions / unaligned regions", "{:,} / {:,}".format(S["complex"], S["unaligned"])),
        ("Coding effects of variants", ", ".join("%s %d" % (k, v) for k, v in S["variant_effects"].most_common()) or "-"),
        ("Genes on inverted blocks", S["genes_inverted"]),
        ("Partial / truncated genes", S["genes_partial"]),
    ]
    if S["uncovered_target"]:
        regs = ", ".join("%s:%d-%d" % (s, a + 1, b) for s, a, b in S["uncovered_target"][:12])
        kv.append(("Target regions &ge;500 bp absent from reference",
                   "%d &mdash; <span class='mono'>%s%s</span>" % (len(S["uncovered_target"]), _e(regs),
                                                                  " ..." if len(S["uncovered_target"]) > 12 else "")))
    for k, v in kv:
        H.append("<tr><td>%s</td><td>%s</td></tr>" % (k, v))
    H.append("</table>")

    H.append("<h2>Transfer by gene type</h2><table><thead><tr><th>type</th><th>in reference</th>"
             "<th>transferred</th><th>failed</th></tr></thead><tbody>")
    kinds = sorted({k for k, _ in S["by_type"]})
    for k in kinds:
        t = S["by_type"].get((k, "transferred"), 0)
        f = S["by_type"].get((k, "failed"), 0)
        H.append("<tr><td>%s</td><td class='n'>%d</td><td class='n'>%d</td><td class='n'>%d</td></tr>"
                 % (_e(k), t + f, t, f))
    H.append("</tbody></table>")

    H.append("<h2>Failed genes (%d)</h2>" % len(failed))
    if failed:
        H.append("<div class='wrap'><table><thead><tr><th>gene</th><th>locus tag</th><th>type</th>"
                 "<th>product</th><th>reference location</th><th>length</th><th>best coverage %</th>"
                 "<th>best identity %</th><th>reason</th></tr></thead><tbody>")
        for lr in failed:
            r = lr.rep
            H.append("<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td class='mono'>%s</td>"
                     "<td class='n'>%d</td><td class='n'>%s</td><td class='n'>%s</td><td>%s</td></tr>" % (
                         _e(lr.locus.name), _e(lr.locus.locus_tag), _e(lr.locus.kind), _e(lr.locus.product),
                         _e(fmt_loc(lr.locus.seqid, lr.span.feat.parts)), r.ref_len,
                         ("%.1f" % (100 * r.coverage)) if r.tgt_parts else "0",
                         ("%.1f" % (100 * r.identity)) if r.tgt_parts else "-", _e(lr.reason)))
        H.append("</tbody></table></div>")
    else:
        H.append("<p>All reference genes were transferred.</p>")

    H.append("<h2>Differences between reference and target genes (%d transferred)</h2>" % len(ok))
    H.append("<div class='sub'>Metrics are for the gene product (CDS or RNA). nt identity = identical "
             "bases / alignment columns; coverage = reference bases aligned. Amino-acid changes use "
             "reference numbering. Click a column header to sort.</div>")
    H.append("<div class='ctrl'><input type='search' placeholder='Filter genes...' "
             "oninput=\"filt('gt',this.value,document.getElementById('od').checked)\">"
             "<label><input type='checkbox' id='od' checked onchange=\"filt('gt',"
             "document.querySelector('input[type=search]').value,this.checked)\"> only genes with "
             "differences (%d)</label></div>" % S["genes_with_differences"])
    H.append("<div class='wrap'><table id='gt'><thead><tr><th>gene</th><th>type</th>"
             "<th>location (reference / target)</th><th>length<br>ref / target</th>"
             "<th>nt id %</th><th>cov %</th><th>mism.</th><th>ins / del<br>events (bp)</th>"
             "<th>aa id %</th><th>non-syn / syn</th><th>flags</th><th>amino-acid changes</th>"
             "</tr></thead><tbody>")
    for lr in ok:
        r = lr.rep
        diff = 0 if r.is_identical else 1
        flags = "".join("<span class='tag%s'>%s</span>" % (" r" if fl.startswith(("frameshift", "premature"))
                                                            else "", _e(fl)) for fl in r.cds_flags)
        if lr.chain is not None and lr.chain.orient == -1:
            flags += "<span class='tag'>inverted</span>"
        if r.parts_lost:
            flags += "<span class='tag'>%d part(s) lost</span>" % r.parts_lost
        for nt in r.notes:
            if "terminal deletion" in nt:
                flags += "<span class='tag'>%s</span>" % _e(nt)
        ch = r.aa_changes
        ch_txt = _e(", ".join(ch[:20])) + (" &hellip; (+%d)" % (len(ch) - 20) if len(ch) > 20 else "")
        is_cds = r.feat.ftype == "CDS"
        indel = "%d (%d) / %d (%d)" % (r.ins_events, r.ins_bases, r.del_events, r.del_bases)
        H.append("<tr data-diff='%d'%s><td><b>%s</b><div class='sub'>%s</div></td><td>%s</td>"
                 "<td class='mono loc'>%s<br>%s</td><td class='n' data-v='%d'>%d / %d<div class='sub'>%+d</div></td>"
                 "<td class='n'>%.2f</td><td class='n'>%.1f</td><td class='n'>%d</td>"
                 "<td class='n' data-v='%d'>%s</td><td class='n'>%s</td><td class='n' data-v='%s'>%s</td>"
                 "<td>%s</td><td class='mono aa'>%s</td></tr>" % (
                     diff, "" if diff else " style='display:none'", _e(lr.locus.name),
                     _e(lr.locus.product[:60]), _e(lr.locus.kind),
                     _e(fmt_loc(lr.locus.seqid, r.feat.parts)), _e(fmt_loc(r.tgt_id, r.tgt_parts)),
                     r.tgt_len - r.ref_len, r.ref_len, r.tgt_len, r.tgt_len - r.ref_len,
                     100 * r.identity, 100 * r.coverage, r.mismatches,
                     r.ins_bases + r.del_bases, indel,
                     ("%.2f" % (100 * r.aa_identity)) if r.aa_identity is not None else "-",
                     r.nonsyn if is_cds else -1, ("%d / %d" % (r.nonsyn, r.syn)) if is_cds else "-",
                     flags, ch_txt))
    H.append("</tbody></table></div>")

    H.append("<h2>Parameters and output files</h2><table class='kv'>")
    params = [("min identity", args.min_identity), ("min coverage", args.min_coverage),
              ("k-mer size", args.kmer_used), ("max gap in block", args.max_gap),
              ("min block score", args.min_chain_score), ("default translation table", args.transl_table),
              ("target topology", args.topology)]
    for k, v in params:
        H.append("<tr><td>%s</td><td>%s</td></tr>" % (_e(k), _e(v)))
    for desc, fn in files:
        H.append("<tr><td>%s</td><td class='mono'>%s</td></tr>" % (_e(desc), _e(os.path.basename(fn))))
    H.append("</table></main><script>%s</script></body></html>" % _JS)
    with open(path, "w") as fh:
        fh.write("\n".join(H))


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def parse_args(argv=None):
    p = argparse.ArgumentParser(
        prog="transfer_annotation.py",
        description="Transfer genome annotation from a reference (GenBank or GFF3) to a target "
                    "genome (FASTA) through orthologous syntenic regions, and report the results.",
        epilog="example: python transfer_annotation.py -r ref.gb -t new_strain.fasta -o out")
    io_ = p.add_argument_group("input / output")
    io_.add_argument("-r", "--reference", required=True,
                     help="reference annotation: GenBank (.gb/.gbk) or GFF3")
    io_.add_argument("-t", "--target", required=True,
                     help="target genome to annotate (FASTA, one or more sequences)")
    io_.add_argument("--ref-fasta", help="reference sequence (FASTA) if the GFF3 has no ##FASTA section")
    io_.add_argument("-o", "--outdir", default="annotation_transfer",
                     help="output directory [%(default)s]")
    io_.add_argument("-p", "--prefix", help="output file prefix [target file name]")
    tr = p.add_argument_group("transfer")
    tr.add_argument("--min-identity", type=float, default=0.5,
                    help="minimum nucleotide identity of a transferred feature, 0-1 [%(default)s]")
    tr.add_argument("--min-coverage", type=float, default=0.5,
                    help="minimum fraction of the reference feature that must align, 0-1 [%(default)s]")
    tr.add_argument("--types", help="only transfer these feature types (comma separated)")
    tr.add_argument("--exclude-types", default="",
                    help="feature types not to transfer (comma separated), e.g. variation,misc_difference")
    tr.add_argument("--transl-table", type=int, default=1,
                    help="genetic code for CDS without a /transl_table qualifier [%(default)s]")
    tr.add_argument("--topology", choices=["auto", "linear", "circular"], default="auto",
                    help="target topology; auto = circular when the reference is circular and the "
                         "target is a complete genome (same number of sequences, similar length) "
                         "[%(default)s]")
    al = p.add_argument_group("orthologous blocks")
    al.add_argument("-k", "--kmer", type=int, default=0,
                    help="anchor k-mer size, 9-31; 0 = automatic from genome size [%(default)s]")
    al.add_argument("--max-gap", type=int, default=30000,
                    help="maximum gap between anchors inside one block, bp [%(default)s]")
    al.add_argument("--min-chain-score", type=float, default=50,
                    help="minimum block score (about the number of matched bases) [%(default)s]")
    ot = p.add_argument_group("output options")
    ot.add_argument("--copy-organism", action="store_true",
                    help="copy /organism from the reference into the target GenBank")
    ot.add_argument("--gff-fasta", action="store_true", help="append the target sequence (##FASTA) to the GFF3")
    ot.add_argument("--no-figure", action="store_true", help="skip the genome map figure and HTML report")
    ot.add_argument("-q", "--quiet", action="store_true", help="only print warnings and the summary")
    p.add_argument("--version", action="version", version="%(prog)s " + __version__)
    a = p.parse_args(argv)
    if not (0 <= a.min_identity <= 1 and 0 <= a.min_coverage <= 1):
        p.error("--min-identity and --min-coverage must be between 0 and 1")
    if a.kmer and not (9 <= a.kmer <= 31):
        p.error("--kmer must be between 9 and 31")
    return a


def main(argv=None):
    args = parse_args(argv)
    logging.basicConfig(level=logging.WARNING if args.quiet else logging.INFO,
                        format="[%(asctime)s] %(message)s", datefmt="%H:%M:%S")
    t_start = time.time()
    for fn in (args.reference, args.target) + ((args.ref_fasta,) if args.ref_fasta else ()):
        if not os.path.isfile(fn):
            raise SystemExit("ERROR: file not found: %s" % fn)
    os.makedirs(args.outdir, exist_ok=True)
    prefix = args.prefix or re.sub(r"\.(fa|fasta|fna|fas|fsa|seq|txt)(\.gz)?$", "",
                                   os.path.basename(args.target), flags=re.I)
    out = lambda suffix: os.path.join(args.outdir, prefix + suffix)

    fmt = detect_format(args.reference)
    LOG.info("Reading reference %s (%s)", args.reference, fmt.upper())
    refs, feats = read_genbank(args.reference) if fmt == "genbank" else read_gff3(args.reference, args.ref_fasta)
    if args.types:
        keep = {x.strip() for x in args.types.split(",") if x.strip()}
        feats = [f for f in feats if f.ftype in keep]
    excl = {x.strip() for x in args.exclude_types.split(",") if x.strip()}
    if excl:
        feats = [f for f in feats if f.ftype not in excl]
    if not feats:
        raise SystemExit("ERROR: no features to transfer in the reference")
    LOG.info("  %d sequence(s), %s bp, %d features", len(refs), "{:,}".format(sum(r.length for r in refs)), len(feats))
    LOG.info("Reading target %s", args.target)
    tgts = read_fasta(args.target)
    if args.topology == "auto":
        # circular only for a complete genome of a circular reference (not for draft contigs)
        ratio = sum(t.length for t in tgts) / max(1, sum(r.length for r in refs))
        tgt_circ = all(r.circular for r in refs) and len(tgts) == len(refs) and 0.8 <= ratio <= 1.25
        args.topology = "circular (auto)" if tgt_circ else "linear (auto)"
    else:
        tgt_circ = args.topology == "circular"
    for t in tgts:
        t.circular = tgt_circ
    LOG.info("  %d sequence(s), %s bp", len(tgts), "{:,}".format(sum(t.length for t in tgts)))

    loci = build_loci(feats, fmt)
    LOG.info("  %d loci (%d genes, %d other features)", len(loci),
             sum(1 for l in loci if l.is_gene), sum(1 for l in loci if not l.is_gene))

    k = args.kmer or auto_kmer(max(sum(r.length for r in refs), sum(t.length for t in tgts)))
    args.kmer_used = k
    LOG.info("Finding orthologous blocks (k=%d)", k)
    chains = build_blocks(refs, tgts, k, args.max_gap, args.min_chain_score)
    if not chains:
        LOG.warning("No orthologous blocks found - are reference and target related? "
                    "Try a smaller --kmer.")

    LOG.info("Transferring annotation")
    lifter = Lifter(refs, tgts, chains, args)
    results = [lifter.lift_locus(l) for l in loci]
    variants = call_variants(chains, refs)
    annotate_variants(variants, lifter, results)

    files = [("annotated GenBank", out(".gbk")), ("annotated GFF3", out(".gff3")),
             ("gene table", out("_genes.tsv")), ("feature table", out("_features.tsv")),
             ("variants", out("_variants.tsv")), ("orthologous blocks", out("_blocks.tsv"))]
    LOG.info("Writing outputs to %s", args.outdir)
    ref_by_id = OrderedDict((r.sid, r) for r in refs)
    write_genbank(out(".gbk"), tgts, results, ref_by_id, args)
    write_gff3(out(".gff3"), tgts, results, args.gff_fasta)
    write_tsv(out("_genes.tsv"), GENE_HEADER, [gene_row(lr) for lr in results if lr.locus.is_gene])
    write_tsv(out("_features.tsv"), FEATURE_HEADER, [feature_row(fr) for lr in results for fr in lr.results])
    write_tsv(out("_variants.tsv"),
              ["ref_seqid", "ref_pos", "ref_allele", "alt_allele", "type", "length", "target_seqid",
               "target_pos", "block", "genes", "effect"],
              [[v["ref"], v["pos"], v["ref_allele"], v["alt"], v["type"], v.get("length", 1),
                v["tgt"], v["tpos"], v["chain"], v["genes"], v["effect"]] for v in variants])
    blk_rows = []
    for c in chains:
        for (rs, re_, ts, te) in c.segments():
            blk_rows.append([c.cid, c.ref.sid, rs + 1, re_, c.tgt.sid, ts + 1, te,
                             "+" if c.orient == 1 else "-", "%.2f" % (100 * c.identity), int(c.score)])
    write_tsv(out("_blocks.tsv"), ["block", "ref_seqid", "ref_start", "ref_end", "target_seqid",
                                   "target_start", "target_end", "orientation", "identity_pct",
                                   "anchor_score"], blk_rows)
    S = summarise(results, variants, chains, refs, tgts)
    if not args.no_figure:
        LOG.info("Drawing genome map")
        ref_name = refs[0].sid if len(refs) == 1 else os.path.basename(args.reference)
        tgt_name = tgts[0].sid if len(tgts) == 1 else os.path.basename(args.target)
        draw_genome_map(out("_genome_map.png"), out("_genome_map.svg"), refs, tgts, results, chains,
                        ref_name, tgt_name)
        files += [("genome map (PNG)", out("_genome_map.png")), ("genome map (SVG)", out("_genome_map.svg"))]
    files += [("HTML report", out("_report.html")), ("text summary", out("_summary.txt"))]
    if not args.no_figure:
        write_html_report(out("_report.html"), S, results, variants, chains, args,
                          out("_genome_map.png"), files, refs, tgts)
    else:
        files = [f for f in files if f[0] != "HTML report"]
    txt = write_summary_text(out("_summary.txt"), S, results, args, files)
    lines = txt.splitlines()
    if len(lines) > 90:                       # keep the terminal readable for big genomes
        head_part = [l for l in lines[:70]]
        tail_start = next((i for i, l in enumerate(lines) if l.startswith("OUTPUT FILES")), len(lines))
        lines = head_part + ["  ... (%d more lines - full list in %s)" % (tail_start - 70, out("_summary.txt")), ""] \
            + lines[tail_start:]
    print("\n".join(lines))
    LOG.info("Done in %.1f s", time.time() - t_start)
    return 0


if __name__ == "__main__":
    sys.exit(main())
