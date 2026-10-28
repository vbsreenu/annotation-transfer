#!/usr/bin/env python3
"""Simulate a target genome from a reference GenBank with known edits, run
transfer_annotation.py, and verify every transferred feature location against
the truth coordinate map."""
import argparse, random, subprocess, sys, os, json

HERE = os.path.dirname(os.path.abspath(__file__))
PROGRAM = os.path.join(HERE, "..", "transfer_annotation.py")
from Bio import SeqIO
from Bio.Seq import Seq

COMP = str.maketrans("ACGTN", "TGCAN")


def simulate(ref_seq, ops, rng, protected=()):
    # each target base: [base, ref_pos or None, strand(+1/-1)]
    t = [[b, i, 1] for i, b in enumerate(ref_seq)]
    log = []
    for op in ops:
        kind = op[0]
        if kind == "snp":
            rate = op[1]
            n = int(len(t) * rate)
            for p in rng.sample(range(len(t)), n):
                b = t[p][0]
                t[p][0] = rng.choice([x for x in "ACGT" if x != b])
            log.append("snps:%d" % n)
        elif kind == "indel":
            for _ in range(op[1]):
                p = rng.randrange(10, len(t) - 10)
                L = rng.randint(1, op[2])
                if rng.random() < 0.5:
                    del t[p:p + L]
                else:
                    t[p:p] = [[rng.choice("ACGT"), None, 1] for _ in range(L)]
            log.append("indels:%d" % op[1])
        elif kind == "delete":              # delete ref interval [s,e)
            s, e = op[1], op[2]
            t = [x for x in t if x[1] is None or not (s <= x[1] < e)]
            log.append("delete:%d-%d" % (s, e))
        elif kind == "insert":              # insert random seq after ref pos
            pos, L = op[1], op[2]
            idx = next(i for i, x in enumerate(t) if x[1] == pos)
            t[idx + 1:idx + 1] = [[rng.choice("ACGT"), None, 1] for _ in range(L)]
            log.append("insert:%d@%d" % (L, pos))
        elif kind == "invert":              # invert target segment covering ref [s,e)
            s, e = op[1], op[2]
            i0 = next(i for i, x in enumerate(t) if x[1] is not None and x[1] >= s)
            i1 = max(i for i, x in enumerate(t) if x[1] is not None and x[1] < e) + 1
            seg = t[i0:i1][::-1]
            for x in seg:
                x[0] = x[0].translate(COMP)
                x[2] = -x[2]
            t[i0:i1] = seg
            log.append("invert:%d-%d" % (s, e))
        elif kind == "rotate":
            r = op[1] % len(t)
            t = t[r:] + t[:r]
            log.append("rotate:%d" % r)
        elif kind == "revcomp":
            t = t[::-1]
            for x in t:
                x[0] = x[0].translate(COMP)
                x[2] = -x[2]
            log.append("revcomp")
    seq = "".join(x[0] for x in t)
    truth = {}
    for i, x in enumerate(t):
        if x[1] is not None:
            truth[x[1]] = (i, x[2])
    return seq, truth, log


def expected_parts(parts, truth, tlen):
    """Expected target parts for ref parts via truth map (None if deleted)."""
    out = []
    for s, e, st in parts:
        pos = [truth[p] for p in range(s, e) if p in truth]
        if not pos:
            out.append(None)
            continue
        strands = {o for _, o in pos}
        if len(strands) > 1:
            out.append("split")
            continue
        o = strands.pop()
        ts = [q for q, _ in pos]
        lo, hi = min(ts), max(ts) + 1
        if hi - lo > 0.5 * tlen and len(ts) < 0.5 * tlen:   # wraps the target origin
            out.append(("wrap", st * o))
            continue
        out.append((lo, hi, st * o))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ref")
    ap.add_argument("outdir")
    ap.add_argument("--ops", required=True, help="json list of ops")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--tol", type=int, default=6)
    ap.add_argument("--split", type=int, default=0, help="split target into N contigs")
    ap.add_argument("--extra", default="", help="extra args for transfer_annotation.py")
    ap.add_argument("--ref-arg", default=None, help="override -r (e.g. gff3) keeping truth from genbank")
    a = ap.parse_args()
    rng = random.Random(a.seed)
    rec = SeqIO.read(a.ref, "genbank")
    ref_seq = str(rec.seq).upper()
    seq, truth, log = simulate(ref_seq, json.loads(a.ops), rng)
    os.makedirs(a.outdir, exist_ok=True)
    tfa = os.path.join(a.outdir, "target.fasta")
    contig_of = None
    with open(tfa, "w") as fh:
        if a.split:
            cuts = sorted(rng.sample(range(1000, len(seq) - 1000), a.split - 1))
            bounds = [0] + cuts + [len(seq)]
            for i in range(a.split):
                fh.write(">contig%d\n%s\n" % (i + 1, seq[bounds[i]:bounds[i + 1]]))
            contig_of = bounds
        else:
            fh.write(">target\n%s\n" % seq)
    print("simulation:", ", ".join(log), "| target length", len(seq))
    cmd = [sys.executable, PROGRAM, "-r", a.ref_arg or a.ref, "-t", tfa,
           "-o", a.outdir, "-p", "out", "-q"] + a.extra.split()
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print(res.stdout[-3000:], res.stderr[-5000:])
        sys.exit("RUN FAILED")
    summary = res.stdout
    # parse output genbank
    out_feats = []
    for orec in SeqIO.parse(os.path.join(a.outdir, "out.gbk"), "genbank"):
        off = 0
        if contig_of:
            off = contig_of[int(orec.id.replace("contig", "")) - 1]
        for f in orec.features:
            if f.type == "source":
                continue
            note = " ".join(f.qualifiers.get("note", []))
            out_feats.append((f, off, orec, note))
    # index reference features by location string in note
    import re
    n_ok = n_bad = n_skip = 0
    bad = []
    ref_feats = [f for f in rec.features if f.type != "source"]
    from_note = {}
    for f, off, orec, note in out_feats:
        m = re.search(r"transferred from (.+?) \(nt identity", note)
        if m:
            from_note.setdefault(re.sub(r"\s+", "", m.group(1)), []).append((f, off, orec))
    sys.path.insert(0, os.path.dirname(PROGRAM))
    import transfer_annotation as ta
    seen_missing = []
    tgt_from_ref = {q: (p, o) for p, (q, o) in truth.items()}
    for rf in ref_feats:
        parts = [(int(p.start), int(p.end), p.strand or 0) for p in rf.location.parts]
        key = ta.fmt_loc(rec.id, parts)
        hits = from_note.get(re.sub(r"\s+", "", key))
        exp_set, exp_strands = set(), set()
        for s0, e0, st in parts:
            for p in range(s0, e0):
                if p in truth:
                    exp_set.add(truth[p][0]); exp_strands.add((st or 1) * truth[p][1])
        if not hits:
            seen_missing.append((rf.type, rf.qualifiers.get("gene", rf.qualifiers.get("product", ["?"]))[0], key, len(exp_set)))
            continue
        hit = [h for h in hits if h[0].type == rf.type]
        if not hit:
            n_skip += 1
            continue
        of, off, orec = hit[0]
        got_set, got_strands = set(), set()
        for prt in of.location.parts:
            for q in range(int(prt.start) + off, int(prt.end) + off):
                if q in tgt_from_ref:          # ignore bases inserted in the target
                    got_set.add(q)
            got_strands.add(prt.strand or 1)
        native_mixed = len({p[2] for p in parts}) > 1
        if len(exp_strands) > 1 and not native_mixed:
            n_skip += 1                      # feature physically split by a rearrangement
            continue
        sym = len(exp_set ^ got_set)
        ok = sym <= a.tol * 2 * len(parts) and (got_strands == exp_strands or len(exp_strands) != 1)
        if ok:
            n_ok += 1
        else:
            n_bad += 1
            bad.append((rf.type, key, "exp %d got %d symdiff %d" % (len(exp_set), len(got_set), sym),
                        sorted(exp_strands), sorted(got_strands), str(of.location)))
    print("location check: %d correct, %d wrong, %d split-by-rearrangement, %d not transferred"
          % (n_ok, n_bad, n_skip, len(seen_missing)))
    for b in bad[:15]:
        print("  WRONG", b)
    for m in seen_missing[:20]:
        print("  not transferred:", m[0], m[1], m[2], "expected:", m[3])
    # summary lines
    for line in summary.splitlines():
        if any(k in line for k in ("transferred successfully", "failed  ", "identical", "with differences",
                                   "disruptive", "orthologous blocks", "SNPs", "partial")):
            print("  |", line.strip())


if __name__ == "__main__":
    main()
