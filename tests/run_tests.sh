#!/bin/bash
# Regression tests: simulate targets with known edits from real reference
# genomes, run transfer_annotation.py and check every transferred feature
# location against the simulated truth.
#   bash tests/run_tests.sh            (takes ~30 s)
HERE="$(cd "$(dirname "$0")" && pwd)"
V="$HERE/simulate_and_verify.py"
PX="$HERE/data/NC_001422.gbk"   # phiX174, 5.4 kb circular, overlapping + origin-spanning genes
CP="$HERE/data/NC_000932.gb"    # Arabidopsis chloroplast, 154 kb, inverted repeats, introns
YP="$HERE/data/NC_005816.gb"    # Yersinia pestis plasmid pPCP1, 9.6 kb
OUT="${1:-$HERE/test_output}"
mkdir -p "$OUT"; cd "$OUT"
run() {
  local name=$1 ref=$2 ops=$3; shift 3
  local out
  out=$(python3 "$V" "$ref" "r_$name" --ops "$ops" "$@" 2>&1)
  printf "%-22s %s\n" "$name" "$(echo "$out" | grep -o 'location check:.*')"
  echo "$out" | grep -qE "Traceback|RUN FAILED" && echo "$out" | tail -20
  return 0
}
echo "Expected: 0 wrong, except where noted"
run phix_self          "$PX" '[]'
run phix_rotated       "$PX" '[["snp",0.01],["indel",5,4],["rotate",1000]]'
run phix_revcomp       "$PX" '[["snp",0.02],["rotate",2500],["revcomp"]]'
run phix_30pct_div     "$PX" '[["snp",0.30],["indel",10,3]]'
run cp_self            "$CP" '[]'
run cp_rearranged      "$CP" '[["snp",0.02],["indel",30,6],["delete",380,1450],["insert",50000,1500],["invert",60000,70000],["rotate",20000]]'
echo "   (cp_rearranged: psbA is deleted -> not transferred; genes cut by the inversion are partial)"
run cp_6_contigs       "$CP" '[["snp",0.01],["indel",10,3]]' --split 6 --extra "--topology linear"
echo "   (cp_6_contigs: 'wrong' = genes cut by contig ends, transferred as partial)"
run cp_30pct_div       "$CP" '[["snp",0.30],["indel",30,3]]'
run plasmid_20pct_div  "$YP" '[["snp",0.20],["indel",10,3]]'
