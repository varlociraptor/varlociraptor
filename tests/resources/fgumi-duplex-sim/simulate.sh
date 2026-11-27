#!/usr/bin/env bash
# Regenerate the data of the test_fgumi_duplex_* testcases from simulated duplex reads.
# Usage: FGUMI=<fgumi 0.9.0 built with --features simulate> pixi run simulate-duplex-testcases
# An existing testcase.yaml is kept, since its expected values and comments are written by hand;
# delete it to get a fresh one.
set -euo pipefail

sim=$(cd "$(dirname "$0")" && pwd)
testcases=$sim/../testcases
varlociraptor=${VARLOCIRAPTOR:-$sim/../../../target/debug/varlociraptor}
fgumi=${FGUMI:?set FGUMI to an fgumi binary built with --features simulate}
"$fgumi" --version | grep -qx "fgumi 0.9.0"

cd "$(mktemp -d)"

# random 2 kb reference with a G at the variant site
python -c '
import random
rng = random.Random(1)
seq = [rng.choice("ACGT") for _ in range(2000)]
seq[999] = "G"
print(">chr1")
print("".join(seq))
' > ref.fa
samtools faidx ref.fa
samtools dict ref.fa -o ref.dict
bwa index ref.fa 2> /dev/null
printf '##fileformat=VCFv4.2\n##contig=<ID=chr1,length=2000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t1000\t.\tG\tT\t.\t.\t.\n' > candidates.vcf

"$fgumi" simulate grouped-reads --duplex --reference ref.fa --seed 1 \
    --num-molecules 1500 --num-positions 500 --umis-per-position 3 \
    --error-rate 0.001 --truth truth.tsv --output grouped.bam

# standard (non-UMI) normal sample: independent molecules, one read pair each
"$fgumi" simulate grouped-reads --reference ref.fa --seed 2 \
    --num-molecules 1500 --num-positions 500 --umis-per-position 3 \
    --error-rate 0.001 --truth normal.truth.tsv --output normal.grouped.bam
samtools view -b -e 'qname =~ "_read0000$"' normal.grouped.bam | samtools sort -o normal.bam
samtools index normal.bam
"$varlociraptor" preprocess variants ref.fa --candidates candidates.vcf \
    --bam normal.bam --output normal.bcf

# install_testcase SRC DST: the data files always, testcase.yaml only when DST has none
install_testcase() {
    local src=$1 dst=$2
    mkdir -p "$dst"
    cp "$src"/*.bam "$src"/*.vcf "$src"/*.fa "$src"/scenario.yaml "$dst/"
    [ -e "$dst/testcase.yaml" ] || cp "$src/testcase.yaml" "$dst/"
}

# testcase NAME MIN_READS FILTER_ARGS CALL_ARGS INJECT_ARGS...
testcase() {
    local name=$1 min_reads=$2 filter_args=$3 call_args=$4
    shift 4
    python "$sim/inject.py" grouped.bam "$name.grouped.bam" --pos 1000 --alt T "$@"
    "$fgumi" duplex --input "$name.grouped.bam" --output "$name.consensus.bam" \
        --min-reads "$min_reads"
    "$fgumi" fastq --input "$name.consensus.bam" \
        | bwa mem -p -Y ref.fa - 2> /dev/null \
        | "$fgumi" zipper --unmapped "$name.consensus.bam" --reference ref.fa \
            --tags-to-reverse Consensus --tags-to-revcomp Consensus --output "$name.mapped.bam"
    "$fgumi" filter --input "$name.mapped.bam" --output "$name.filtered.bam" --ref ref.fa \
        --min-reads "$min_reads" $filter_args
    samtools sort -o "$name.bam" "$name.filtered.bam"
    samtools index "$name.bam"
    "$varlociraptor" preprocess variants ref.fa --candidates candidates.vcf \
        --bam "$name.bam" --output "$name.bcf"
    "$varlociraptor" call variants --testcase-locus chr1:1000 --testcase-prefix "$name.testcase" \
        $call_args generic --scenario "$sim/scenario.yaml" --obs "duplex=$name.bcf" normal=normal.bcf
    install_testcase "$name.testcase" "$testcases/$name"
}

# true variant, all consensus records in F1R2 orientation
testcase test_fgumi_duplex_orientation_skew 1,1,1 "" "" --op duplex:6:F1R2
# true variant in a single duplex molecule
testcase test_fgumi_duplex_singleton 1,1,1 "" "" --op duplex:1:any
# single-strand damage, masked by the consensus filter
testcase test_fgumi_duplex_single_strand_damage 1,1,1 --require-single-strand-agreement "" \
    --op a:6:any
# single-strand damage in F2R1 families, winning the strand disagreement in the consensus
# (strand bias omitted, to test read orientation bias alone)
testcase test_fgumi_duplex_strand_disagreement 1,1,1 "" --omit-strand-bias --other-qual 12 \
    --op a:5:F2R1 --op b:5:F2R1

# the strand disagreement without the normal: the duplex sample alone has to show the orientation artifact
tumor_only() {
    local src=$testcases/test_fgumi_duplex_strand_disagreement dst=$testcases/test_fgumi_duplex_strand_disagreement_tumor_only
    mkdir -p "$dst"
    cp "$src"/{candidates.vcf,duplex.bam,ref.fa} "$dst/"
    cp "$sim/scenario_duplex_only.yaml" "$dst/scenario.yaml"
    [ -e "$dst/testcase.yaml" ] || python - "$src/testcase.yaml" > "$dst/testcase.yaml" <<'PY'
import re, sys
# drop the normal sample: its entry is indented by two spaces, its fields by four
print(re.sub(r"\n  normal:\n(    .*\n)+", "\n", open(sys.argv[1]).read()), end="")
PY
}
tumor_only
