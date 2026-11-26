"""Inject an SNV into strand families of fgumi-simulated grouped duplex reads.

Each --op KIND:N:ORIENTATION mutates N randomly chosen molecules covering --pos:
  duplex both families (true variant)
  a      family A only (single-strand damage)
  b      family B only
ORIENTATION (F1R2, F2R1, any) restricts the pair orientation of the first mutated family.
With --other-qual, the base quality of the other family is lowered, such that the
mutated family wins the disagreement in the duplex consensus.
"""

import argparse
import random
from collections import defaultdict

import pysam

MUTATED = {"duplex": "AB", "a": "A", "b": "B"}
OTHER = {"a": "B", "b": "A"}


def orientation(reads):
    r1 = next(r for r in reads if r.is_read1)
    return "F2R1" if r1.is_reverse else "F1R2"


def qpos(read, pos):
    return next(
        (q for q, r in read.get_aligned_pairs(matches_only=True) if r == pos), None
    )


def set_base(read, pos, base=None, qual=None):
    q = qpos(read, pos)
    if q is not None:
        quals = read.query_qualities
        if base is not None:
            seq = read.query_sequence
            read.query_sequence = seq[:q] + base + seq[q + 1 :]
        if qual is not None:
            quals[q] = qual
        read.query_qualities = quals


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("input")
parser.add_argument("output")
parser.add_argument("--pos", type=int, required=True, help="1-based position")
parser.add_argument("--alt", required=True)
parser.add_argument("--op", action="append", required=True)
parser.add_argument("--other-qual", type=int)
parser.add_argument("--seed", type=int, default=0)
args = parser.parse_args()
pos = args.pos - 1

with pysam.AlignmentFile(args.input) as bam:
    header = bam.header
    reads = list(bam)

families = defaultdict(lambda: defaultdict(list))
for read in reads:
    mol, family = read.get_tag("MI").split("/")
    families[mol][family].append(read)

eligible = sorted(
    mol
    for mol, fams in families.items()
    if set(fams) == {"A", "B"}
    and all(any(qpos(r, pos) is not None for r in fams[f]) for f in "AB")
)
random.Random(args.seed).shuffle(eligible)

for op in args.op:
    kind, n, orient = op.split(":")
    chosen = [
        mol
        for mol in eligible
        if orient == "any" or orientation(families[mol][MUTATED[kind][0]]) == orient
    ][: int(n)]
    assert len(chosen) == int(n), f"not enough eligible molecules for {op}"
    for mol in chosen:
        eligible.remove(mol)
        for family in MUTATED[kind]:
            for read in families[mol][family]:
                set_base(read, pos, base=args.alt)
        if kind in OTHER and args.other_qual is not None:
            for read in families[mol][OTHER[kind]]:
                set_base(read, pos, qual=args.other_qual)
    print(op, " ".join(chosen))

with pysam.AlignmentFile(args.output, "wb", header=header) as out:
    for read in reads:
        out.write(read)
