#!/usr/bin/env python3
"""Generate the realignment benchmark data: a 50 kb random reference carrying ~175 heterozygous
variants (SNVs, indels of 1-12 bp and clustered SNVs), and 100x paired-end 150 bp reads with
0.2 % substitution errors and varying base qualities. Deterministic (seed 7).

Requires minimap2, samtools and freebayes:

    python3 generate.py && \
    minimap2 -ax sr -R '@RG\tID:s\tSM:sample' ref.fa reads_1.fq reads_2.fq | samtools sort -o alignment.bam - && \
    samtools index alignment.bam && samtools faidx ref.fa && \
    freebayes -f ref.fa --pooled-continuous --min-alternate-count 2 --min-alternate-fraction 0.05 alignment.bam > candidates.vcf

The reads and the aligner are only a means to obtain a realistic BAM; the benchmark measures
`varlociraptor preprocess variants` on the candidates.
"""
import random

random.seed(7)
L, READ, DEPTH, INSERT = 50_000, 150, 100, 350
ref = "".join(random.choice("ACGT") for _ in range(L))

edits = []
pos = 200
while pos < L - 200:
    kind = random.random()
    if kind < 0.5:
        edits.append((pos, 1, random.choice([b for b in "ACGT" if b != ref[pos - 1]])))
    elif kind < 0.7:
        edits.append((pos, random.randint(1, 12), ""))
    elif kind < 0.85:
        edits.append((pos, 0, "".join(random.choice("ACGT") for _ in range(random.randint(1, 12)))))
    else:
        d = random.randint(1, 4)
        edits.append((pos, 1, random.choice([b for b in "ACGT" if b != ref[pos - 1]])))
        edits.append((pos + d, 1, random.choice([b for b in "ACGT" if b != ref[pos + d - 1]])))
    pos += random.randint(150, 400)
edits.sort()
out, cur = [], 0
for p, rlen, alt in edits:
    if p - 1 < cur:
        continue
    out.append(ref[cur:p - 1])
    out.append(alt)
    cur = p - 1 + rlen
out.append(ref[cur:])
alt_seq = "".join(out)

COMP = str.maketrans("ACGT", "TGCA")

def read(seq, start):
    r = seq[start:start + READ]
    r = "".join(random.choice("ACGT") if random.random() < 0.002 else b for b in r)
    q = "".join(chr(33 + random.choice([37, 37, 37, 32, 25, 12])) for _ in range(READ))
    return r, q

with open("ref.fa", "w") as f:
    f.write(f">ref\n{ref}\n")
n = (L * DEPTH // READ) // 4  # pairs per haplotype
with open("reads_1.fq", "w") as f1, open("reads_2.fq", "w") as f2:
    for tag, seq in (("ref", ref), ("alt", alt_seq)):
        for k in range(n):
            ins = max(READ + 20, int(random.gauss(INSERT, 40)))
            s = random.randrange(0, len(seq) - ins + 1)
            r1, q1 = read(seq, s)
            r2, q2 = read(seq, s + ins - READ)
            r2 = r2.translate(COMP)[::-1]
            f1.write(f"@{tag}_{k}/1\n{r1}\n+\n{q1}\n")
            f2.write(f"@{tag}_{k}/2\n{r2}\n+\n{q2[::-1]}\n")
print(f"{len(edits)} variants, {2 * n} read pairs")
