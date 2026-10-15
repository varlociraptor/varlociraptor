#!/usr/bin/env bash
# Generate the realignment benchmark data with mason (SeqAn): a 50 kb random reference, two
# haplotypes carrying SNVs and indels of 1-12 bp, and 100x paired-end 150 bp Illumina reads
# simulated from them and aligned with bwa. The candidates are the simulated variants, split
# into biallelic records. Deterministic (seed 7). All tools are in the pixi environment:
#
#     pixi run bash tests/resources/benchmark_realignment/generate.sh
set -euo pipefail
cd "$(dirname "$0")"
SEED=7
mason_genome -q -s $SEED -l 50000 -o ref.fa
samtools faidx ref.fa
mason_variator -q -s $SEED -ir ref.fa -n 2 \
    --snp-rate 0.002 --small-indel-rate 0.0012 --min-small-indel-size 1 --max-small-indel-size 12 \
    --sv-indel-rate 0 --sv-inversion-rate 0 --sv-translocation-rate 0 --sv-duplication-rate 0 \
    -ov variants.vcf
mason_simulator -q --seed $SEED -ir ref.fa -iv variants.vcf -n 16667 \
    --illumina-read-length 150 --fragment-mean-size 350 --fragment-size-std-dev 40 \
    -o reads_1.fq -or reads_2.fq
bwa index ref.fa 2> /dev/null
bwa mem -R '@RG\tID:s\tSM:s' ref.fa reads_1.fq reads_2.fq 2> /dev/null | samtools sort -o alignment.bam -
samtools index alignment.bam
bcftools view -G variants.vcf | bcftools annotate -x INFO | bcftools norm -m -any -o candidates.vcf
rm -f variants.vcf reads_1.fq reads_2.fq ref.fa.amb ref.fa.ann ref.fa.bwt ref.fa.pac ref.fa.sa
