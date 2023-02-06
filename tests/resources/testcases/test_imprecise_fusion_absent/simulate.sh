#!/bin/sh

mason_simulator -ir ref.fa -o left.fq -or right.fq --illumina-read-length 20 -n 100 --fragment-mean-size 70 --fragment-min-size 70 --fragment-max-size 70
bwa index ref.fa
bwa mem -T 10 -k 10 ref.fa left.fq right.fq | samtools sort -O BAM - > mapped.bam
