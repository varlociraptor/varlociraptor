#!/bin/sh

mason_simulator -ir ref_fusion.fa -o left.fq -or right.fq --illumina-read-length 20 -n 100 --fragment-mean-size 70 --fragment-min-size 70 --fragment-max-size 70
bwa mem ref.fa left.fq right.fq | samtools view -b - > mapped.bam
