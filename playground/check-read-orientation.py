import pysam
import sys
from collections import defaultdict

orientations = defaultdict(int)

with pysam.AlignmentFile(sys.argv[1]) as bam:
    for i, rec in enumerate(bam):
        if rec.is_paired and rec.is_proper_pair and rec.reference_id == rec.next_reference_id:
            strand = "R" if rec.is_reverse else "F"
            mate_strand = "R" if rec.mate_is_reverse else "F"
            read = "1"
            mate_read = "2"
            if not rec.is_read1:
                read, mate_read = mate_read, read

            orientation = [strand + read, mate_strand + mate_read]

            if rec.reference_start > rec.next_reference_start:
                orientation = orientation[::-1]

            orientations["".join(orientation)] += 1

        if i == 100000:
            break

print(orientations)
