#!/bin/bash

igv.sh --genome assembly/assembly.fasta --batch /dev/stdin <<EOF
load short_reads/short_reads_to_assembly.cram
colorBy READ_STRAND
load long_reads/long_reads_to_assembly.cram
colorBy READ_STRAND
#load ./pipeline/ref/hsv1.gff
#expand hsv1.gff
#snapshot out/snapshot.png
EOF

