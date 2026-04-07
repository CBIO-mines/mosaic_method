#!/usr/bin/env sh

# $1 is the name of the genome to mask
# $2 is the genome directory
# $3 is the output directory

cat   "$2/$1" \
    | ./tools/fasta_fragments.py --fragment=200 --step=100 \
    | lastz "$2/$1"[multiple,unmask,nameparse=darkspace] /dev/stdin --masking=6 \
        --progress+masking=10K \
        --format=none --outputmasking+:soft="$3/$1.masked_intervals.dat" \
        --notransition


cat  "$2/$1" \
    | ./tools/fasta_softmask_intervals.py --origin=1 "$3/$1.masked_intervals.dat" \
    > "$3/$1"
