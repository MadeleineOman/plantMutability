#!bin/bash
species=$1

grep -v '^#' data/gffs/"$species".gff | awk '{if ($3=="gene")print $0}' > data/gffs/"$species"_allGenes.gff
sort -k1,1 -k4,4n data/gffs/"$species"_allGenes.gff > data/gffs/"$species"_allGenes_sorted.gff