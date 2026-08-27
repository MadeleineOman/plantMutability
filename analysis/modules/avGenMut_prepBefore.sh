#!/bin/bash
#sed 's/"//g' ../data/data_dataframes/ProbabilityDf{wildcards.species}.txt | tail -n +2 >             data/data_dataframes/Probability{wildcards.species}_prepped.txt 

species="$1"
replicate="$2"

sed -i 's/,/\t/g' data/data_dataframes/"$species"/ProbabilityDf"$species"_"$replicate".csv #csv to tsv
sed -i 's/"//g' data/data_dataframes/"$species"/ProbabilityDf"$species"_"$replicate".csv #removing the quotes that have been fucking everything 
tail -n +2 data/data_dataframes/"$species"/ProbabilityDf"$species"_"$replicate".csv |awk '{print $0"\t"($NF+1)}' > data/data_dataframes/"$species"/ProbabilityDf"$species"_prepped_"$replicate".txt

sort -k38,38 -k39,39n  data/data_dataframes/"$species"/ProbabilityDf"$species"_prepped_"$replicate".txt > data/data_dataframes/"$species"/ProbabilityDf"$species"_prepped_sorted_"$replicate".txt 

rm data/data_dataframes/"$species"/ProbabilityDf"$species"_prepped_"$replicate".txt
rm data/data_dataframes/"$species"/ProbabilityDf"$species"_"$replicate".csv