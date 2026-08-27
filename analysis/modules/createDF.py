# imports ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import random
import pandas as pd
import numpy as np 
from numpy.random import choice
import collections
from Bio import AlignIO
from Bio import SeqIO

import pysam 
from datetime import datetime
import gzip
import multiprocessing
import sys 
import json 
import time # for timing the loop 

#home-made modules : 
# sys.path.append('/research/projects/hsapiens/mutability/analysis/global/track_data/annotation/') 
# import annotation_handling
# sys.path.append('/research/projects/hsapiens/mutability/analysis/modules/closest_value/') 
sys.path.append('/node1nfs2/archive/research/projects/hsapiens/mutability/mutability/analysis/modules/closest_value/')
import closest_val
# sys.path.append('/research/projects/hsapiens/mutability/analysis/modules/closest_tss/')
# import dist2tss

tmp_file_dir =""
# species = "A.halleri"
# replicate = "1000"
replicate = (sys.argv[2])
species = sys.argv[1]

list_of_surrounding_contexts = [0,100,10000]






# all_gene_lines = open(tmp_file_dir+"data/gffs/processed/{s}/{s}_allGenes_sorted_{r}.gff".format(s=species,r=replicate)).readlines()
all_geneProm_lines = open(tmp_file_dir+"data/gffs/promoters/processed/{s}/{s}Prom_{r}.bed".format(s=species,r=replicate)).readlines()


sites = []
for line in all_geneProm_lines: 
    chrom= line.split()[0]
    start= int(line.split()[1])
    end= int(line.split()[2])
    for i in range(start,end):
        sites.append([chrom,i])


chrom_labels = []
for line in all_geneProm_lines: 
    chrom=line.split("\t")[0]
    if chrom not in chrom_labels: 
        chrom_labels.append(chrom)
        

        
        
# genral declarations before the big function ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
distance_max = max(list_of_surrounding_contexts)




fastas_dict = {}   # creating dictionary with fasta alignment, length of seq, 
print("making the fastas dictionary")
with gzip.open(tmp_file_dir+"data/genomes/{s}.fasta.gz".format(s=species), "rt") as handle:
    for seq in SeqIO.parse(handle, "fasta"):
        if seq.id in chrom_labels: 
            fastas_dict[seq.id] = [seq.seq,len(str(seq.seq))]
print(len(fastas_dict))




# filename_df = tmp_file_dir+'data/data_dataframes/{s}/{s}_predictorDf_{r}.txt'.format(s=species,r=replicate)
# with open(filename_df, 'a') as f:
# for line in sites: 
scaffolds_notin_fasta = []
def predictor_rowString(line): 
#     line = sites[0]
    row = []
    chrom,pos=line
    row.extend([chrom,pos])
    try: fastas_dict[chrom] 
    except: 
        scaffolds_notin_fasta.append(chrom)
    else: 
        alignment = fastas_dict[chrom] 
        triplet = str(alignment[0][pos-1:pos+2])
        row.append(triplet)
        for distance in list_of_surrounding_contexts:         
            if triplet != '': 
                min_pos = max([0,pos-distance])
                seq_around = str(alignment[0][min_pos:pos+distance+1])
                Acount = seq_around.count('a')+seq_around.count("A")
                Gcount = seq_around.count('g')+seq_around.count("G")
                Ccount = seq_around.count('c')+seq_around.count("C")
                Tcount = seq_around.count('t')+seq_around.count("T")
                if len(seq_around) != 0:
                    Apercent = Acount/len(seq_around)
                    Gpercent = Gcount/len(seq_around)
                    Cpercent = Ccount/len(seq_around)
                    Tpercent = Tcount/len(seq_around)
                    row.extend([Apercent, Gpercent, Cpercent, Tpercent])
                else: 
                    row.extend(['NA','NA','NA','NA'])

            else: 
                row.extend(['NA','NA','NA','NA'])
                #list_no_seq_at_site.append(site)

        row_string = str()
        for i in range(0,len(row)): 
            row_string = row_string+str(row[i])+"\t"
        row_string = row_string.rstrip("\t") # dont need to add the "\n" here as it is added below int he f.write 
        return(row_string)


#WRITE THE MODEL #
filename_df = tmp_file_dir+'data/data_dataframes/{s}/{s}_predictorDf_{r}.txt'.format(s=species,r=replicate)
start_time = time.time()
print("starting big loop")
timestamp = datetime.now().strftime("%Y_%m_%d")
error_log = str("df created on "+timestamp+"\n")


def rowString_handler():
    p = multiprocessing.Pool(10)
    with open(filename_df, 'a') as f:
        for result in p.imap(predictor_rowString, sites):
            f.write('%s\n' % result)

if __name__=='__main__':
    rowString_handler()

error_log += (("creating the df loop took "+str(time.time()-start_time)[0:4]+" seconds\n"))
scaffolds_notin_fasta_uniq = list(set(scaffolds_notin_fasta))
print(len(scaffolds_notin_fasta_uniq)," uniq scaffolds not present in the fasta")
print(scaffolds_notin_fasta_uniq)
error_log+= "scaffolds not in the fasta: "+ str(scaffolds_notin_fasta_uniq)+"\n"

#writing error log to file 
with open(tmp_file_dir+"data/data_dataframes/{s}/{s}_predictorDf_errorlog_{r}.txt".format(s=species,r=replicate),"w") as f: 
      f.write(error_log)




