# ReadME

This is the repo containing all the code used to predict mutability for the plant diversity project available at: https://github.com/eliciapauls/plant_diversity/tree/main. 

This pipleline takes the genome fastas (`data/genomes/`) and gffs (`data/gffs/raw/`), applies the mutabiltiy model generated from a _Arabidopsis thaliana_ mutation accumulation experiment, and creates mutability values for every gene provided in the species' gff (example output in `data/mut_tables/A.halleri/A.halleri_geneMutAv_1000.txt`). 

This pipeline is managed through SnakeMake. The rulegrpah.png and filegraph.png display the sequence steps for the analysis scripts in `analysis/modules/`. 

