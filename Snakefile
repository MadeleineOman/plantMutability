rule all: 
    input: 
        "data/data_dataframes/A.thaliana/ProbabilityDfA.thaliana_prepped_sorted_16000.txt.gz.tbi"

        

wildcard_constraints:
    replicate="[0-9]+"    
    #species="[A-Z]\.[A-Za-z]+"
    
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rule prepGFFforCreateDF:
    input: 
        "data/gffs/{species}.gff"
    output: 
        "/data/gffs/{species}_allGenes_sorted.gff"
    threads: 10
    conda: "conda_envs/conda_createDF.yml"
    shell: 
        "bash analysis/modules/prepGffForCreateDf.sh {wildcards.species}"

rule createDF: 
    input: 
        "data/genomes/{species}.fasta.gz", 
        "data/gffs/processed/{species}/{species}_allGenes_sorted_{replicate}.gff"
    output: "data/data_dataframes/{species}/{species}_predictorDf_{replicate}.txt"
    threads: 5
    conda: "conda_envs/conda_createDF.yml"
    shell: 
        "python3 analysis/modules/createDF.py {wildcards.species} {wildcards.replicate}"

        
rule dataPrep: 
    input: 
        "data/data_dataframes/{species}/{species}_predictorDf_{replicate}.txt"
    output: "data/data_dataframes/{species}/all_data_readyForPrediction_{species}_{replicate}.csv"
    threads: 2
    conda: "conda_envs/conda_RcreateDfModel_env.yml"
    shell: "Rscript --vanilla analysis/modules/dataPrepWrangle.R {wildcards.species} {wildcards.replicate} ;"
            "rm data/data_dataframes/{wildcards.species}/{wildcards.species}_predictorDf_{wildcards.replicate}.txt"
    
    
rule predict: 
    input: 
        "data/data_dataframes/{species}/all_data_readyForPrediction_{species}_{replicate}.csv"
    output: "data/data_dataframes/{species}/ProbabilityDf{species}_{replicate}.csv"
    
    threads: 10
    conda: "conda_envs/conda_RcreateDfModel_env.yml"
    shell: "Rscript --vanilla analysis/modules/predict.R {wildcards.species} {wildcards.replicate} ;"
            "rm data/data_dataframes/{wildcards.species}/all_data_readyForPrediction_{wildcards.species}_{wildcards.replicate}.csv"
  
rule prepFor_avMutgene: 
    input: 
         "data/data_dataframes/{species}/ProbabilityDf{species}_{replicate}.csv"
    output:  "data/data_dataframes/{species}/ProbabilityDf{species}_prepped_sorted_{replicate}.txt"
    threads: 1
    conda: "conda_envs/conda_createDF.yml"
    shell: 
        "bash analysis/modules/avGenMut_prepBefore.sh {wildcards.species} {wildcards.replicate}"    
        
        

rule prepFor_avMutgene_bzip: 
    input: 
         "data/data_dataframes/{species}/ProbabilityDf{species}_prepped_sorted_{replicate}.txt"
    output:  "data/data_dataframes/{species}/ProbabilityDf{species}_prepped_sorted_{replicate}.txt.gz"
    threads: 1
    conda: "conda_envs/conda_snakeSomMut_env.yml"
    shell: 
        "bgzip -f data/data_dataframes/{wildcards.species}/ProbabilityDf{wildcards.species}_prepped_sorted_{wildcards.replicate}.txt "  
  
  

rule prepFor_avMutgene_Tabix: 
    input: 
         "data/data_dataframes/{species}/ProbabilityDf{species}_prepped_sorted_{replicate}.txt.gz"
    output:  "data/data_dataframes/{species}/ProbabilityDf{species}_prepped_sorted_{replicate}.txt.gz.tbi"
    threads: 1
    conda: "conda_envs/conda_snakeSomMut_env.yml"
    shell: 
        "tabix -s 38 -b 39 -e 40  data/data_dataframes/{wildcards.species}/ProbabilityDf{wildcards.species}_prepped_sorted_{wildcards.replicate}.txt.gz"    
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
rule confirm_avMutgene: 
    input: 
        "data/data_dataframes/{species}/ProbabilityDf{species}_1000.csv",
         "data/data_dataframes/{species}/ProbabilityDf{species}_2000.csv",
         "data/data_dataframes/{species}/ProbabilityDf{species}_3000.csv",
         "data/data_dataframes/{species}/ProbabilityDf{species}_4000.csv",
         "data/data_dataframes/{species}/ProbabilityDf{species}_5000.csv",
         "data/data_dataframes/{species}/ProbabilityDf{species}_6000.csv",
         "data/data_dataframes/{species}/ProbabilityDf{species}_7000.csv",
         "data/data_dataframes/{species}/ProbabilityDf{species}_8000.csv",
         "data/data_dataframes/{species}/ProbabilityDf{species}_9000.csv",
         "data/data_dataframes/{species}/ProbabilityDf{species}_10000.csv",
         "data/data_dataframes/{species}/ProbabilityDf{species}_11000.csv",
         "data/data_dataframes/{species}/ProbabilityDf{species}_12000.csv",
         "data/data_dataframes/{species}/ProbabilityDf{species}_13000.csv",
         "data/data_dataframes/{species}/ProbabilityDf{species}_14000.csv",
         "data/data_dataframes/{species}/ProbabilityDf{species}_15000.csv",
         "data/data_dataframes/{species}/ProbabilityDf{species}_16000.csv",
         "data/data_dataframes/{species}/ProbabilityDf{species}_17000.csv",
         "data/data_dataframes/{species}/ProbabilityDf{species}_18000.csv",
         "data/data_dataframes/{species}/ProbabilityDf{species}_19000.csv",
         "data/data_dataframes/{species}/ProbabilityDf{species}_20000.csv",
         "data/data_dataframes/{species}/ProbabilityDf{species}_21000.csv",
         "data/data_dataframes/{species}/ProbabilityDf{species}_22000.csv",
         "data/data_dataframes/{species}/ProbabilityDf{species}_23000.csv",
         "data/data_dataframes/{species}/ProbabilityDf{species}_24000.csv",
         "data/data_dataframes/{species}/ProbabilityDf{species}_25000.csv",
         "data/data_dataframes/{species}/ProbabilityDf{species}_26000.csv",
         "data/data_dataframes/{species}/ProbabilityDf{species}_27000.csv",
         "data/data_dataframes/{species}/ProbabilityDf{species}_28000.csv",
         "data/data_dataframes/{species}/ProbabilityDf{species}_29000.csv",
         "data/data_dataframes/{species}/ProbabilityDf{species}_30000.csv",
         "data/data_dataframes/{species}/ProbabilityDf{species}_31000.csv",
         "data/data_dataframes/{species}/ProbabilityDf{species}_32000.csv",
         "data/data_dataframes/{species}/ProbabilityDf{species}_33000.csv",
         "data/data_dataframes/{species}/ProbabilityDf{species}_34000.csv",
         "data/data_dataframes/{species}/ProbabilityDf{species}_35000.csv",
         "data/data_dataframes/{species}/ProbabilityDf{species}_36000.csv",
         "data/data_dataframes/{species}/ProbabilityDf{species}_37000.csv",
         "data/data_dataframes/{species}/ProbabilityDf{species}_38000.csv",
         "data/data_dataframes/{species}/ProbabilityDf{species}_39000.csv",
         "data/data_dataframes/{species}/ProbabilityDf{species}_40000.csv",
         "data/data_dataframes/{species}/ProbabilityDf{species}_41000.csv",
         "data/data_dataframes/{species}/ProbabilityDf{species}_42000.csv"
    output:  "data/data_dataframes/ProbabilityDf{species}_confimations.csv"
    threads: 1
    conda: "conda_envs/conda_createDF.yml"
    shell: 
        "echo 'got em boys' > data/data_dataframes/{wildcards.species}/ProbabilityDf_{wildcards.species}_confimations.csv;"
        "bgzip data/data_dataframes/{species}/ProbabilityDf{species}_*"
        
        
        
        
        

    
rule avMutgene: 
    input: 
         "data/data_dataframes/{species}/ProbabilityDf{species}_prepped_sorted_{replicate}.txt.gz"
    output:  "data/mut_tables/{species}/{species}_geneMutAv_{replicate}.txt"
    threads: 10
    conda: "conda_envs/conda_createDF.yml"
    shell: 
        "python3 analysis/modules/avMutGenes.py {wildcards.species} {wildcards.replicate};"
        "rm data/data_dataframes/{wildcards.species}/ProbabilityDf{wilcards.species}_prepped_sorted_{wildcards.replicate}.txt.gz"
        
