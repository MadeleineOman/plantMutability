library(dplyr)
library(stringr)
library(stringi)
library(corrplot)
library(rlang)
library(ggplot2)



equiv_toLowest=FALSE
exclude_CpG=FALSE
exclude_TCX_CCX=FALSE
exclude_triplet=FALSE
NA_omit=TRUE
null_model =FALSE


args = commandArgs(trailingOnly=TRUE)
species = args[1]
replicate = args[2]
# species = "A.thaliana"
# replicate = "1000"
do_selfDiag = FALSE
doThatRobThang = FALSE
isModel=FALSE


options(warn=1)

model_desc_modify = ""
tmp_file_path=""
start_time = Sys.time()

#import the data 
if(isModel==TRUE){
    input_filePath = paste(tmp_file_path,"data/data_dataframes/",species,"/",species,"_predictorDf_model.txt",sep="")
}else{input_filePath = paste(tmp_file_path,"data/data_dataframes/",species,"/",species,"_predictorDf_",replicate,".txt",sep="")}
all_data <- read.table(input_filePath, header = FALSE,sep="\t")



if (doThatRobThang==FALSE){
    if (isModel==TRUE){
        colnames(all_data) <- c("Chromosome","site","triplet","mutation_status",
                           "Apercent.0","Gpercent.0","Cpercent.0","Tpercent.0",
                           "Apercent.100","Gpercent.100","Cpercent.100","Tpercent.100",
                           "Apercent.10000","Gpercent.10000","Cpercent.10000","Tpercent.10000")
    }else{
        colnames(all_data) <- c("Chromosome","site","triplet",
                           "Apercent.0","Gpercent.0","Cpercent.0","Tpercent.0",
                           "Apercent.100","Gpercent.100","Cpercent.100","Tpercent.100",
                           "Apercent.10000","Gpercent.10000","Cpercent.10000","Tpercent.10000")
    }
    
}else{
if (isModel==FALSE){
    colnames(all_data) <- c("Chromosome","site","base","seqId","cds_pos","strand","frame","codon","aminoacid","degen","triplet",
                           "Apercent.0","Gpercent.0","Cpercent.0","Tpercent.0",
                           "Apercent.100","Gpercent.100","Cpercent.100","Tpercent.100",
                           "Apercent.10000","Gpercent.10000","Cpercent.10000","Tpercent.10000")
}else{
    colnames(all_data) <- c("Chromosome","site","triplet","mutation_status",
                           "Apercent.0","Gpercent.0","Cpercent.0","Tpercent.0",
                           "Apercent.100","Gpercent.100","Cpercent.100","Tpercent.100",
                           "Apercent.10000","Gpercent.10000","Cpercent.10000","Tpercent.10000",
                            "base","seqId", # why arent these inthe a.thaliana table 
                            "cds_pos","strand","frame","codon","aminoacid","degen")
}}






# if (isModel==TRUE){ 
#     mut_file_path = paste(tmp_file_path,"data/mutations/Unique_mutations_SNMs_final.txt",sep="")
#     mutations_df <- read.table(mut_file_path,header=TRUE,col.names = c('Chromosome', 'site', 'MA_Line', 'Type', 'ref', 'alt', 'Genomic_region',
#     'TAIR_ID'))


#     #merging on the shared columns (chromosome and site) 
#     mutations_df$site <- mutations_df$site-1
#     all_data <- merge(x = all_data, y = mutations_df, all.x = TRUE)

#     #MAKING SURE THE reference base matches the middle of the triplet (all tissues) 
#     all_data$triplet_middle <- str_sub(all_data$triplet,2,2) 
#     all_data_muts = filter(all_data,all_data$mutation_status==1)
#     all_data_muts = filter(all_data_muts,nchar(all_data_muts$ref)==1)#removing the rows where the "ref" column has more than 1 element (deletion) 
#     stopifnot(all_data_muts$triplet_middle == all_data_muts$ref) #another test for alignment to reference. 
#     all_data <- all_data[,!(names(all_data) %in% c('ref','alt',"triplet_middle"))]

#     #removing the mutation-df associated columns as the non-mutant sites will have NAs 

#     all_data <- all_data[,!(names(all_data) %in% c('notes','Genomic_region','TAIR_ID','MA_Line', 'Type'))]
# }






# print("where")
# print(dim(all_data))
# if (doThatRobThang==TRUE){
#     model_desc_modify = paste(model_desc_modify,"_doThatRobThang",sep="")}
# if (exclude_CpG==TRUE){
#     all_data <-all_data[!str_detect(all_data$triplet,"CG"),]
#     model_desc_modify = paste(model_desc_modify,"_noCpG",sep="")}
# if (exclude_TCX_CCX==TRUE){
#     model_desc_modify = paste(model_desc_modify,"_noTCX_CCX",sep="")}

# print(summary(as.factor(all_data$triplet)))
# if (exclude_triplet==TRUE){
#     all_data <- all_data[,!(names(all_data) %in% c('triplet'))]
#     model_desc_modify = paste(model_desc_modify,"_noTriplets",sep="")}
# # print(colnames(all_data))

#setting up the error output file 
error_output_file = paste(tmp_file_path,"data/model/",species,"/create_model_text_output",model_desc_modify,"_",replicate,".txt",sep="")
cat("output for the create_model notebook",file=error_output_file,sep="\n")


print("exactly")
print(dim(all_data))
#handling gc content"
all_data <- all_data %>%
    mutate(GC_content.0 = Gpercent.0+Cpercent.0) %>%
    mutate(GC_content.100 = Gpercent.100+Cpercent.100) %>% 
    mutate(GC_content.10000 = Gpercent.10000+Cpercent.10000)
all_data <- all_data[,!(names(all_data) %in% c('Apercent.0','Gpercent.0','Cpercent.0','Tpercent.0',
                                       'Apercent.100','Gpercent.100','Cpercent.100','Tpercent.100',
                                       'Apercent.10000','Gpercent.10000','Cpercent.10000','Tpercent.10000'))]




#editing the triplets 
# if (exclude_triplet==FALSE){
all_data$triplet <- toupper(all_data$triplet)
#filter for rows that dont have NNN as the triplet --> write the details to file 

cat(paste(sum(str_detect(all_data$triplet,"[N,R,Y,M,K,S,W,H,B,V,D]")), "rows removed due to N in triplet, ",sep=" "))
all_data <- all_data[(str_detect(all_data$triplet,"[N,R,Y,M,K,S,W,H,B,V,D]",negate=TRUE)),]
cat(paste(nrow(all_data),"rows left",sep=" "),file=error_output_file,sep="\n",append=TRUE)
# }





print("is")
print(dim(all_data))
#converting 64-->32 triplets 

rc_removeAG <- function(dna){
    middle_base = substr(dna, 2, 2)
    if(middle_base %in% c("A","G")){
        dna <- stri_reverse(chartr("acgtACGT", "tgcaTGCA", dna))}
    return(dna)
}#substring slicing https://www.johnmyleswhite.com/notebook/2009/02/25/text-processing-in-r/
all_data$triplet <- unlist(lapply(as.character(all_data$triplet),rc_removeAG)) #need unlistto turn the list into a vector 



#check for problems with levels / column tye etc. 

all_data$triplet <- as.character(all_data$triplet) #need to convert to char so the as.factor properly reducs to 64 levels after the removal of 'triplets" longer than 3
string_to_print = paste(nrow(all_data[nchar(as.character(all_data$triplet))!=3,]),"rows removed due to triplet larger than 3 in length ",sep=" ")
cat(string_to_print,file=error_output_file,sep="\n",append=TRUE)
all_data <- all_data[nchar(all_data$triplet)==3,] #make sure only including triplets 
string_to_print=paste(nrow(all_data)," rows left", sept = " ")
cat(string_to_print,file=error_output_file,sep="\n",append=TRUE)


print("the")
print(dim(all_data))

#factorizing the columns 
all_data$triplet <- as.factor(as.character(all_data$triplet))
all_data$Chromosome <- as.factor(all_data$Chromosome)

print("problem")
print(dim(all_data))

#checking that the factor variables are correct (mutation status, triplets, chroms) --> will raise error if not true 
stopifnot(length(levels(all_data$triplet))==32)

print("linda")
print(dim(all_data))


#na omit and checking how many rows were removed 
string_to_print = paste(nrow(all_data)- nrow(na.omit(all_data)), "lost to NA values, ", nrow(na.omit(all_data))," rows remain after",sep=" ")
cat(string_to_print,file=error_output_file,sep="\n",append=TRUE)
if (NA_omit ==TRUE){
    all_data <- na.omit(all_data)
    }

#randomize the mutation if null 
if (null_model==TRUE){
    model_desc_modify = paste(model_desc_modify,"_null_model",sep="")
    all_data$mutation_status <- sample(x =c(0,1), replace=TRUE,size=nrow(all_data))}

print(".... like")
print(dim(all_data))


#printing the total rows included to file 
string_to_print = paste("total nrow", nrow(all_data),sep=" ")
cat(string_to_print,file=error_output_file,sep="\n",append=TRUE)

print ("just get")
print(dim(all_data))


cnames_toInclude = c("triplet","GC_content.0","GC_content.100","GC_content.10000")
cnames_toInclude = append(cnames_toInclude, c("Chromosome","site"))

all_data <- all_data[,(names(all_data) %in% cnames_toInclude)]



print("the")
print(dim(all_data))
all_data$dummy <- 1#need a dumym column that the model.matrix can remove (as it has to remove a column aparently)
non_num_preds = c("Chromosome","site")

print("job")
print(dim(all_data))



preds_to_standardize<-colnames(all_data)[!(names(all_data) %in% non_num_preds)]

chrom_col <- all_data$Chromosome
site_col <- all_data$site
all_data <-(data.frame((model.matrix(dummy~., all_data[, preds_to_standardize]) )))
all_data <- cbind(all_data,site_col)



print("done")
print(dim(all_data))

#print("standardizing")
#STANDARDIZING~~~~~~~~~~~~~~~~~~~~~~~~~~~~
non_num_preds = c("mutation_status","Chromosome","site","chrom_col","site_col")
preds_to_standardize<-!(names(all_data) %in% non_num_preds)#yes i have to do this again because i chnaged the order of the mutant column

 
all_data[, preds_to_standardize]<- (as.data.frame(scale(all_data[, preds_to_standardize],center=TRUE,scale=TRUE)))



print("already")
print(dim(all_data))
                          


#saving the model variables 
filename =paste(tmp_file_path,"data/data_dataframes/",species,"/all_data_readyForPrediction",model_desc_modify,"_",species,"_",replicate,"_bareMin.csv", sep="")


write.csv(all_data,filename)    
print(paste("whole data prep took",Sys.time()-start_time,sep=" "))