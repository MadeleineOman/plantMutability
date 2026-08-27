library(dplyr)
library(stringr)


equiv_toLowest=FALSE
exclude_triplet = FALSE
exclude_CpG= FALSE
fullModel = TRUE
null_model = FALSE

args = commandArgs(trailingOnly=TRUE)

# species_predON = "C.rubella"
# replicate = "1000"
do_selfDiag = FALSE
species_predON = args[1]
replicate = args[2]

tmp_file_path=""


model_desc_modify = ""


# if (fullModel==TRUE){
#     model_desc_modify = paste(model_desc_modify,"_fullModel",sep="")
# }
# if (equiv_toLowest==TRUE){
#     model_desc_modify = paste(model_desc_modify,"_equiv_toLowest",sep="")}
# if (exclude_CpG==TRUE){
#     model_desc_modify = paste(model_desc_modify,"_noCpG",sep="")}
# if (exclude_triplet==TRUE){
#     model_desc_modify = paste(model_desc_modify,"_noTriplets",sep="")}
# if (null_model==TRUE){
#     model_desc_modify = paste(model_desc_modify,"_null_model",sep="")}


load(paste(tmp_file_path,"data/model/model.RData",sep=""))#model
# load(paste(tmp_file_path,"data/model/samples_sites_test.RData",sep=""))#sample_sites_test
all_data <- read.csv(paste(tmp_file_path,"data/data_dataframes/",species_predON,"/all_data_readyForPrediction_",species_predON,"_",replicate,".csv",sep=""),header=TRUE)


#rerplacing with "not_transcribed" is fine as they are the last level before non_transcribed in the annotation module     
probs <- predict.glm(model, all_data, type="response")
probs_df <- data.frame(x = probs)
#binding the predicted proabilities to the OG data (but the test sites only)
probs_df <- probs_df %>%
    bind_cols(all_data)
#renaming the glm_probs coloumn
colnames(probs_df ) <- replace(colnames(probs_df ), 1, "glm_probs")
#writing to file 
filename = paste(tmp_file_path,"data/data_dataframes/",species_predON,"/ProbabilityDf",species_predON,"_",replicate,".csv",sep="")#this sep is for the filename string
write.csv(probs_df ,filename, row.names = FALSE)

