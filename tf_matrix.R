# Load packages
library(dplyr)
library(reshape2)

# Variables
animal<-"OY"
input<-paste0("/Users/kavyashah/Harvard Drive/Harvard/Rubin Lab/R/top_regulator_files/by_animal/", animal)
output<-"/Users/kavyashah/Harvard Drive/Harvard/Rubin Lab/R/TF_matrices"

# Load data (SCENIC top regulator files)
# Inputs:
# cell lineage
# animal type
# cluster type (big or small)
load_data <- function(lineage, animal, cluster_type){
  file_name <- paste0(lineage, "_", animal, "_", cluster_type, "_top_regulators.v1.txt")
  read.delim(file_name, sep="\t", header=T)
}

clean_regulons <- function(data) {
  # Clean up the Regulon names
  cleaned<-sub(" .*", "",data$Regulon)
  cleaned<-sub("_[^_]+$", "",cleaned)
  data$Regulon<-cleaned
  
  # Turn the regulons for each animal into character vectors
  regs<-data$Regulon
  return(regs)
}

# Function that returns master matrix of regulon scores by cluster
# Inputs:
# animal type
# input directory
# output directory
make_tf_matrix <- function(animal, input, output){
  # Set input directory for loading in data
  setwd(input)
  
  # Load in data
  vasc <- load_data("VASC", animal, "big")
  immune <- load_data("IMMUNE", animal, "big")
  neuron <- load_data("NEURON", animal, "big")
  asc_epc <- load_data("ASC_EPC", animal, "big")
  # olg <- load_data("OLG", animal, "big")
  
  # Set output directory for final matrix
  setwd(output)
  marker<-read.delim("marker_file.txt", sep="\t", header=T)
  
  # subset out OLG cell types (for now)
  marker<-marker[!(marker$Cluster=="OPC" | marker$Cluster=="OLG" | marker$Cluster=="OEG"),]
  
  # Bind data
  all_lin <- bind_rows(vasc, immune, neuron, asc_epc)
  
  # Clean regulon names
  all_lin$Regulon <- clean_regulons(all_lin)
  
  # Change df dimensions
  all_lin<- dcast(all_lin, Regulon ~ CellType, value.var = "RelativeActivity")
  # rownames(all_lin)<-all_lin$Regulon
  # all_lin$Regulon<-NULL
  
  # Count regulon frequency and store as col
  a<-1:nrow(all_lin)
  counts <- lapply(a, function(x) length(which(!is.na(all_lin[x,]))) - 1)
  counts<-as.numeric(counts)
  all_lin$counts<-counts
  
  # Sort df by counts
  all_lin<-all_lin[order(counts, decreasing=TRUE),]
  
  # Reorder cols by marker file
  all_lin<-all_lin %>% select(marker$Cluster)
  
  # Save matrix
  file_name<-paste0(animal, "_Regulon_Activity_Matrix.txt")
  write.table(all_lin, file=file_name, sep="\t", row.names=F, quote=F)
  
  # Reset wd
  setwd(input)
}
