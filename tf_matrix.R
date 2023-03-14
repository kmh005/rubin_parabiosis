# Kavya Shah
# Last update: 08.20.2020 5:25 PM
#
# This script organizes the SCENIC top regulators data into a parseable format.
# The output is, for each animal type, a single matrix. Each row displays activity values
# for a TF by cell type. This matrix is useful for comparing TF activity across cell types.
#
# Usage instructions:
# 1. Load the libraries
# 2. Load all three functions (load_data, clean_regulons, make_tf_matrix)
# 3. Set up variables for directories and load them (example set of vars provided)
# 5. Execute all_animals at the bottom and run the lapply function to get all the matrices

# Load packages
library(dplyr)
library(reshape2)
library(tibble)

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
  olg <- load_data("OLG", animal, "big")
  
  # Set output directory for final matrix
  setwd(output)
  marker<-read.delim("marker_file.txt", sep="\t", header=T)
  # marker[nrow(marker) + 1,] = c("counts","counts")
  
  
  # HypENC and TNC Cell types don't exist for ASC_EPC YX and OO
  if (animal == "OO" | animal == "YX") {
    marker<-marker[!(marker$Cluster=="HypEPC" | marker$Cluster=="TNC"),]
  }
  
  # Bind data
  all_lin <- bind_rows(vasc, immune, neuron, asc_epc, olg)
  
  # Clean regulon names
  all_lin$Regulon <- clean_regulons(all_lin)
  
  # Change df dimensions
  all_lin<- dcast(all_lin, Regulon ~ CellType, value.var = "RelativeActivity")
  rownames(all_lin)<-all_lin$Regulon
  all_lin$Regulon<-NULL
  
  # Reorder cols by marker file
  all_lin<-all_lin %>% rownames_to_column("Regulon") %>% select(c("Regulon", marker$Cluster))

  # Count regulon frequency and store as col
  counts<- sapply(1:nrow(all_lin), function(x) length(which(!is.na(all_lin[x,]))) - 1)
  all_lin$counts<-counts
  
  # Sort df by counts
  all_lin<-all_lin[order(counts, decreasing=TRUE),]

  # Save matrix
  file_name<-paste0(animal, "_Regulon_Activity_Matrix.txt")
  write.table(all_lin, file=file_name, sep="\t", row.names=F, quote=F)
  
  return(all_lin)
  # Reset wd
  setwd(input)
}
                  
# Variables
input<-"/Users/kavyashah/Harvard Drive/Harvard/Rubin Lab/R/top_regulator_files/all/"
output<-"/Users/kavyashah/Harvard Drive/Harvard/Rubin Lab/R/TF_matrices"
                  
# Run script
all_animals<-c("OY", "OX", "OO", "YO", "YX", "YY")
all_results<-lapply(all_animals, function(animal) make_tf_matrix(animal, input, output))
