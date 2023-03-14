# Kavya Shah
# Last update: 08.20.2020 5:25 PM
#
# This script several post-processing functions using SCENIC top regulators data
# and ParaCool DGE data, including making a Venn diagram of the top regulons from 
# SCENIC, comparing bidirectionality of expression in OY and YO animals, and finding
# TFs that are DGE and SCENIC-significant.
#
# Usage instructions:
# 1. Load the libraries
# 2. Load all three functions (load_data, get_regulons, post_processing)
# 3. Set up variables and load them (example set of vars provided)
# 4. At the bottom, un-comment out the all_clusters list for the lineage you are 
# interested in
# 5. Run the lapply function to get output for all cell clusters within a lineage

# Load libraries
library(ggplot2)
library(pheatmap)
library(openxlsx)
library(VennDiagram)
library(dplyr)
library(futile.logger)

# Load data (SCENIC top regulator files)
# Inputs:
# cell lineage
# animal type
# cluster type (big or small)
load_data <- function(lineage, animal, cluster_type){
  file_name <- paste0(lineage, "_", animal, "_", cluster_type, "_top_regulators.v1.txt")
  read.delim(file_name, sep="\t", header=T)
}

# Function to get regulons for each animal type and clean up their names
# Inputs: 
# Scenic regulon data
# Cluster we want
get_regulons <- function(data, cluster) {
  # Subset cluster from Scenic Data
  subset<-subset(data, data$CellType %in% cluster)
  
  # Clean up the Regulon names
  cleaned<-sub(" .*", "",subset$Regulon)
  cleaned<-sub("_[^_]+$", "",cleaned)
  subset$Regulon<-cleaned
  
  # Turn the regulons for each animal into character vectors
  regs<-subset$Regulon
  return(regs)
}

# Wrapper function for post-processing
# Inputs:
# SCENIC top regulators data
# cell lineage
# paraCoolO file name
# paraCoolY file name
# file output directory folder
# cluster name (like EC)
# p-value
post_processing <- function(cluster, oy, oo, ox, yx, yo, yy, lineage, paracoolo, paracooly, output, animals, pval, input){
  
  # Get list of clean regulon names for Scenic runs for each animal
  oy_regs<-get_regulons(oy, cluster)
  oo_regs<-get_regulons(oo, cluster)
  ox_regs<-get_regulons(ox, cluster)
  yo_regs<-get_regulons(yo, cluster)
  yy_regs<-get_regulons(yy, cluster)
  yx_regs<-get_regulons(yx, cluster)
  
  # upload ParaCoolO and ParaCoolY files
  paracoolo<-read.xlsx(paracoolo, sheet=cluster)
  paracooly<-read.xlsx(paracooly, sheet=cluster)
  
  # Create directory
  dir_name<-paste0(output, cluster, "/")
  dir.create(dir_name, showWarnings=FALSE)
  setwd(dir_name)
  
  # Make regulon lists for Venn diagrams
  paracoolo_reg_list<-list(OY=oy_regs, OO=oo_regs, OX=ox_regs, YX=yx_regs)
  paracooly_reg_list<-list(YO=yo_regs, YY=yy_regs, OX=ox_regs, YX=yx_regs)
  oy_yo_reg_list<-list(OY=oy_regs, YO=yo_regs)
  
  # Plot Venn diagrams
  # ParaCoolY Venn 
  vd_name<-paste0(lineage, "_regulons_ParaCoolY.png")
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  vd<- venn.diagram(lwd= 5, paracooly_reg_list, vd_name, main="Significant regulons identified by SCENIC", 
                    sub=paste0(lineage," Lineage"), main.fontface = "bold", main.fontfamily = "Helvetica",sub.fontfamily = "Helvetica", 
                    cat.fontface = "bold", cat.fontfamily = "Helvetica", fill=c("blue", "red", "green", "yellow"), 
                    fontfamily="Helvetica", alpha=rep(0.4,4), imagetype = "png")
  
  # ParaCoolO Venn 
  vd_name<-paste0(lineage, "_regulons_ParaCoolO.png")
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  vd<- venn.diagram(lwd= 5, paracoolo_reg_list, vd_name, main="Significant regulons identified by SCENIC", 
                    sub=paste0(lineage," Lineage"), main.fontface = "bold", main.fontfamily = "Helvetica",sub.fontfamily = "Helvetica", 
                    cat.fontface = "bold", cat.fontfamily = "Helvetica", fill=c("blue", "red", "green", "yellow"), 
                    fontfamily="Helvetica", alpha=rep(0.4,4), imagetype = "png")
  
  # OY-YO comparison Venn
  vd_name<-paste0(lineage, "_OY_YO.png")
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  vd<- venn.diagram(lwd= 5, oy_yo_reg_list, vd_name, main="Significant regulons in OY and YO Mice", sub=paste0(lineage," Lineage"), 
                    main.fontface = "bold", main.fontfamily = "Helvetica",sub.fontfamily = "Helvetica", cat.fontface = "bold", 
                    cat.fontfamily = "Helvetica", fill=c("blue", "red"), fontfamily="Helvetica", alpha=rep(0.4,2), imagetype="png")
  
  # Find the TFs that only intersect among OY and YX but not in OO and OX 
  # "Rejuvenation" genes
  oy_yx_overlap<-intersect(oy_regs,yx_regs)
  oy_yx_no_ox<-setdiff(oy_yx_overlap,ox_regs)
  oy_yx_no_ox_no_oo<-setdiff(oy_yx_no_ox, oo_regs)
  write.table(oy_yx_no_ox_no_oo, file=paste0(cluster, "_rejuvenation_genes.txt"), sep="\t", row.names=F, quote=F)
  
  # Find the TFs that only intersect among YO and OX but not in YY and YX 
  # "Aging acceleration" genes
  yo_ox_overlap<-intersect(yo_regs,ox_regs)
  yo_ox_no_yx<-setdiff(yo_ox_overlap,yx_regs)
  yo_ox_no_yx_no_yy<-setdiff(yo_ox_no_yx, yy_regs)
  write.table(yo_ox_no_yx_no_yy, file=paste0(cluster, "_aging_acceleration_genes.txt"), sep="\t", row.names=F, quote=F)
  
  # Plot directionality heatmap of TFs significant in OY and YO
  # Identify the 101 TFs Significant in OY and YO
  oy_yo_overlap<-intersect(oy_regs,yo_regs)
  
  # Make a directionality table for OY and YO
  pruned_paracoolo<-subset(paracoolo, paracoolo$gene %in% oy_yo_overlap)
  pruned_paracooly<-subset(paracooly, paracooly$gene %in% oy_yo_overlap)
  
  pco_logfc<-data.frame(gene=pruned_paracoolo$gene, ParaCoolO=pruned_paracoolo$logFC)
  pcy_logfc<-data.frame(gene=pruned_paracooly$gene, ParaCoolY=pruned_paracooly$logFC)
  
  # Merge on gene
  df<-merge(pco_logfc, pcy_logfc,by="gene")
  
  # Strip gene names
  genes<-df$gene
  df$gene<-NULL
  
  # make numeric
  df$ParaCoolO<-as.numeric(df$ParaCoolO)
  df$ParaCoolY<-as.numeric(df$ParaCoolY)
  
  # add in gene names as row names 
  row.names(df)<-genes
  mtable<-as.matrix(df)
  
  ph<-pheatmap(mtable, fontsize_row=7,
               color=colorRampPalette(c("blue","white","deeppink4"))(100),      
               scale="none",
               treeheight_row=10, treeheight_col=10, 
               border_color=NA, silent=TRUE, clustering_distance_rows="euclidean",
               clustering_method="ward.D2", cluster_row=T)
  
  pdf(paste0(cluster, "_OY_YO_bidirectionality.pdf"), width=8.5, height=11)
  grid::grid.newpage()
  grid::grid.draw(ph$gtable)
  
  #closes the plot
  dev.off()

  # Restrict DGEs to adjusted p-vals <= 0.05
  sig_dge_pco<-subset(paracoolo, p_adj.loc<=pval | grepl("-", p_adj.loc))
  sig_dge_pcy<-subset(paracooly, p_adj.loc<=pval | grepl("-", p_adj.loc))
  
  # Intersect DGEs and Scenic Data
  # ParaCoolO
  oy_pco<-intersect(sig_dge_pco$gene, oy_regs)
  write.table(oy_pco, file=paste0(cluster, "_ParaCoolO_OY_comparison.txt"), sep="\t", row.names=F, quote=F)
  
  oo_pco<-intersect(sig_dge_pco$gene, oo_regs)
  write.table(oo_pco, file=paste0(cluster, "_ParaCoolO_OO_comparison.txt"), sep="\t", row.names=F, quote=F)
  
  ox_pco<-intersect(sig_dge_pco$gene, ox_regs)
  write.table(ox_pco, file=paste0(cluster, "_ParaCoolO_OX_comparison.txt"), sep="\t", row.names=F, quote=F)
  
  yx_pco<-intersect(sig_dge_pco$gene, yx_regs)
  write.table(yx_pco, file=paste0(cluster, "_ParaCoolO_YX_comparison.txt"), sep="\t", row.names=F, quote=F)
  
  # ParaCoolY
  yo_pcy<-intersect(sig_dge_pcy$gene, yo_regs)
  write.table(yo_pcy, file=paste0(cluster, "_ParaCoolY_YO_comparison.txt"), sep="\t", row.names=F, quote=F)
  
  yy_pcy<-intersect(sig_dge_pcy$gene, yy_regs)
  write.table(yy_pcy, file=paste0(cluster, "_ParaCoolY_YY_comparison.txt"), sep="\t", row.names=F, quote=F)
  
  yx_pcy<-intersect(sig_dge_pcy$gene, yx_regs)
  write.table(yx_pcy, file=paste0(cluster, "_ParaCoolY_YX_comparison.txt"), sep="\t", row.names=F, quote=F)
  
  ox_pcy<-intersect(sig_dge_pcy$gene, ox_regs)
  write.table(ox_pcy, file=paste0(cluster, "_ParaCoolY_OX_comparison.txt"), sep="\t", row.names=F, quote=F)
  
  # Combine all intersects into a data frame
  intersects_pco<-list(oy=oy_pco, ox=ox_pco, yx=yx_pco, oo=oo_pco)
  intersects_pco<-lapply(intersects_pco, function(gene) as.data.frame(gene))
  intersects_pco<-bind_rows(intersects_pco, .id="mouse")
  write.table(intersects_pco, file="ParaCoolO_scenic_intersects.txt", sep="\t", row.names=F, quote=F)
  
  intersects_pcy<-list(oo=yo_pcy, ox=ox_pcy, yx=yx_pcy, yy=yy_pcy)
  intersects_pcy<-lapply(intersects_pcy, function(gene) as.data.frame(gene))
  intersects_pcy<-bind_rows(intersects_pcy, .id="mouse")
  write.table(intersects_pcy, file="ParaCoolY_scenic_intersects.txt", sep="\t", row.names=F, quote=F)
  
  # plot 
  title<-paste0("ParaCoolO and Scenic comparison, ", lineage, " lineage, ", cluster, " cluster")
  g<-ggplot(intersects_pco, aes(x=gene, fill=mouse)) + geom_bar() + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))
  ggsave("ParaCoolO_Scenic_comparison.png", height=7, width=7, g)
  
  title<-paste0("ParaCoolY and Scenic comparison, ", lineage, " lineage, ", cluster, " cluster")
  g<-ggplot(intersects_pcy, aes(x=gene, fill=mouse)) + geom_bar() + ggtitle(title)
  # try out + labs(title="Blah", x="blah")
  ggsave("ParaCoolY_Scenic_comparison.png", height=7, width=7, g)

  # Reset directory
  setwd(input)
}
                         
# Load variables
input<-"/Users/kavyashah/Harvard Drive/Harvard/Rubin Lab/R/top_regulator_files/by_lineage/ASC_EPC/"
output<-"/Users/kavyashah/Harvard Drive/Harvard/Rubin Lab/R/post_processing_files/ASC_EPC/"
oy<-load_data("ASC_EPC", "OY", "big")
oo<-load_data("ASC_EPC", "OO", "big")
ox<-load_data("ASC_EPC", "OX", "big")
yx<-load_data("ASC_EPC", "YX", "big")
yo<-load_data("ASC_EPC", "YO", "big")
yy<-load_data("ASC_EPC", "YY", "big")
lineage<-"ASC_EPC"
cluster_type<-"big"
paracoolo<-"ALL_big_paracoolo.v1.xlsx"
paracooly<-"ALL_big_paracooly.v1.xlsx"
cluster<-"EC"
pval<-0.05
                         
# Run script
# IMMUNE clusters
# all_clusters<-c("MG", "MAC", "MNC", "T_cell", "DC", "NK", "B_cell", "NEUT")
# OLG clusters
# all_clusters<-c("OPC", "OLG", "OEG")
# mNEUR small clusters
# all_clusters<-c("GABA", "DOPA", "GLUT", "CHOL")
# NEURON clusters
# all_clusters<-c("GABA", "DOPA", "GLUT", "CHOL", "NendC", "ImmN", "NRP")
# ASC_EPC clusters
# all_clusters<-c("ASC", "NSC", "ARP", "CPC", "EPC")
# VASC clusters
# all_clusters<-c("EC", "PC", "VSMC", "VLMC", "Hb_VC", "ABC")
all_results<-lapply(all_clusters, function(x) post_processing(x, oy, oo, ox, yx, yo, yy, lineage, paracoolo, paracooly, output, animals, pval, input))
