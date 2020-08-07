# import libraries
library(ggplot2)
library(pheatmap)
library(openxlsx)
library(VennDiagram)
library(dplyr)
library(futile.logger)

# load data (SCENIC top regulator files)
# Inputs:
# cell lineage
# animal type
# cluster type (big or small)
load_data <- function(lineage, animal, cluster_type){
  file_name <- paste0(lineage, "_", animal, "_", cluster_type, "_top_regulators.v1.txt")
  read.delim(file_name, sep="\t", header=T)
}

# Variables
oy<-load_data("VASC", "OY", "big")
oo<-load_data("VASC", "OO", "big")
ox<-load_data("VASC", "OX", "big")
yx<-load_data("VASC", "YX", "big")
yo<-load_data("VASC", "YO", "big")
yy<-load_data("VASC", "YY", "big")
lineage<-"VASC"
cluster_type<-"big"
paracoolo<-"ALL_big_paracoolo.v1.xlsx"
paracooly<-"ALL_big_paracooly.v1.xlsx"
output<-"/Users/kavyashah/Harvard Drive/Harvard/Rubin Lab/R/VASC/"
# cluster<-"PC"
pval<-0.05
input<-"/Users/kavyashah/Harvard Drive/Harvard/Rubin Lab/R/DGE_Scenic_comparison"

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
                    sub="VASC Lineage", main.fontface = "bold", main.fontfamily = "Helvetica",sub.fontfamily = "Helvetica", 
                    cat.fontface = "bold", cat.fontfamily = "Helvetica", fill=c("blue", "red", "green", "yellow"), 
                    fontfamily="Helvetica", alpha=rep(0.4,4), imagetype = "png")
  
  # ParaCoolO Venn 
  vd_name<-paste0(lineage, "_regulons_ParaCoolO.png")
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  vd<- venn.diagram(lwd= 5, paracoolo_reg_list, vd_name, main="Significant regulons identified by SCENIC", 
                    sub="VASC Lineage", main.fontface = "bold", main.fontfamily = "Helvetica",sub.fontfamily = "Helvetica", 
                    cat.fontface = "bold", cat.fontfamily = "Helvetica", fill=c("blue", "red", "green", "yellow"), 
                    fontfamily="Helvetica", alpha=rep(0.4,4), imagetype = "png")
  
  # OY-YO comparison Venn
  vd_name<-paste0(lineage, "_OY_YO.png")
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  vd<- venn.diagram(lwd= 5, oy_yo_reg_list, vd_name, main="Significant regulons in OY and YO Mice", sub="VASC Lineage", 
                    main.fontface = "bold", main.fontfamily = "Helvetica",sub.fontfamily = "Helvetica", cat.fontface = "bold", 
                    cat.fontfamily = "Helvetica", fill=c("blue", "red"), fontfamily="Helvetica", alpha=rep(0.4,2), imagetype="png")
  
  # Find the TFs that only intersect among OY and YX but not in OO and OX 
  # "Rejuvenation" genes
  oy_yx_overlap<-intersect(oy_regs,yx_regs)
  oy_yx_no_ox<-setdiff(oy_yx_overlap,ox_regs)
  oy_yx_no_ox_no_oo<-setdiff(oy_yx_no_ox, oo_regs)
  file_name<-file(paste0(cluster, "_rejuvenation_genes.txt"))
  writeLines(oy_yx_no_ox_no_oo, file_name)
  close(file_name)
  
  # Find the TFs that only intersect among YO and OX but not in YY and YX 
  # "Aging acceleration" genes
  yo_ox_overlap<-intersect(yo_regs,ox_regs)
  yo_ox_no_yx<-setdiff(yo_ox_overlap,yx_regs)
  yo_ox_no_yx_no_yy<-setdiff(yo_ox_no_yx, yy_regs)
  file_name<-file(paste0(cluster, "_aging_acceleration_genes.txt"))
  writeLines(yo_ox_no_yx_no_yy, file_name)
  close(file_name)
  
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
  file_name<-file(paste0(cluster, "_ParaCoolO_OY_comparison.txt"))
  writeLines(oy_pco, file_name)
  close(file_name)
  
  oo_pco<-intersect(sig_dge_pco$gene, oo_regs)
  file_name<-file(paste0(cluster, "_ParaCoolO_OO_comparison.txt"))
  writeLines(oo_pco, file_name)
  close(file_name)
  
  ox_pco<-intersect(sig_dge_pco$gene, ox_regs)
  file_name<-file(paste0(cluster, "_ParaCoolO_OX_comparison.txt"))
  writeLines(ox_pco, file_name)
  close(file_name)
  
  yx_pco<-intersect(sig_dge_pco$gene, yx_regs)
  file_name<-file(paste0(cluster, "_ParaCoolO_YX_comparison.txt"))
  writeLines(yx_pco, file_name)
  close(file_name)
  
  # ParaCoolY
  yo_pcy<-intersect(sig_dge_pcy$gene, yo_regs)
  file_name<-file(paste0(cluster, "_ParaCoolY_YO_comparison.txt"))
  writeLines(yo_pcy, file_name)
  close(file_name)
  
  yy_pcy<-intersect(sig_dge_pcy$gene, yy_regs)
  file_name<-file(paste0(cluster, "_ParaCoolY_YY_comparison.txt"))
  writeLines(yy_pcy, file_name)
  close(file_name)
  
  yx_pcy<-intersect(sig_dge_pcy$gene, yx_regs)
  file_name<-file(paste0(cluster, "_ParaCoolY_YX_comparison.txt"))
  writeLines(yx_pcy, file_name)
  close(file_name)
  
  ox_pcy<-intersect(sig_dge_pcy$gene, ox_regs)
  file_name<-file(paste0(cluster, "_ParaCoolY_OX_comparison.txt"))
  writeLines(ox_pcy, file_name)
  close(file_name)

  # Reset directory
  setwd(input)
}

# Run script
# all_clusters<-c("MG", "MAC", "MNC", "T_cell", "DC", "NK", "B_cell", "NEUT")
all_clusters<-c("EC", "PC", "VSMC", "VLMC", "Hb_VC", "ABC")
all_results<-lapply(all_clusters, function(x) post_processing(x, oy, oo, ox, yx, yo, yy, lineage, paracoolo, paracooly, output, animals, pval, input))
