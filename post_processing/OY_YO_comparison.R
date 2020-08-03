# load packages 

library(ggplot2)
library(pheatmap)
library(openxlsx)

# Read in the top regulators
oy<-read.delim("VASC_OY_big_top_regulators.v1.txt", sep="\t", header=T)
yo<-read.delim("VASC_YO_big_top_regulators.v1.txt", sep="\t", header=T)

paracoolo<-read.xlsx("ALL_big_paracoolo.v1.xlsx", sheet="EC")
paracooly<-read.xlsx("ALL_big_paracooly.v1.xlsx", sheet="EC")

# Subset ECs from Scenic Data
oy_ec<-subset(oy, oy$CellType %in% "EC")
yo_ec<-subset(yo, yo$CellType %in% "EC")

# Clean up the Regulon names
oy_ec_clean<-sub(" .*", "",oy_ec$Regulon)
oy_ec_clean<-sub("_[^_]+$", "",oy_ec_clean)
oy_ec$Regulon<-oy_ec_clean

yo_ec_clean<-sub(" .*", "", yo_ec$Regulon)
yo_ec_clean<-sub("_[^_]+$", "", yo_ec_clean)
yo_ec$Regulon<-yo_ec_clean

# Turn the regulons for each animal into character vectors
oy_r=oy_ec$Regulon
yo_r=yo_ec$Regulon

# Put the character vectors into a list
reg_list<-list(OY=oy_r, YO=yo_r)

# Plot Venn Diagram
vd<- venn.diagram(lwd= 5, reg_list, "VASC_big_OY_YO", main="Significant regulons in OY and YO Mice", sub="VASC Lineage", main.fontface = "bold", main.fontfamily = "Helvetica",sub.fontfamily = "Helvetica", cat.fontface = "bold", cat.fontfamily = "Helvetica", fill=c("blue", "red"), fontfamily="Helvetica", alpha=rep(0.4,2), imagetype="png")

# Identify the 101 TFs Significant in OY and YO
oy_yo_overlap<-intersect(oy_r,yo_r)

# Make a directionality table for OY and YO
pruned_paracoolo<-subset(paracoolo, paracoolo$gene %in% oy_yo_overlap)
pruned_paracooly<-subset(paracooly, paracooly$gene %in% oy_yo_overlap)

oy_yo_overlap=sort(oy_yo_overlap)
pruned_paracoolo<-pruned_paracoolo[order(pruned_paracoolo$gene),]
pruned_paracooly<-pruned_paracooly[order(pruned_paracooly$gene),]

pco_logfc<-pruned_paracoolo$logFC
pco_logfc<-as.numeric(pco_logfc)
pcy_logfc<-pruned_paracooly$logFC
pcy_logfc<-as.numeric(pcy_logfc)


table<-cbind(row.names(oy_yo_overlap),pco_logfc,pcy_logfc)
row.names(table)<-oy_yo_overlap

colnames(table)[1] = "ParaCoolO"
colnames(table)[2] = "ParaCoolY"

# Make heatmap for directionality table
mtable<-as.matrix(table)

ph<-pheatmap(mtable, #fontsize_row=3,
             color=colorRampPalette(c("blue","white","deeppink4"))(100),      
             scale="none",
             treeheight_row=10, treeheight_col=10, 
             border_color=NA, silent=TRUE, clustering_distance_rows="euclidean",
             clustering_method="ward.D2", cluster_row=T)

pdf("OY_YO_bidirectionality.pdf", width=8, height=15)
grid::grid.newpage()
grid::grid.draw(ph$gtable)

#closes the plot
dev.off()
