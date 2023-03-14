# Load libraries
library(VennDiagram)

# Function to get regulons
get_regulons <- function(file_name, animal) {
  # Read in the Top regulators data
  data<-read.delim(file_name, sep="\t", header=T)
  
  # Subset ECs from Scenic Data
  subset<-subset(data, data$CellType %in% "EC")
  
  # Clean up the Regulon names
  cleaned<-sub(" .*", "",subset$Regulon)
  cleaned<-sub("_[^_]+$", "",cleaned)
  subset$Regulon<-cleaned
  
  # Turn the regulons for each animal into character vectors
  regs<-subset$Regulon
  return(regs)
}

# run function on animal types
oy_r<-get_regulons("VASC_OY_big_top_regulators.v1.txt", "OY")
oo_r<-get_regulons("VASC_OO_big_top_regulators.v1.txt", "OO")
ox_r<-get_regulons("VASC_OX_big_top_regulators.v1.txt", "OX")
yo_r<-get_regulons("VASC_YO_big_top_regulators.v1.txt", "YO")
yy_r<-get_regulons("VASC_YY_big_top_regulators.v1.txt", "YY")
yx_r<-get_regulons("VASC_YX_big_top_regulators.v1.txt", "YX")

# Put the character vectors into a list
paracoolo_reg_list<-list(OY=oy_r, OO=oo_r, OX=ox_r, YX=yx_r)
paracooly_reg_list<-list(YO=yo_r, YY=yy_r, OX=ox_r, YX=yx_r)

# ParaCoolY Venn 
vd<- venn.diagram(lwd= 5, paracooly_reg_list, "VASC_big_regulons_ParaCoolY.png", main="Significant regulons identified by SCENIC", sub="VASC Lineage", main.fontface = "bold", main.fontfamily = "Helvetica",sub.fontfamily = "Helvetica", cat.fontface = "bold", cat.fontfamily = "Helvetica", fill=c("blue", "red", "green", "yellow"), fontfamily="Helvetica", alpha=rep(0.4,4), imagetype = "png")


# ParaCoolO Venn 
vd<- venn.diagram(lwd= 5, paracoolo_reg_list, "VASC_big_regulons_ParaCoolO.png", main="Significant regulons identified by SCENIC", sub="VASC Lineage", main.fontface = "bold", main.fontfamily = "Helvetica",sub.fontfamily = "Helvetica", cat.fontface = "bold", cat.fontfamily = "Helvetica", fill=c("blue", "red", "green", "yellow"), fontfamily="Helvetica", alpha=rep(0.4,4), imagetype = "png")

# Find the TFs that only intersect among OY and YX but not in OO and OX 
oy_yx_overlap<-intersect(oy_r,yx_r)
oy_yx_no_ox<-setdiff(oy_yx_overlap,ox_r)
oy_yx_no_ox_no_oo<-setdiff(oy_yx_no_ox, oo_r)

# Find the TFs that only intersect among YO and OX but not in YY and YX 
yo_ox_overlap<-intersect(yo_r,ox_r)
yo_ox_no_yx<-setdiff(yo_ox_overlap,yx_r)
yo_ox_no_yx_no_yy<-setdiff(yo_ox_no_yx, yy_r)
