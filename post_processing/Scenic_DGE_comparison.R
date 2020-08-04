library(openxlsx)

oy<-read.delim("VASC_OY_big_top_regulators.v1.txt", sep="\t", header=T)
oo<-read.delim("VASC_OO_big_top_regulators.v1.txt", sep="\t", header=T)
ox<-read.delim("VASC_OX_big_top_regulators.v1.txt", sep="\t", header=T)
yx<-read.delim("VASC_YX_big_top_regulators.v1.txt", sep="\t", header=T)
yo<-read.delim("VASC_YO_big_top_regulators.v1.txt", sep="\t", header=T)
yy<-read.delim("VASC_YY_big_top_regulators.v1.txt", sep="\t", header=T)
scenic<-rbind(oy, oo, ox, yx, yo, yy)

ec<-read.xlsx("ALL_big_paracoolo.v1.xlsx", sheet="EC")

# Subset ECs from Scenic data
scenic_ec<-subset(scenic, scenic$CellType %in% "EC")

# clean up regulon names
scenic_ec_clean<-sub(" .*", "", scenic_ec$Regulon)
scenic_ec_clean<-sub("_[^_]+$", "", scenic_ec_clean)
scenic_ec$Regulon<-scenic_ec_clean

# Restrict DGEs to adjusted p-vals <= 0.05
sig_ec<-subset(ec, p_adj.loc<=0.1 | grepl("-", p_adj.loc))

# Intersect DGEs and Scenic Data
tf<-intersect(sig_ec$gene, scenic_ec$Regulon)
