library(SCENIC)
library(reshape2)
library(pheatmap)
library(htmlwidgets)
library(Seurat)


save_pheatmap_pdf <- function(x, filename, width=8, height=30) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   pdf(filename, width=width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}

loadRData <- function(fileName){
    load(fileName)
    get(ls()[ls() != "fileName"])
}


runScenic<-function(seur_input, file_output, fileloc_cells, animals, clusters, database, cores){

    stopifnot(!missing(seur_input))
    stopifnot(!missing(file_output))
    stopifnot(!missing(fileloc_cells))
    stopifnot(!missing(animals))
    stopifnot(!missing(clusters))
    stopifnot(!missing(cores))

    dir.create(file_output, showWarnings=FALSE)
    setwd(file_output)

    seur_temp<-loadRData(seur_input)

    seur_temp$animal_type<-substr(seur_temp@meta.data$sample_order,1,2)

    if (length(animals) == 0){
      seur_temp=seur_temp
      }else{
      seur_temp<-subset(seur_temp, animal_type %in% animals)
      }

    if (length(clusters) == 0){
      seur_temp=seur_temp
    }else{
      seur_temp<-SubsetData(seur_temp, ident.use=clusters)
    }


    exprMat<-as.matrix(GetAssayData(seur_temp, slot="counts"))
    cellInfo<-data.frame(seuratCluster=Idents(seur_temp))


    scenicOptions <- initializeScenic(org="mgi", dbDir=database, nCores=cores)

    genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                               minCountsPerGene=3*.01*ncol(exprMat),
                               minSamples=ncol(exprMat)*.01)
    exprMat_filtered <- exprMat[genesKept, ]
    runCorrelation(exprMat_filtered, scenicOptions)
    exprMat_filtered_log <- log2(exprMat_filtered+1)
    runGenie3(exprMat_filtered_log, scenicOptions)


    ### Build and score the GRN
    exprMat_log <- log2(exprMat+1)
    runSCENIC_1_coexNetwork2modules(scenicOptions)
    runSCENIC_2_createRegulons(scenicOptions)
    runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
    runSCENIC_4_aucell_binarize(scenicOptions)

    motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")

    saveWidget(viewMotifs(motifEnrichment_selfMotifs_wGenes) , paste0(fileloc_cells, "_motifs.v1.html"))

    regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
    write.table(as.data.frame(regulonTargetsInfo), file=paste0(fileloc_cells, "_regulon_targets.v1.txt"), sep="\t", quote=F)
    saveWidget(viewMotifs(regulonTargetsInfo), paste0(fileloc_cells, "_regulon_targets.v1.html"))

    regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
    regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
    regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$seuratCluster),
                                         function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
    regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

    ph<-pheatmap(regulonActivity_byCellType_Scaled, #fontsize_row=3,
                       color=colorRampPalette(c("blue","white","deeppink4"))(100), breaks=seq(-3, 3, length.out = 100),
                       treeheight_row=10, treeheight_col=10, border_color=NA, silent=TRUE)

    save_pheatmap_pdf(ph, paste0(fileloc_cells, "_regulon_heatmap.v1.pdf"))

    topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
    colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
    topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]

    write.table(topRegulators, file=paste0(fileloc_cells, "_top_regulators.v1.txt"),sep="\t", quote=F)


    minPerc <- .7
    binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
    cellInfo_binarizedCells <- cellInfo[which(rownames(cellInfo)%in% colnames(binaryRegulonActivity)),, drop=FALSE]
    regulonActivity_byCellType_Binarized <- sapply(split(rownames(cellInfo_binarizedCells), cellInfo_binarizedCells$seuratCluster),
                                                   function(cells) rowMeans(binaryRegulonActivity[,cells, drop=FALSE]))
    binaryActPerc_subset <- regulonActivity_byCellType_Binarized[which(rowSums(regulonActivity_byCellType_Binarized>minPerc)>0),]

    pb<-pheatmap::pheatmap(binaryActPerc_subset, # fontsize_row=5,
                       color = colorRampPalette(c("white","pink","deeppink4"))(100), breaks=seq(0, 1, length.out = 100),
                       treeheight_row=10, treeheight_col=10, border_color=NA, silent=TRUE)

    save_pheatmap_pdf(pb, file=paste0(fileloc_cells, "_binary_regulon_heatmap_v1.pdf"))

    topRegulatorsBinary <- reshape2::melt(regulonActivity_byCellType_Binarized)
    colnames(topRegulatorsBinary) <- c("Regulon", "CellType", "RelativeActivity")
    topRegulatorsBinary <- topRegulatorsBinary[which(topRegulatorsBinary$RelativeActivity>minPerc),]

    write.table(topRegulatorsBinary, file=paste0(fileloc_cells, "_top_binary_regulators.v1.txt"), sep="\t", quote=F)

    save.image(paste0(fileloc_cells, ".v1.RData"))

}
