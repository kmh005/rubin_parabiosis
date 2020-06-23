# load packages
library(SCENIC)
library(reshape2)
library(pheatmap)
library(htmlwidgets)
library(Seurat)

#Example on defining variables

# I used the following input: seur_input<-"/net/rcstorenfs02/ifs/rc_labs/HSCI/Users/kavyams/20200427_Seur_R12345_Parabiosis_VASC_Lin_Small_no_doublet_reclustered_rename.Robj"
seur_input<-"/net/rcstorenfs02/ifs/rc_labs/HSCI/Users/kholton/parabiosis/R_Files/VASC_Lin_Final/20200427_Seur_R12345_Parabiosis_VASC_Lin_Small_no_doublet_reclustered_rename.Robj"
# I used the following output: "/net/rcstorenfs02/ifs/rc_labs/HSCI/Users/kavyams/scenic/VASC_Lin_Final/code_scratch/"
file_output="/net/rcstorenfs02/ifs/rc_labs/HSCI/Users/kholton/parabiosis/scenic/VASC_Lin_Final/code_scratch/"
fileloc_cells="VASC_OY_big_scratch"
animals="OY"
clusters=NULL
cores=16
# I used the following database path: database<-"/net/rcstorenfs02/ifs/rc_labs/HSCI/Users/kavyams/scenic_db"
database<-"/net/rcstorenfs02/ifs/rc_labs/HSCI/Users/kholton/parabiosis/scenic/cisTarget_databases"

#pheatmap saving function
save_pheatmap_pdf <- function(x, filename, width=8, height=30) {
    # ensures that we aren't missing function arguments
    # x contains a pheatmap object, which in turn contains a gtable heatmap
   stopifnot(!missing(x))
   stopifnot(!missing(filename))

    # creates pheatmap pdf
   pdf(filename, width=width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)

    #closes the plot
   dev.off()
}


#Robj loading function
loadRData <- function(fileName){
    # loads the file
    load(fileName)

    #returns a list of all the local variables except the file name (?)
    get(ls()[ls() != "fileName"])
}

#scenic function to annotate


runScenic<-function(seur_input, file_output, fileloc_cells, animals, clusters, database, cores){

    # ensures that we have all of the function arguments
    stopifnot(!missing(seur_input))
    stopifnot(!missing(file_output))
    stopifnot(!missing(fileloc_cells))
    stopifnot(!missing(animals))
    stopifnot(!missing(clusters))
    stopifnot(!missing(cores))

    # creates directory for file output to specified var
    # main output is in the "output" folder
    dir.create(file_output, showWarnings=FALSE)
    # sets the working directory to where we want the output to be saved
    setwd(file_output)

    # loads Seurat object with all the data
    seur_temp<-loadRData(seur_input)

    # sets the animal type (Q: what is it setting it to?)
    seur_temp$animal_type<-substr(seur_temp@meta.data$sample_order,1,2)

    # if no animal type is specified, Seurat object remains unchanged
    # if animal type is specified, subset the data to include only that type
    if (length(animals) == 0){
      seur_temp=seur_temp
    }else{
      seur_temp<-subset(seur_temp, animal_type %in% animals)
    }

    # if no cluster is specified, Seurat object remains unchanged
    # else, subset the data to only include the chosen clusters
    if (length(clusters) == 0){
      seur_temp=seur_temp
    }else{
      seur_temp<-SubsetData(seur_temp, ident.use=clusters)
    }

    # pulls the "counts" data matrix and returns it as an expression matrix (input for SCENIC)
    exprMat<-as.matrix(GetAssayData(seur_temp, slot="counts"))
    # holds cell annotation data as a data.frame; rows = cells, cols = annotation vars
    cellInfo<-data.frame(seuratCluster=Idents(seur_temp))

    # initializes SCENIC settings with mouse org
    scenicOptions <- initializeScenic(org="mgi", dbDir=database, nCores=cores)

    # gene selection based on reads per gene and num cells that have gene
    genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                               minCountsPerGene=3*.01*ncol(exprMat),
                               minSamples=ncol(exprMat)*.01)

    #filter expression matrix
    exprMat_filtered <- exprMat[genesKept, ]
    # calculate correlations between TFs and targets
    runCorrelation(exprMat_filtered, scenicOptions)
    # normalize expression matrix data
    exprMat_filtered_log <- log2(exprMat_filtered+1)
    # run GENIE3, which creates coexpression modules between TFs and targets (Takes a while to run)
    runGenie3(exprMat_filtered_log, scenicOptions)


    ### Build and score the GRN
    exprMat_log <- log2(exprMat+1)
    # Get coexpression modules
    runSCENIC_1_coexNetwork2modules(scenicOptions)
    # Get regulons with RcisTarget
    runSCENIC_2_createRegulons(scenicOptions)
    # Score regulons with AUCell
    runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
    # Cluster cells by activity/cell state
    # (Before running this, we may want to optimize thresholds for binary matrix)
    runSCENIC_4_aucell_binarize(scenicOptions)

    # Note: We aren't visualizing with tSNE or exporting to loom

    # Loads list of TF motifs that support the regulons
    motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")

    # saves a widget of the TF motifs that support the regulons, saves to "_motifs.v1.html"
    saveWidget(viewMotifs(motifEnrichment_selfMotifs_wGenes) , paste0(fileloc_cells, "_motifs.v1.html"))

    # loads the regulon targets, formats as table, saves as widget
    regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
    write.table(as.data.frame(regulonTargetsInfo), file=paste0(fileloc_cells, "_regulon_targets.v1.txt"), sep="\t", quote=F)
    saveWidget(viewMotifs(regulonTargetsInfo), paste0(fileloc_cells, "_regulon_targets.v1.html"))

    # loads regulon activity data and scales it by cell type
    regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
    regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
    regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$seuratCluster),
                                         function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
    regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

    # creates heat map of regulon acticity by cell type (scaled)
    ph<-pheatmap(regulonActivity_byCellType_Scaled, #fontsize_row=3,
                       color=colorRampPalette(c("blue","white","deeppink4"))(100), breaks=seq(-3, 3, length.out = 100),
                       treeheight_row=10, treeheight_col=10, border_color=NA, silent=TRUE)

    # saves heatmap as a PDF
    save_pheatmap_pdf(ph, paste0(fileloc_cells, "_regulon_heatmap.v1.pdf"))

    # writes a file of the most active genes in regulons (?)
    topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
    colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
    topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]

    write.table(topRegulators, file=paste0(fileloc_cells, "_top_regulators.v1.txt"),sep="\t", quote=F)

    # threshold for binary matrix is 70% of the cells must have active regulon to be "on"

    # creates binary regulon matrix
    minPerc <- .7
    binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
    cellInfo_binarizedCells <- cellInfo[which(rownames(cellInfo)%in% colnames(binaryRegulonActivity)),, drop=FALSE]
    regulonActivity_byCellType_Binarized <- sapply(split(rownames(cellInfo_binarizedCells), cellInfo_binarizedCells$seuratCluster),
                                                   function(cells) rowMeans(binaryRegulonActivity[,cells, drop=FALSE]))
    binaryActPerc_subset <- regulonActivity_byCellType_Binarized[which(rowSums(regulonActivity_byCellType_Binarized>minPerc)>0),]

    # stores the binary regulon matrix as a heatmap
    pb<-pheatmap::pheatmap(binaryActPerc_subset, # fontsize_row=5,
                       color = colorRampPalette(c("white","pink","deeppink4"))(100), breaks=seq(0, 1, length.out = 100),
                       treeheight_row=10, treeheight_col=10, border_color=NA, silent=TRUE)

    save_pheatmap_pdf(pb, file=paste0(fileloc_cells, "_binary_regulon_heatmap_v1.pdf"))

    # writes a file of the most active genes in the regulon, based on the binarized regulon activity
    topRegulatorsBinary <- reshape2::melt(regulonActivity_byCellType_Binarized)
    colnames(topRegulatorsBinary) <- c("Regulon", "CellType", "RelativeActivity")
    topRegulatorsBinary <- topRegulatorsBinary[which(topRegulatorsBinary$RelativeActivity>minPerc),]

    write.table(topRegulatorsBinary, file=paste0(fileloc_cells, "_top_binary_regulators.v1.txt"), sep="\t", quote=F)

    # save data
    save.image(paste0(fileloc_cells, ".v1.RData"))

}
