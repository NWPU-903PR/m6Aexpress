Get_express_data <- function(INPUT_BAM, TREATED_INPUT_BAM=character(0),annot_file="hg19",
                             species="human",isPairedEnd=FALSE,GENE_ANNO_GTF = NULL,isGTFAnnotationFile=FALSE){
  Input_data <- c(INPUT_BAM, TREATED_INPUT_BAM)
  gene_count <- featureCounts(Input_data,useMetaFeatures=TRUE, annot.inbuilt=annot_file, 
                              isPairedEnd=isPairedEnd, annot.ext = GENE_ANNO_GTF,isGTFAnnotationFile=isGTFAnnotationFile)
  counts_data <- gene_count$counts
  countMatrix <- sapply(as.matrix(counts_data), as.numeric)
  gene_countmatrix <- matrix(countMatrix, nrow=nrow(counts_data), ncol = ncol(counts_data))
  colnames(gene_countmatrix) <- paste0("sample_",1:length(Input_data))
  rownames(gene_countmatrix) <- rownames(counts_data)
  conds <- factor(colnames(counts_data))
  cds <- newCountDataSet(gene_countmatrix, conds )
  size_factor <-  sizeFactors(estimateSizeFactors( cds ))
  gene_ID <- rownames(gene_countmatrix)
  if(isGTFAnnotationFile==FALSE){
    if(species=="human"){
      org_db <- org.Hs.eg.db
      tans_name <- select(org.Hs.eg.db, keys=gene_ID, columns = c("SYMBOL"),keytype= "ENTREZID")
    }
    if(species=="mouse"){
      org_db <- org.Mm.eg.db
      tans_name <- select(org_db, keys=gene_ID, columns = c("SYMBOL"),keytype= "ENTREZID")
    }
    if(species=="yeast"){
      org_db <- org.Sc.sgd.db
      KeyType="ORF"
      tans_name <- select(org_db, keys=gene_ID, columns = c("GENENAME"),keytype= "ORF")
    }
  gene_name <- as.character(tans_name[,2])
  }

  gene_name <- gene_ID
  rownames(gene_countmatrix)  <- gene_name
  gene_countdata <- as.data.frame(gene_countmatrix)
  gene_countdata <- na.omit(gene_countdata)
  Gene_count_infor <- list(gene_countdata, size_factor)
  return(Gene_count_infor)
}
