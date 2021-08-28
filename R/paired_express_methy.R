##match gene count and methylation intensity
match_expr_methy <- function(gene_count_infor, gene_methy_infor,OUTPUT_DIR=NA){
  ##select DE gene data
  gene_count <- gene_count_infor[[1]]
  gene_count <- gene_count[which(rowSums(gene_count)>10),]
  for (i in 1:nrow(gene_count)) {
    for (j in 2:ncol(gene_count)) {
      if(gene_count[i,j]==0){

        gene_count[i,j] <- NA
      }

    }

  }
  select_genecount <- na.omit(gene_count)

  size_factor <- gene_count_infor[[2]]
  gene_name <- (as.character(rownames(select_genecount)))
  select_genecount <- as.data.frame(cbind(gene_name,select_genecount))
  rownames(select_genecount) <- NULL
  intersect_gene <- intersect(select_genecount$gene_name, gene_methy_infor$gene_name)
  match_data <- data.frame()
  for (i in 1:length(intersect_gene)) {

    methy_name <- intersect_gene[i]
    match_gene <- as.data.frame(select_genecount[which(!is.na(match(select_genecount$gene_name,  methy_name))),])
    select_methy <- as.data.frame(gene_methy_infor[which(!is.na(match(gene_methy_infor$gene_name, methy_name))),])
    match_methy_expr <- cbind(match_gene, select_methy[,-1])
    match_data <- rbind(match_data, match_methy_expr)

  }
  if(is.na(OUTPUT_DIR)){
    OUTPUT_DIR=getwd()
  }
  dir.create(paste0(OUTPUT_DIR,"/m6Aexpress_result"))
  write.table(match_data, file = paste0(OUTPUT_DIR,"/m6Aexpress_result/","expr_methy.tab"), quote = FALSE, row.names = FALSE)
  out_putfile <- paste0(OUTPUT_DIR,"/m6Aexpress_result/","expr_methy.tab")
  m6A_express_input <- list(out_putfile,size_factor)
  names(m6A_express_input) <- c("gene_express_methy_dir","size_factor")
  return(m6A_express_input)
}
