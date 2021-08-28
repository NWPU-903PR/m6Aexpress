add_LFC_DDM <- function(expre_methyre, DE_gene, methy_intensity,num_cond1,OUTPUT_DIR=NA){
  DE_infor <- DE_gene[[3]]
  DE_infor <- data.frame(gene_name=as.character(rownames(DE_infor)),DE_infor)
  rownames(DE_infor) <- NULL
  control_methyintensity <- as.data.frame(methy_intensity[,2:(num_cond1+1)]) 
  treated_mehtyintensity <- as.data.frame(methy_intensity[,-c(1,2:(num_cond1+1))])
  DDM <- rowMeans(treated_mehtyintensity)-rowMeans(control_methyintensity)
  add_methyDDM <- data.frame(methy_intensity, DDM) 
  m6A_reg_expr_gene <- as.character(expre_methyre$gene_name)
  m6Aexpree_gene_DDMLFC <- data.frame()
  for (i in 1:length(m6A_reg_expr_gene)) {
    one_DEDM <- as.data.frame(cbind(expre_methyre[i,],DE_infor[DE_infor$gene_name==m6A_reg_expr_gene[i],]$log2FoldChange,
                                    add_methyDDM[add_methyDDM$gene_name==m6A_reg_expr_gene[i],ncol(add_methyDDM)]))
    colnames(one_DEDM)[(ncol(one_DEDM)-1):ncol(one_DEDM)] <- c("DE_log2foldchange","DM_intensity")
    m6Aexpree_gene_DDMLFC <- rbind(m6Aexpree_gene_DDMLFC,one_DEDM)
  }
  return(m6Aexpree_gene_DDMLFC)
  if(is.na(OUTPUT_DIR)){
    OUTPUT_DIR <- getwd()
    write.table(m6Aexpree_gene_DDMLFC,file=paste(OUTPUT_DIR,"m6Aexpress_result", "m6A_reg_exp_LFC_DDM.xls",sep="/"), sep="\t",row.names =FALSE,quote = FALSE)
  }
  if(!is.na(OUTPUT_DIR)){
    write.table(m6Aexpree_gene_DDMLFC,file=paste(OUTPUT_DIR,sep="/", "m6A_reg_exp_LFC_DDM.xls"), sep="\t",row.names =FALSE,quote = FALSE)
  }
}
