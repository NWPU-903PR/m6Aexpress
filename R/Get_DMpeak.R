## DE methy for each gene
## select DE significant site for each gene and add the log2OR
DM_detect <- function(peak_inform,DM_CUTOFF_TYPE,num_ctl,
                            diff_peak_fdr,diff_peak_pvalue){
  peak_infor <- peak_inform[[1]]
  libsizes <- peak_inform[[2]]
  methy_read <- peak_infor[,grep("IP", colnames(peak_infor))]
  ##update
  CTL_methy_read <- as.data.frame(methy_read[,1:num_ctl])
  Treated_methy_read <- methy_read[,(num_ctl+1):ncol(methy_read)]
  unmethy_read <- peak_infor[,grep("Input", colnames(peak_infor))]
  CTL_unmethy_read <- as.data.frame(unmethy_read[,1:num_ctl])
  Treated_unmethy_read <- unmethy_read[,(num_ctl+1):ncol(unmethy_read)]
  gene_name <- peak_infor$gene_name
  peak_width <- peak_infor$width
  ##get size factor
  read_count <- cbind(methy_read, unmethy_read)
  size_factor=as.numeric(libsizes/exp(mean(log(libsizes))))
  ##QNB test methylation site for each gene
  gene_site_read <- cbind(gene_name,peak_infor$seqnames, peak_infor$start, peak_infor$end,peak_width,peak_infor$strand, CTL_methy_read, Treated_methy_read, CTL_unmethy_read, Treated_unmethy_read )
  colnames(gene_site_read)[2:6] <- c("seqnames", "start", "end", "width", "strand")
  qnb_test <- qnbtest(CTL_methy_read, Treated_methy_read,CTL_unmethy_read,Treated_unmethy_read,mode="per-condition")
  QNB_test <- cbind(gene_site_read,qnb_test[,c(4,5,7)])
  if (DM_CUTOFF_TYPE =="padj") {select_DMpeak_site <- QNB_test[which((QNB_test$padj<diff_peak_fdr)==TRUE),]}
  if (DM_CUTOFF_TYPE =="pvalue") {select_DMpeak_site <- QNB_test[which((QNB_test$pvalue<diff_peak_pvalue)==TRUE),]}
  select_DMpeak <- list(select_DMpeak_site, size_factor)
  DM_peak_infor <- select_DMpeak_site[,c(1:4,6,(ncol(select_DMpeak_site)-1):(ncol(select_DMpeak_site)-2))]
  names(select_DMpeak) <- c("DM_peak_infor", "library_sizefactor")
  return(select_DMpeak)
}
