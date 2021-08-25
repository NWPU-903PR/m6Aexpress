##peak methy-level distance decay
##caculate all differential methylation peak site normalize dist decay methy level
gene_methy_intensity <- function(peak_inform,txdbinfor, GENE_ANNO_GTF, species="human"){
  ## caculate the site distance to stop codon

##methy level decay
.methy_level <-  function(IP_Input_read,size_factor){
  IP_site_read <- IP_Input_read[,grep("IP",colnames(IP_Input_read))]
  Input_site_read <- IP_Input_read[,(grep("Input",colnames(IP_Input_read)))]
  IP_Input <- cbind(IP_site_read, Input_site_read)
  IP_Input_norm <- as.data.frame(t(t(IP_Input)/size_factor))
  IP_norm_site <- IP_Input_norm[,1:(ncol(IP_Input_norm)/2)]
  Input_norm_site <- IP_Input_norm[,((ncol(IP_Input_norm)/2)+1):(ncol(IP_Input_norm))]
  methy_level <- log2((IP_norm_site+0.01)/(Input_norm_site+0.01))
  methy_level_infor <- cbind(IP_Input_read[,1:6], methy_level )
  for (i in 1:nrow(methy_level_infor)) {
    for(j in 7:ncol(methy_level_infor)){

      if(methy_level_infor[i,j]<0){

        methy_level_infor[i,j] <- 0.001
      }
    }

  }
  return(methy_level_infor)
}
  peak_infor <- peak_inform[[1]]
  total_reads <- peak_inform[[2]]
  size_factor <- as.numeric(total_reads/exp(mean(log(total_reads))))
  norm_methy_level <- .methy_level(peak_infor,size_factor)
  dist_decay <- .dist_stopcodon(norm_methy_level, txdbinfor, GENE_ANNO_GTF, species)
  select_peak_methy <- norm_methy_level[which(!is.na(match(as.character(norm_methy_level$gene_name), as.character(dist_decay$gene_name)))),]
  rownames(select_peak_methy) <- NULL
  peaklevel_decay <- data.frame()
  table_genename <- as.data.frame(table(as.character(select_peak_methy$gene_name)))
  colnames(table_genename) <- c("gene_name","freq")
  select_gene_name <- unique(as.character(select_peak_methy$gene_name))
  for (i in 1:length(select_gene_name)) {
    one_gene <- table_genename[table_genename$gene_name==select_gene_name[i],]
    if(as.numeric(one_gene$freq)==1){
      add_decay_methy <- select_peak_methy[which(!is.na(match(select_peak_methy$gene_name, select_gene_name[i]))),7:ncol(select_peak_methy)]*
        dist_decay[dist_decay$gene_name==select_gene_name[i],-1]
      all_peak_methydecay <- colSums(add_decay_methy)
      peaklevel_decay <- rbind(peaklevel_decay, all_peak_methydecay)
    }
    if(as.numeric(one_gene$freq)>1){
      add_decay_methy <- select_peak_methy[which(!is.na(match(select_peak_methy$gene_name, select_gene_name[i]))),7:ncol(select_peak_methy)]*
        dist_decay[dist_decay$gene_name==select_gene_name[i],-1]
      all_peak_methydecay <- colSums(add_decay_methy)
      peakcount_decay <- unique(colSums(dist_decay[dist_decay$gene_name==select_gene_name[i],-1]))
      norm_peak_methydecay <- round((all_peak_methydecay/peakcount_decay),2)
      peaklevel_decay <- rbind(peaklevel_decay, norm_peak_methydecay)

    }
  }

  new_peakdecaylevel <- cbind(select_gene_name, peaklevel_decay)
  colnames(new_peakdecaylevel) <- c("gene_name", colnames(select_peak_methy)[7:ncol(select_peak_methy)])
  ##select methy-level
  select_peakdecaylevel <- new_peakdecaylevel[rowSums(new_peakdecaylevel[,-1])>0,]
  last_select_peakdecaylevel <- select_peakdecaylevel[rowSums(round(select_peakdecaylevel[,-1],2))>0,]
  last_select_peak <- na.omit(last_select_peakdecaylevel)
  last_gene_name <- as.character(last_select_peak$gene_name)
  last_selectpeakdecaylevel <- cbind(last_gene_name,round(last_select_peak[,-1],2))
  colnames(last_selectpeakdecaylevel)[1] <- "gene_name"
  rownames(last_selectpeakdecaylevel) <- NULL
  peakdecaylevel <- last_selectpeakdecaylevel[,-1]
  meanmethy <- rowMeans(peakdecaylevel)
  select_label <- which(meanmethy>0.1)
  last_peakdecaylevel <- last_selectpeakdecaylevel[select_label, ]
  return(last_peakdecaylevel)
}

