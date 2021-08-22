##stop codon distance decay
.dist_stopcodon <- function(methy_level_infor, txdbinfor,GENE_ANNO_GTF,species){
  ##construct GRange object
  seqname <- as.character(methy_level_infor$seqnames)
  start_site <- as.numeric(as.character(methy_level_infor$start))
  end_site <- as.numeric(as.character(methy_level_infor$end))
  strand_site <- as.character(methy_level_infor$strand)
  gene_name <- as.character(methy_level_infor$gene_name)
  gene_methy_GR <- GRanges(seqnames = seqname, ranges = IRanges(start_site, end_site, names =gene_name)
                           , strand =strand_site)
  if(is.null(GENE_ANNO_GTF)){
    txdb <- loadDb(txdbinfor)
    CDS <- cdsBy(txdb, by="tx",use.names=TRUE)
    cds_name <- names(CDS)
    if(species=="human"){
      org_db <- org.Hs.eg.db
      tans_cdsname <- select(org.Hs.eg.db, keys=cds_name, columns = c("SYMBOL"),keytype="UCSCKG")
      select_cds <- tans_cdsname[!is.na(tans_cdsname$SYMBOL),]
      select_cdsname <- as.character(select_cds$UCSCKG)
      CDS_select <- CDS[select_cdsname]
      stop_codon <- resize(CDS_select,3,fix = "end")
      ##stop codon region
      stopcodon_trans <- select(org_db, keys=names(stop_codon), columns = c("SYMBOL"),keytype=("UCSCKG"))
      stopcodon_trans <- na.omit(stopcodon_trans)
      
      gene_nameselect <- intersect(unique(as.character(stopcodon_trans$SYMBOL)), unique(names(gene_methy_GR)))
      select_stopcodon_trans <- stopcodon_trans[which(!is.na(match(stopcodon_trans$SYMBOL,gene_nameselect))),]
      
      select_stopcodon <- stop_codon[as.character(select_stopcodon_trans$UCSCKG)]
      stopcodon_select <- select(org_db, keys=names(select_stopcodon), columns = c("SYMBOL"),keytype=("UCSCKG"))
      names(select_stopcodon) <- as.character(stopcodon_select$SYMBOL)
      select_gene_methy_GR <- gene_methy_GR[which(!is.na(match(names(gene_methy_GR),gene_nameselect)))]
    }
    if(species=="mouse"){
      org_db <- org.Mm.eg.db
      tans_cdsname <- select(org_db, keys=cds_name, columns = c("SYMBOL"),keytype="UCSCKG")
      select_cds <- tans_cdsname[!is.na(tans_cdsname$SYMBOL),]
      select_cdsname <- as.character(select_cds$UCSCKG)
      CDS_select <- CDS[select_cdsname]
      stop_codon <- resize(CDS_select,3,fix = "end")
      ##stop codon region
      stopcodon_trans <- select(org_db, keys=names(stop_codon), columns = c("SYMBOL"),keytype=("UCSCKG"))
      stopcodon_trans <- na.omit(stopcodon_trans)
      
      gene_nameselect <- intersect(unique(as.character(stopcodon_trans$SYMBOL)), unique(names(gene_methy_GR)))
      select_stopcodon_trans <- stopcodon_trans[which(!is.na(match(stopcodon_trans$SYMBOL,gene_nameselect))),]
      
      select_stopcodon <- stop_codon[as.character(select_stopcodon_trans$UCSCKG)]
      stopcodon_select <- select(org_db, keys=names(select_stopcodon), columns = c("SYMBOL"),keytype=("UCSCKG"))
      names(select_stopcodon) <- as.character(stopcodon_select$SYMBOL)
      select_gene_methy_GR <- gene_methy_GR[which(!is.na(match(names(gene_methy_GR),gene_nameselect)))]
      
    }
    if(species=="yeast"){
      org_db <- org.Sc.sgd.db
      tans_cdsname <- select(org_db, keys=cds_name, columns = "GENENAME",keytype="ORF")
      select_cds <- tans_cdsname[!is.na(tans_cdsname$GENENAME),]
      select_cdsname <- as.character(select_cds$ORF)
      CDS_select <- CDS[select_cdsname]
      stop_codon <- resize(CDS_select,3,fix = "end")
      ##stop codon region
      stopcodon_trans <- select(org_db, keys=names(stop_codon), columns = c("GENENAME"),keytype=("ORF"))
      stopcodon_trans <- na.omit(stopcodon_trans)
      gene_nameselect <- intersect(unique(as.character(stopcodon_trans$GENENAME)), unique(names(gene_methy_GR)))
      select_stopcodon_trans <- stopcodon_trans[which(!is.na(match(stopcodon_trans$GENENAME,gene_nameselect))),]
      
      select_stopcodon <- stop_codon[as.character(select_stopcodon_trans$ORF)]
      stopcodon_select <- select(org_db, keys=names(select_stopcodon), columns = c("GENENAME"),keytype=("ORF"))
      names(select_stopcodon) <- as.character(stopcodon_select$GENENAME)
      select_gene_methy_GR <- gene_methy_GR[which(!is.na(match(names(gene_methy_GR),gene_nameselect)))]
    }

  }

  if (suppressWarnings(!is.null(GENE_ANNO_GTF) & is.na(txdbinfor))) {
    op <- options(warn = (-1))
    txdb <- makeTxDbFromGFF(GENE_ANNO_GTF, format = "gtf")
    options(op)
    CDS <- cdsBy(txdb, by="gene",use.names=FALSE)
    stop_codon <- resize(CDS,3,fix = "end")
    select_gene_methy_GR <- gene_methy_GR[which(!is.na(match(names(gene_methy_GR),names(stop_codon))))]
    select_stopcodon <- stop_codon[which(!is.na(match(names(stop_codon),names(gene_methy_GR))))]
  }



  
  ## caculate the distance to stop codon
  ## caculate the distance
  dist_peaksite_stop <- vector()
  for (i in 1:length(select_gene_methy_GR)) {
    ##one stopcodon
    one_STD_GR <- select_stopcodon[which(!is.na(match(names(select_stopcodon),names(select_gene_methy_GR)[i])))]
    if(length(one_STD_GR)==1){
      if(as.character(strand(unlist(one_STD_GR)))[1]=="+"){
        STDMid_site <- max(start(one_STD_GR))+1
      }
      if(as.character(strand(unlist(one_STD_GR)))[1]=="-"){
        STDMid_site <- min(start(one_STD_GR))+1
      }
    }
    if(length(one_STD_GR)>1){
      STDMid_site <- vector()
      for (j in 1:length(one_STD_GR)) {
        if(as.character(strand(unlist(one_STD_GR[j])))[1]=="+"){
          STDMid_site[j] <- max(start(one_STD_GR[j]))+1
        }
        if(as.character(strand(unlist(one_STD_GR[j])))[1]=="-"){
          STDMid_site[j] <- min(start(one_STD_GR[j]))+1
        }
      }
    }
    ##one DE methy
    one_gene_GR <- select_gene_methy_GR[i]
    gene_start_site <- rep(start(one_gene_GR), length(STDMid_site))
    gene_end_site <- rep(end(one_gene_GR), length(STDMid_site))
    gene_center_site <- (gene_start_site+rep(((width(one_gene_GR)-1)/2),length(STDMid_site)))
    ##caculate the mine distance between two center
    dist_peaksite_stop[i] <- min(round(abs(gene_center_site - STDMid_site)))
  }
  peak_dist_STDM<- dist_peaksite_stop
  ##select the 75% quantile distance
  d0 <- round(quantile( peak_dist_STDM,0.75))
  ##caculate the peak site decay coefficient
  decay_dist <- data.frame()
  for (j in 1:length(peak_dist_STDM)) {
    decay_fun <- rep(round(exp(-(peak_dist_STDM[j]/d0)),4),(ncol(methy_level_infor)-6))
    decay_dist <- rbind(decay_dist, decay_fun)
  }
  gene_name <- as.character(names(select_gene_methy_GR))
  dist_decay_coeff <- cbind(gene_name, decay_dist)
  return(dist_decay_coeff)
}



