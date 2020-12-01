##Get peak site and peak information for each gene
Get_peakinfor <- function(IP_BAM, INPUT_BAM,TREATED_IP_BAM=character(0),
                          TREATED_INPUT_BAM=character(0), species="human",
                          GENOME = NA, UCSC_TABLE_NAME = "knownGene", GENE_ANNO_GTF=NULL, TXDB=NA, OUTPUT_DIR= NA){
  IP_bam <- c(IP_BAM, TREATED_IP_BAM)
  INPUT_bam <- c(INPUT_BAM, TREATED_INPUT_BAM)
  # Get the annotation file
  if (suppressWarnings((!is.na(GENOME)) & (!is.na(UCSC_TABLE_NAME)) &
                       is.na(TXDB) & is.null(GENE_ANNO_GTF))) {
    op <- options(warn = (-1))
    txdb = makeTxDbFromUCSC(genome = GENOME, tablename = UCSC_TABLE_NAME)
    KeyType="ENTREZID"
    options(op)
  }
  if (suppressWarnings(!is.null(GENE_ANNO_GTF) & is.na(TXDB))) {
    op <- options(warn = (-1))
    txdb <- makeTxDbFromGFF(GENE_ANNO_GTF, format = "gtf")
    options(op)
  }

  # use provided annotation data file
  if (suppressWarnings(!is.na(TXDB))) {
    txdb <- loadDb(TXDB)
    KeyType="ENTREZID"
  }

  ##Get peak site from bam file
  if(is.na(OUTPUT_DIR)){
    OUTPUT_DIR=getwd()
  }
  result <- exomepeak(IP_BAM=IP_bam, INPUT_BAM=INPUT_bam,GENE_ANNO_GTF=GENE_ANNO_GTF,
                      GENOME =  GENOME, UCSC_TABLE_NAME = "knownGene",TXDB=TXDB,
                      OUTPUT_DIR=OUTPUT_DIR)
  load(paste(OUTPUT_DIR,"exomePeak_output", "exomePeak.Rdata",sep="/"))
  peak_data <- tmp_rs
  consisten_peak <- paste(OUTPUT_DIR,"exomePeak_output","con_peak.bed",sep="/")
  read_peak <- import.bed(consisten_peak)
  read_peak <- as.data.frame(read_peak)
  peak_name <- as.character(read_peak$name)
  if(is.null(GENE_ANNO_GTF)){
    if(species=="human"){
      org_db <- org.Hs.eg.db
      tans_name <- select(org.Hs.eg.db, keys=peak_name, columns = c("SYMBOL"),keytype= KeyType)
    }
    if(species=="mouse"){
      org_db <- org.Mm.eg.db
      tans_name <- select(org_db, keys=peak_name, columns = c("SYMBOL"),keytype= KeyType)
    }
    if(species=="yeast"){
      org_db <- org.Sc.sgd.db
      KeyType="ORF"
      tans_name <- select(org_db, keys=peak_name, columns = c("GENENAME"),keytype= KeyType)
    }
    select_peak <-cbind( as.character(read_peak$seqnames) , read_peak$start, read_peak$end,read_peak$width, as.character(read_peak$strand), as.character(tans_name[,2]))
  }
  select_peak <-cbind( as.character(read_peak$seqnames) , read_peak$start, read_peak$end,read_peak$width, as.character(read_peak$strand),as.character(read_peak$name))
  colnames(select_peak) <- c("seqnames", "start", "end", "width", "strand", "gene_name")
  select_peak <- as.data.frame(select_peak)
  peak <- peak_data[["PEAK"]]
  READS_COUNT <- peak_data[["READS_COUNT"]]
  ## get size factor
  reads_count <- READS_COUNT[,-((ncol(READS_COUNT)-1):ncol(READS_COUNT))]
  totalreads <- colSums(reads_count)
  # get number of consistent peaks
  no_peak=length(peak$loci2peak_consistent[,1])
  # peak_reads_count
  peak_reads_count = READS_COUNT[1:no_peak,]
  no_sample=length(peak_reads_count[1,])-2
  # cut the unnecessary information
  peak_reads_count = READS_COUNT[1:no_peak,1:no_sample]
  # count
  for (ipeak in 1:no_peak) {
    temp=peak$loci2peak_consistent[ipeak,]
    temp2=colSums(READS_COUNT[temp[1]:temp[2],1:no_sample])
    peak_reads_count[ipeak,1:no_sample]=temp2
  }
  # remove the overlapping window effects
  peak_reads_count = round (peak_reads_count * 30 / 200);
  if(length(TREATED_INPUT_BAM)==0){
    colnames(peak_reads_count) <- c(paste0("IP",1:length(IP_BAM)),paste0("Input",1:length(INPUT_BAM)))
  }
  if(length(TREATED_INPUT_BAM)!=0){
    colnames(peak_reads_count) <- c(paste0("IP",1:length(IP_BAM)),paste0("Treated_IP",1:length(TREATED_IP_BAM)),
                                    paste0("Input",1:length(INPUT_BAM)), paste0("Treated_Input",1:length(TREATED_INPUT_BAM)))
  }
  ## get fdr
  log_fdr <- peak$PW$log.fdr
  log_fc <- peak$PW$log.fc
  consisten_peak <- peak$Consistent
  consit_log_fdr <- log_fdr[consisten_peak]


  # initializa the peak reporting
  no_peak=length(peak$loci2peak_consistent[,1])
  # get peak
  peak_report <- data.frame()
  for (i in 1:no_peak) {
    peak_row_id=peak$loci2peak_consistent[i,]

    # batch id
    batch_id=unique(READS_COUNT$batch_id[peak_row_id])
    lg.p=min(peak$PW$log.p[peak_row_id[1]:peak_row_id[2]])/log(10)
    lg.fdr=min(peak$PW$log.fdr[peak_row_id[1]:peak_row_id[2]])/log(10)
    fold_enrchment=exp(max(peak$PW$log.fc[peak_row_id[1]:peak_row_id[2]]))

    # get sig digits
    lg.p=signif(lg.p, digits = 3)
    lg.fdr=signif(lg.fdr, digits = 3)
    fold_enrchment=signif(fold_enrchment, digits = 3)

    test_result <- data.frame(lg.p=lg.p,lg.fdr=lg.fdr, fold_enrchment=fold_enrchment)
    peak_report=rbind(peak_report,test_result)
  }
  ## peak site read and log(fdr)
  peak_site_reads <- cbind(select_peak, peak_reads_count, peak_report[,2:3])
  peak_site_reads <- peak_site_reads[which(!is.na(peak_site_reads$gene_name)), ]
  peak_site_infor <- list(peak_site_reads, totalreads)
  return(peak_site_infor)

}
