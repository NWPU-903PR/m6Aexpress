m6Aexpress <- function(express_data, treated_express_data=character(0),
                        IP_BAM, INPUT_BAM,TREATED_IP_BAM=character(0),TREATED_INPUT_BAM=character(0),
                        annot_type="hg19",species="human",isPairedEnd=FALSE,
                        isGTFAnnotationFile=FALSE,
                        GENE_ANNO_GTF=NULL,
                        GENOME = NA, UCSC_TABLE_NAME = "knownGene", model="basic",
                        pvalue=NA,
                        FDR=0.05,
                        CV_values=0.3,
                        num_sample_subgroup,
                        diff_gene_pvalue=NA,
                        diff_gene_fdr=0.05,
                        diff_peak_pvalue=NA,
                        diff_peak_fdr=0.05,
                        TXDB=NA, OUTPUT_DIR= NA){
  ##get gene read count
  print("Get reads count for each gene. It may spend some minutes.")
  gene_express_data <- Get_express_data(INPUT_BAM=express_data, TREATED_INPUT_BAM=treated_express_data,
                                        annot_file=annot_type,species=species,isPairedEnd=isPairedEnd,
                                        GENE_ANNO_GTF=GENE_ANNO_GTF,isGTFAnnotationFile=isGTFAnnotationFile)
  ##get peak sites infor
  print("Do peak calling. It may spend some hours.")
  if (is.na(OUTPUT_DIR)) {OUTPUT_DIR=getwd()} else {OUTPUT_DIR=OUTPUT_DIR}
  get_peak_site <- Get_peakinfor(IP_BAM=IP_BAM, INPUT_BAM=INPUT_BAM,TREATED_IP_BAM=TREATED_IP_BAM,
                                 TREATED_INPUT_BAM=TREATED_INPUT_BAM, species=species,
                                 GENOME = GENOME, UCSC_TABLE_NAME = UCSC_TABLE_NAME, GENE_ANNO_GTF=GENE_ANNO_GTF, TXDB=TXDB, OUTPUT_DIR=OUTPUT_DIR)
  if(model=="basic"){
    ##methylation intensity
    print("Calculate methylation intensity for each gene. It may spend some minutes.")
    gene_methyintensity <- gene_methy_intensity(peak_inform=get_peak_site,txdbinfor=TXDB,GENE_ANNO_GTF=GENE_ANNO_GTF, species=species)

    ##match expression and methylation intensity
    print("Obtain methylation regulated expression gene by m6A-express model")
    paired_expr_methy <- match_expr_methy(gene_count_infor=gene_express_data, gene_methy_infor=gene_methyintensity,OUTPUT_DIR=OUTPUT_DIR)
    if (is.na(pvalue)) {CUTOFF_TYPE="padj"} else  {CUTOFF_TYPE="pvalue"}
    m6Areg_expr_gene <- m6A_Express_model(Input_file=paired_expr_methy,CUTOFF_TYPE=CUTOFF_TYPE,pvalue=pvalue,
                                          FDR=FDR, out_dir=OUTPUT_DIR)
    return(m6Areg_expr_gene)
  }
    if(model=="HVP"){
     ###obtain the high variable Peak sites
     get_HVP_peak <- obtain_HVP_sites(peak_inform=get_peak_site,CV_values=0.3,num_sample_subgroup=num_sample_subgroup)
    ##methylation intensity
    print("Calculate methylation intensity for each gene with high variable peak. It may spend some minutes.")
    gene_methyintensity <- gene_methy_intensity(peak_inform=get_HVP_peak,txdbinfor=TXDB,GENE_ANNO_GTF=GENE_ANNO_GTF, species=species)

    ##match expression and methylation intensity
    print("Obtain methylation regulated expression gene by m6A-express model")
    paired_expr_methy <- match_expr_methy(gene_count_infor=gene_express_data, gene_methy_infor=gene_methyintensity,OUTPUT_DIR=OUTPUT_DIR)
    if (is.na(pvalue)) {CUTOFF_TYPE="padj"} else  {CUTOFF_TYPE="pvalue"}
    m6Areg_expr_gene <- m6A_Express_model(Input_file=paired_expr_methy,CUTOFF_TYPE=CUTOFF_TYPE,pvalue=pvalue,
                                          FDR=FDR, out_dir=OUTPUT_DIR)
    return(m6Areg_expr_gene)
  }

  if(model=="DE-DM"){
    ##get DE gene
    print("Obtain differential expression gene.")
    num_cond1 <- length(express_data)
    num_cond2 <- length(treated_express_data)
    cond1 <- "control"
    cond2 <- "treated"
    if (is.na(diff_gene_pvalue)) {DE_CUTOFF_TYPE="padj"} else  {DE_CUTOFF_TYPE="pvalue"}
    DE_gene <- Select_DEgene(gene_count_infor=gene_express_data,cond1=cond1,cond2=cond2,
                             num_cond1=num_cond1, num_cond2=num_cond2,DE_CUTOFF_TYPE=DE_CUTOFF_TYPE,
                             DIFF_GENE_cutoff_FDR=diff_gene_fdr,DIFF_GENE_CUTOFF_PVALUE=diff_gene_pvalue)
    ##get DM peak site
    print("Obtain differential methylation peak sites")
    num_ctl <- length(IP_BAM)
    if (is.na(diff_peak_pvalue)) {DM_CUTOFF_TYPE="padj"} else  {DM_CUTOFF_TYPE="pvalue"}
    DM_methy <- DM_detect(peak_inform=get_peak_site,DM_CUTOFF_TYPE=DM_CUTOFF_TYPE,num_ctl=num_ctl,
                          diff_peak_fdr=diff_peak_fdr,diff_peak_pvalue=diff_peak_pvalue)
    print("Calculate methylation intensity for differential methylation gene. It may spend some minutes.")
    gene_methyintensity <- gene_methy_intensity(peak_inform=DM_methy,txdbinfor=TXDB,GENE_ANNO_GTF=GENE_ANNO_GTF, species=species)

    ##match expression and methylation intensity
    print("Identify differential expression gene are regulated by differential methylation peak.")
    paired_expr_methy <- match_expr_methy(gene_count_infor=DE_gene, gene_methy_infor=gene_methyintensity,OUTPUT_DIR=OUTPUT_DIR)
    if (is.na(pvalue)) {CUTOFF_TYPE="padj"} else  {CUTOFF_TYPE="pvalue"}
    m6Areg_expr_gene <- m6A_express_model(Input_file=paired_expr_methy,CUTOFF_TYPE=CUTOFF_TYPE,pvalue=pvalue,
                                          FDR=FDR, out_dir=OUTPUT_DIR)

    ##add log2foldchange and difference degree methylation
    m6A_express_addLFC_DDM <- add_LFC_DDM(expre_methyre=m6Areg_expr_gene, DE_gene=DE_gene, methy_distdecay=DM_methy,OUTPUT_DIR=OUTPUT_DIR)
    return(m6A_express_addLFC_DDM)
  }
}
