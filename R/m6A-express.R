m6A_Express_model <- function(Input_file,CUTOFF_TYPE,pvalue, FDR,out_dir=NA){
  py_code <- system.file("extdata", "R_runpython.py", package = "m6Aexpress")
  source_python(py_code)
  fileNameCount=Input_file[[1]]
  librarySizes=as.numeric(Input_file[[2]])
  out_result <- suppressMessages(try(fun_R_call(fileNameCount, librarySizes),silent=TRUE))
  j=0
  while((is.null(nrow(out_result)))&j<10){
    out_result <- suppressMessages(try(fun_R_call(fileNameCount, librarySizes),silent=TRUE))
    j <- j+1
  }
  if(j==10){
    print("There is no significant m6A regulated expression gene")
    }
  if(j<10){
    genecoutmethy <-read.table(fileNameCount,header = T)
  size_factor<-librarySizes
  exprmethyre <- out_result
  match_count_methy <- data.frame()
  for (i in 1:length(exprmethyre$Gene_ID)) {
    one_gene <- genecoutmethy[genecoutmethy$gene_name==as.character(exprmethyre$Gene_ID)[i],]
    match_count_methy <- rbind(match_count_methy,one_gene)
  }

  match_methy <- data.frame()
  for (i in 1:nrow(match_count_methy)) {
    one_methy <- exprmethyre[exprmethyre$Gene_ID==as.character(match_count_methy$gene_name)[i],]
    match_methy <- rbind(match_methy, one_methy)
  }

  genecount <- match_count_methy[,2:((ncol(match_count_methy)+1)/2)]
  methyintensity <- match_count_methy[,(((ncol(match_count_methy)+1)/2)+1):ncol(match_count_methy)]

  alpha <- vector()
  beta <- vector()
  intersecpt <- vector()
  for (i in 1:nrow(match_count_methy)) {
    alpha[i] <- as.numeric(exprmethyre[exprmethyre$Gene_ID==as.character(match_count_methy$gene_name)[i],ncol(exprmethyre)])
    beta[i] <-  as.numeric(exprmethyre[exprmethyre$Gene_ID==as.character(match_count_methy$gene_name)[i],3])
    intersecpt[i] <-as.numeric(exprmethyre[exprmethyre$Gene_ID==as.character(match_count_methy$gene_name)[i],2])
  }

  beta_value <- cbind(intersecpt,beta)
  base_level <- rep(1, ncol(methyintensity))
  gene_name <- as.character(match_count_methy$gene_name)

  ##empirical prior estimate for beta
  ##use bet0 and bet1 prior with mu
  beta0_lower <- quantile(intersecpt,0.1)
  beta0_up <- quantile(intersecpt,0.9)
  beta1_lower <- quantile(beta,0.1)
  beta1_up <- quantile(beta,0.9)
  select_beta0 <- intersecpt[which((intersecpt>beta0_lower)&(intersecpt<beta0_up))]
  select_beta1 <- beta[which((beta>beta1_lower)&(beta<beta1_up))]
  nu=mean(select_beta1)
  iu=mean(select_beta0)
  var_beta<-var(select_beta1)
  var_IN <- var(select_beta0)
  lambda <- (1/var_beta)
  lambda_IN <- (1/var_IN)
  ##final estimate the beta value
  new_beta <- data.frame()
  for (i in 1:nrow(methyintensity)) {
    methy_gene <- t(methyintensity[i,])
    gene_matrix <- cbind( base_level,methy_gene)
    model_matrix <- matrix(cbind(gene_matrix[,1],gene_matrix[,2]),ncol = ncol(gene_matrix))
    mu <- as.numeric(size_factor*t(exp(gene_matrix%*% (beta_value[i,]))))
    gene_count <- as.numeric(genecount[i,])
    disp <- alpha[i]

    beta_mat <- matrix(c(beta_value[i,1],beta_value[i,2]),ncol = ncol(beta_value))
    counts <- matrix(gene_count,nrow = length(gene_count))
    matrix_size_factor <-matrix(size_factor, ncol = length(size_factor))


    m = nrow(model_matrix)
    ridge <- diag(lambda, nrow = ncol(model_matrix), ncol = ncol(model_matrix))
    beta_hat = t(beta_mat)
    ##get loglike hood function of beta
    beta_loglike <- function(beta) {
      mu_row <- as.numeric(t(matrix_size_factor)*exp(model_matrix%*%beta))
      count_row <- as.numeric(counts)
      nbinomLogLike <- sum(dnbinom(count_row,mu=mu_row,size=1/disp,log=TRUE))
      logPrior <- sum(dnorm(beta,mean=c(iu,nu),sd=c(sqrt(1/lambda_IN), sqrt(1/lambda)),log=TRUE))
      beta_LogLike <- -1*(nbinomLogLike+logPrior)

    }
    o <- optim(beta_hat, beta_loglike, method="L-BFGS-B",hessian = TRUE)
    one_new <- t(o$par)
    se.beta.hat <- sqrt(diag(solve(o$hessian)))
    ##wald test

    sigma <- solve(o$hessian)

    wdtest <- wald.test(b =as.numeric(one_new) , Sigma = sigma, Terms = 2)
    pvalues <- as.numeric(wdtest$result$chi2[3])
    one_data <- t(c(as.numeric(one_new),se.beta.hat, pvalues))
    new_beta <- rbind(new_beta, one_data)
  }
  ##get loglike hood function of beta
  gene_name <-as.character(match_methy$Gene_ID)
  adj_beta <- cbind(gene_name, new_beta)
  colnames(adj_beta) <- c("gene_name", "Beta0","Beta1","SE_Beta0","SE_Beat1", "pvalue")
  pvalues <- as.numeric(adj_beta$pvalue)
  padj <- p.adjust(pvalues, method = "BH")
  padj_beta <- cbind(adj_beta, padj)
  if(is.na(out_dir)){
    out_dir = dirname(Input_file[[1]])
  }
   if (CUTOFF_TYPE =="FDR") {
     select_adjbeta <- padj_beta[padj_beta$padj<FDR,]
     write.table(select_adjbeta,file=paste0(out_dir,sep="/","m6A-express_result.xls"), sep="\t",row.names =FALSE,quote = FALSE)

    }
  if (CUTOFF_TYPE =="pvalue") {
    select_adjbeta <- padj_beta[padj_beta$pvalue<pvalue,]
    write.table(select_adjbeta,file=paste0(out_dir,sep="/", "m6A-express_result.xls"), sep="\t",row.names =FALSE,quote = FALSE)

    }
   return(select_adjbeta)
  }

}

