obtain_HVP_sites <- function(peak_inform,CV_values=0.3,num_sample_subgroup){
  ## caculate the site distance to stop codon
  ##Obtain High Variable Peak
.methy_level <-  function(IP_Input_read,size_factor){
  ##get methy level
  IP_site_read <- IP_Input_read[,grep("IP",colnames(IP_Input_read))]
  Input_site_read <- IP_Input_read[,(grep("Input",colnames(IP_Input_read)))]
  IP_Input <- cbind(IP_site_read, Input_site_read)
  size_factor=size_factor
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
  size_factor <- peak_inform[[2]]
  norm_methy_level <- .methy_level(IP_Input_read=peak_infor,size_factor=size_factor)
  methy_level<- norm_methy_level[,-c(1:6)]
  subgroup_methy <- data.frame()
  for (i in 1:length(num_sample_subgroup)) {
    j <- num_sample_subgroup[i]
    one_group_methy <- as.numeric(as.character(rowMeans(as.data.frame(methy_level[,1:j]))))
    subgroup_methy <- as.data.frame(rbind(subgroup_methy,one_group_methy))
    if(i<length(num_sample_subgroup)){
      methy_level <- as.data.frame(methy_level[,-c(1:j)])
    }
    colnames(subgroup_methy) <- NULL
  }
  subgroup_methy <- as.data.frame(t(subgroup_methy))
  colnames(subgroup_methy) <- paste0("subgroup_",1:length(num_sample_subgroup))
  methy_mean<- round(as.numeric(apply(subgroup_methy, 1, mean)),4) 
  methy_sd <- round(as.numeric(apply(subgroup_methy,1,sd)),4)
  cv_value<- methy_sd/methy_mean
  variable_peak_lable<-which(cv_value>CV_values)
  High_variable_peak <- peak_infor[variable_peak_lable,]
  high_variable_peak_infor <- list(High_variable_peak,size_factor)
  names(high_variable_peak_infor) <- c("high_variable_peak","size_factor")
  return(high_variable_peak_infor)
}
