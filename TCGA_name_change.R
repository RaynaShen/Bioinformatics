##author: Rayna
##time: 20180705
##last update: 20180705
##descript: 不同文件之间样品名称的改写对应，等TCGA内需要用到的R function

TMB_info_trans<-function(TMB_info.m,cap=TRUE,pattern="-",length=4){ ##TCGA-02-0047-01
    tmp = c(4,7,12,15)
    if (cap==TRUE){
        changed_name <- sapply(colnames(mut_sum.m),function(x) {substring(x,1,tmp[length])})}
    else {
        changed_name <- sapply(colnames(mut_sum.m),function(x) {tolower(substring(x,1,tmp[length]))})}
    if (pattern == "."){
        changed_name <- sapply(changed_name,function(x) gsub("-",".",x,fixed=TRUE))}
    return(changed_name)
}

cli_info_trans<-function(cli_info.m,cap=FALSE,pattern="-",length=3){ ##tcga-05-4244
    tmp = c(4,7,12)
    if (cap==TRUE){
        changed_name <- sapply(colnames(mut_sum.m),function(x) {toupper(substring(x,1,tmp[length]))})}
    else {
        changed_name <- sapply(colnames(mut_sum.m),function(x) {substring(x,1,tmp[length])})}
    if (pattern == "."){
        changed_name <- sapply(changed_name,function(x) gsub("-",".",x,fixed=TRUE))}
    return(changed_name)
}

mut_sum_trans<-function(mut_sum.m,cap=TRUE,pattern=".",length=4,alpha=TRUE){ ##TCGA.05.4249.01A.01D.1105.08
    ##alpha(TRUE) for 01A, cap(TRUE) for capital TCGA
    tmp = c(4,7,12,15)
    if (cap==TRUE){
        if(alpha==TRUE){
            changed_name <- sapply(colnames(mut_sum.m),function(x) {substring(x,1,tmp[length]+1)})}
        else{changed_name <- sapply(colnames(mut_sum.m),function(x) {substring(x,1,tmp[length])})}
        }
    else {
        if(alpha==TRUE){
            changed_name <- sapply(colnames(mut_sum.m),function(x) {tolower(substring(x,1,tmp[length]+1))})}
        else{changed_name <- sapply(colnames(mut_sum.m),function(x) {tolower(substring(x,1,tmp[length]))})}
        }
    if (pattern == "-"){
        changed_name <- sapply(changed_name,function(x) gsub(".","-",x,fixed=TRUE))
    }
    return(changed_name)
}

RNA_norm_trans<-function(RNA_norm.m,cap=TRUE,pattern=".",length=4,alpha=TRUE){ #TCGA.05.4382.01A.01R.1206.07
    tmp = c(4,7,12,15)
    if (cap==TRUE){
        if(alpha==TRUE){
            changed_name <- sapply(colnames(RNA_norm.m),function(x) {substring(x,1,tmp[length]+1)})}
        else{changed_name <- sapply(colnames(RNA_norm.m),function(x) {substring(x,1,tmp[length])})}
        }
    else {
        if(alpha==TRUE){
            changed_name <- sapply(colnames(RNA_norm.m),function(x) {tolower(substring(x,1,tmp[length]+1))})}
        else{changed_name <- sapply(colnames(RNA_norm.m),function(x) {tolower(substring(x,1,tmp[length]))})}
        }
    if (pattern == "-"){
        changed_name <- sapply(changed_name,function(x) gsub(".","-",x,fixed=TRUE))
    }
    return(changed_name)
} 
