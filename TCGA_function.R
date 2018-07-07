##author: Rayna
##time: 20180705
##last update: 20180705
##descript: TCGA内需要用到的R function


##BothGeneMatrixGeneration 生成gene1 vs gene2矩阵
BGMatrixG<-function(gene_list,colname){
	BGMatrix <- matrix(nrow=length(gene_list)^2,ncol=length(colname))
	colnames(BGMatrix) <- colname
	BGMatrix[,1] <- rep(gene_list,each=length(gene_list))
	BGMatrix[,2] <- rep(gene_list,length(gene_list))
	return(BGMatrix)
}

##对RNAnorm文件第一列的Symbol|Entrez的分割
RNAnormTrans<-function(colname_list,type){
	if (type=="Entrez"){
		sapply(RNA_norm[2:nrow(RNA_norm),1],function(x) {strsplit(as.character(x),split="|",fixed=TRUE)[[1]][2]})
	}
	else if (type=="Symbol"){
		sapply(RNA_norm[2:nrow(RNA_norm),1],function(x) {strsplit(as.character(x),split="|",fixed=TRUE)[[1]][1]})
	}
}



##GeneGene(both,either,neither,only.... sample selection) 双基因样品选择
GeneGeneSample<-function(gene1,gene2,mut_sum.m,both_gene_out.m){
	result.lv <- list()
	gene1_sample.idx <- which(mut_sum.m[which(rownames(mut_sum.m)==gene1),]==1)
	both_sample.idx <- intersect(which(mut_sum.m[which(rownames(mut_sum.m)==gene1),]==1),which(mut_sum.m[which(rownames(mut_sum.m)==gene2),]==1))
    both_sample.v <- colnames(mut_sum.m)[both_sample.idx]
    gene2_sample.idx <- which(mut_sum.m[which(rownames(mut_sum.m)==gene2),]==1)
    either_sample.idx <- union(which(mut_sum.m[which(rownames(mut_sum.m)==gene1),]==1),which(mut_sum.m[which(rownames(mut_sum.m)==gene2),]==1))
    gene1_only_sample.idx <- setdiff(gene1_sample.idx,both_sample.idx)
    gene1_only_sample.v <- colnames(mut_sum.m)[gene1_only_sample.idx]
    gene2_only_sample.idx <- setdiff(gene2_sample.idx,both_sample.idx)
    gene2_only_sample.v <- colnames(mut_sum.m)[gene2_only_sample.idx]
    neither_sample.idx <- intersect(which(mut_sum.m[which(rownames(mut_sum.m)==gene1),]==0),which(mut_sum.m[which(rownames(mut_sum.m)==gene2),]==0))
    neither_sample.v <- colnames(mut_sum.m)[neither_sample.idx]
    row.idx <- intersect(which(both_gene_out.m[,1]==gene1),which(both_gene_out.m[,2]==gene2))
    result.lv <- list(gene1.idx=gene1_only_sample.idx,gene2.idx=gene2_only_sample.idx,both.idx=both_sample.idx,neither.idx=neither_sample.idx,
    	gene1.sample=gene1_only_sample.v,gene2.sample=gene2_only_sample.v,both.sample=both_sample.v,neither.sample=neither_sample.v,row.idx=row.idx)
}

countNA <- function(tmp.v){
    out <- length(which(is.na(tmp.v)));
    return(out);
}

