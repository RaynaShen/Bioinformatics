##author: Rayna
##time: 20180705
##last update: 20180705


MethGeneExtr <- function(meth.m,gene_list.v=FALSE,cutoff.v){ 
##从450K中提取位于选定基因promoter区域的探针向量。meth.m为行名探针的甲基化矩阵。gene_list.v为gene_symbol,为FALSE时表示取450K中全基因;cutoff.v为基因位置,如：c("5'UTR","TSS200","TSS1500")
    load("/work/home/OM/shenyr/gene_gene/exp_mut/probe_gene.Rd") ###450K anno_info: probe_gene.lv
    target.idx <- c()
    if (gene_list.v == FALSE){
        gene_list.v <- vector()
        for (i in 1:length(probe_gene.lv)){
            gene_list.v <- c(gene_list.v,probe_gene.lv[[i]]$gene)
        }
        gene_list.v <- gene_list.v[!duplicated(gene_list.v)]
    }
    for (i in 1:length(probe_gene.lv)){
        if (length(intersect(probe_gene.lv[[i]]$gene,gene_list.v))>=1){
            if (length(probe_gene.lv[[i]]$gene==1)){
                if (probe_gene.lv[[i]]$position %in% cutoff.v){
                    target.idx <- c(i,target.idx)
                    #print(paste(probe_gene.lv[[i]]$gene,probe_gene.lv[[i]]$position,sep=": "))
                }
            }
            else{for (j in 1:length(probe_gene.lv[[i]]$gene)){
                if (probe_gene.lv[[i]]$gene[j] %in% gene_list.v){
                    if (probe_gene.lv[[i]]$position[j] %in% cutoff.v){
                        target.idx <- c(i,target.idx)
                        #print(paste(probe_gene.lv[[i]]$gene[j],probe_gene.lv[[i]]$position[j],sep=": "))
                        }
                    }
                }
            }
        }
        print(paste0("done with probe: ",i))
    }
    probe_sel.v <- vector()
    for (i in 1:length(target.idx)){
        probe_sel.v <- c(probe_sel.v,probe_gene.lv[[i]]$probe)
    }
    tmp.idx <- match(probe_sel.v,rownames(meth.m))
    non_na.idx <- which(is.na(tmp.idx)==FALSE)
    meth_sel.m <- meth.m[tmp.idx[non_na.idx],]
    result.lv <- list(meth_sel=meth_sel.m,probe_sel=probe_sel.v)
    return(result.lv)
}


