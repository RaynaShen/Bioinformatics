##author: yiru
##time: 20180704
##last update: 20180704
##descript；将geo上的series matrix分成pheno文件和exp文件两类。

import os

def get_fplist(dir_list):
    fplist = []
    for folder_dir in dir_list:
        for x in os.listdir(folder_dir):
            fp = folder_dir + x
            fplist.append(fp)
    return(fplist)

def SepPhenoExp(ContentList,OutputFoldPheno,OutputFoldExp,OriFilePath): ##seperate the series matrix data into PhenoData and ExpData
    for line in ContentList:
        if line.startswith("!Series_platform_id"):
            platform = line.split("\t")[1][1:-1]
    head_idx = ContentList.index("")
    head_end_idx = ContentList.index("!series_matrix_table_begin")
    pheno = ContentList[head_idx+1:head_end_idx]
    exp = ContentList[head_end_idx+1:-1]
    OutPath = OutputFoldPheno + str(OriFilePath.split("/")[-1].split("_")[0]) + "_pheno.txt"
    if os.path.exists(OutPath):
        pass
    else:
        fw = open(OutputFoldPheno + str(OriFilePath.split("/")[-1].split("_")[0]) + "_pheno.txt","wb")
        fw.write(platform + "\n")
        for i in pheno:
            fw.write(i +'\n')
        fw.close()
    OutPath = OutputFoldExp + str(OriFilePath.split("/")[-1].split("_")[0]) + "_exp.txt"
    if os.path.exists(OutPath):
        pass
    else:
        fw = open(OutputFoldExp + str(OriFilePath.split("/")[-1].split("_")[0]) + "_exp.txt","wb")
        for i in exp:
            fw.write(i +'\n')
        fw.close()

def PhenoExtract(PhenoFileDir): ##extract the pheno charastics type list for each data set. 
    OutPath = PhenoFileDir[:-10] + "_PhenoList.txt"
    if os.path.exists(OutPath):
        pass
    else:
        with open(PhenoFileDir) as f:
            pheno = [line.rstrip("\n") for line in f]
        PhenoList = []
        for line in pheno[1:]:
            title = line.split("\t")[0]
            if title == "!Sample_characteristics_ch1":
                PhenoList.append(line.split("\t")[1].split(":")[0][1:])
        with open(PhenoFileDir[:-10] + "_PhenoList.txt","wb") as w:
            w.write(pheno[0] + "\n")
            w.write("\t".join(PhenoList) + "\n")
        print(pheno[0])
        print("find pheno type #: %s" %(len(PhenoList)))


def main(OriFilePath):
    with open(OriFilePath) as f:
        tmp = [line.rstrip("\n") for line in f] 
    SepPhenoExp(ContentList=tmp,OutputFoldPheno="/work/home/OM/shenyr/gene_gene/dataset/PhenoFile/",OutputFoldExp="/work/home/OM/shenyr/gene_gene/dataset/ExpFile/",OriFilePath=item)
    PhenoFileDir = "/work/home/OM/shenyr/gene_gene/dataset/PhenoFile/" + str(item.split("/")[-1].split("_")[0]) + "_pheno.txt"
    PhenoExtract(PhenoFileDir = PhenoFileDir)


fplist = get_fplist(["/work/home/OM/shenyr/gene_gene/dataset/"])
outputdir = "/work/home/OM/shenyr/gene_gene/dataset/PhenoFile/"
for item in fplist:
    if item.endswith("_matrix.txt"):
       main(item)




