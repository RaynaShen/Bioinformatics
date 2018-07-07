##author:yiru
##time:20180703
##last update:20180704
##descript：对geo上的表达量文件进行基于基因的平均表达量处理
import json

def PlatformAnno(PlatName,GeneSymbol):
    """
    生成平台json annotation文件,platform_name类似：GPL570-55999, colname_line是在annotation文件中character名称的那一行
    GPL570为16，GPL0379为9;GeneSymbol为在annotation文件中，genesymbol是如何表示的，如：“GeneSymbol”，“Gene Symbol”
    """
    with open('/work/home/OM/shenyr/gene_gene/dataset/platform/'+PlatName+'.txt') as f:
        plat_anno = [line.rstrip("\r\n") for line in f]
    probe_list = []
    anno_info = []
    idx = -1
    for line in plat_anno:
        if line.startswith("ID"):
            start_index = plat_anno.index(line)
    with open('/work/home/OM/shenyr/source/Entrez_GeneSybo.txt') as f:
        Entrez_Symbo = [line.rstrip("\r\n") for line in f]
    Entrez = [x.split(" ")[2][1:-1] for x in Entrez_Symbo[1:]]
    Symbo = [x.split(" ")[1][1:-1] for x in Entrez_Symbo[1:]]
    Entrez_Symbo_dict = dict(zip(Entrez,Symbo))
    for line in plat_anno:
        if line.split("\t")[0] == "ID":
            colnames = line.split("\t")[1:]
        idx += 1
        if PlatName == "GPL22085":
            probe_list.append(line.split("\t")[0])
            tmp_line = line.split("\t")[1:]
            Symbo = Entrez_Symbo_dict[tmp_line[0]]
            tmp_line.append(Symbo)
            anno_info.append(tmp_line)
        else:
            if idx > start_index:
                probe_list.append(line.split("\t")[0])
                anno_info.append(line.split("\t")[1:])
    anno_info_dict = dict(zip(probe_list,anno_info))
    for key in anno_info_dict.keys():
        if PlatName == "GPL22085":
            colnames.append("Gene_Symbol")
        anno_info_dict[key] = dict(zip(colnames,anno_info_dict[key]))
    with open('/work/home/OM/shenyr/gene_gene/dataset/platform/'+PlatName+'.json',"wb") as outfile:
        json.dump(anno_info_dict,outfile)
        outfile.write('\n')
    probe_gene_dict = {}
    for probe in anno_info_dict.keys():
        the_gene = anno_info_dict[probe][GeneSymbol]
        gene_list = the_gene.split("///")
        for gene in gene_list:
            if gene in probe_gene_dict:
                probe_gene_dict[gene].append(probe)
            else:
                probe_gene_dict[gene] = [probe]
    with open('/work/home/OM/shenyr/gene_gene/dataset/platform/'+PlatName+'_PG.json',"wb") as outfile:
        json.dump(probe_gene_dict,outfile)
        outfile.write('\n')

PlatformAnno(PlatName="GPL570-55999",GeneSymbol="Gene Symbol")
PlatformAnno(PlatName="GPL10379",GeneSymbol="GeneSymbol")
PlatformAnno(PlatName="GPL10558-50081",GeneSymbol="ILMN_Gene")
PlatformAnno(PlatName="GPL14550-9757",GeneSymbol="GENE_SYMBOL")
PlatformAnno(PlatName="GPL22085",GeneSymbol="Gene_Symbol")
PlatformAnno(PlatName="GPL6480-9577",GeneSymbol="GENE_SYMBOL")
PlatformAnno(PlatName="GPL6884-11607",GeneSymbol="ILMN_Gene")
PlatformAnno(PlatName="GPL6947-13512",GeneSymbol="ILMN_Gene")
PlatformAnno(PlatName="GPL7280-24714",GeneSymbol="Gene Symbol")
PlatformAnno(PlatName="GPL96-57554",GeneSymbol="Gene Symbol")

def PlatAnnoProbExtr(PlatName,ExpDict,SampList,OutputPath,GenePick=None,):
    '''
    通过annotation文件，得到需要的gene在expression文件中的平均表达量文件
    '''
    with open("/work/home/OM/shenyr/gene_gene/dataset/platform/" + PlatName + "_PG.json",'r') as load_f:
        probe_gene_dict = json.load(load_f)
    gene_sel = probe_gene_dict.keys()
    if GenePick != None:
        gene_sel = GenePick   
    gene_sel_ave_exp_dict = {}
    for gene in gene_sel:
        probes = probe_gene_dict[gene]
        ave_exp = [0] *len(exp_dict[probes[0]])
        for probe in probes:
            sum_exp = [ave_exp[i] + exp_dict[probe][i] for i in range(len(ave_exp))]
        ave_exp = [x/len(probes) for x in sum_exp]
        gene_sel_ave_exp_dict[gene] = ave_exp
    with open(OutputPath,"wb") as w:
        w.write("\t".join(SampList) + "\n")
        for key in gene_sel_ave_exp_dict.keys():
            w.write(key + "\t" + "\t".join([str(x) for x in gene_sel_ave_exp_dict[key]]) + "\n")



########分别处理表达文件和pheno文件#################
########处理目标：得到基于基因的平均表达量矩阵
ID = "GSE8894"
exp_file = "/work/home/OM/shenyr/gene_gene/dataset/ExpFile/" + ID + "_exp.txt"
pheno_file = "/work/home/OM/shenyr/gene_gene/dataset/PhenoFile/" + ID + "_pheno.txt"
pheno_list_file = "/work/home/OM/shenyr/gene_gene/dataset/PhenoFile/" + ID + "_PhenoList.txt"

with open(exp_file) as f:
    exp = [line.rstrip("\n") for line in f]

sample_list = [x[1:-1] for x in exp[0].split("\t")[1:]]
probe_list = [x.split("\t")[0][1:-1] for x in exp[1:]]
exp_list = []
for line in exp[1:]:
    exp_value = [float(x) for x in line.split("\t")[1:]]
    exp_list.append(exp_value)

exp_dict = dict(zip(probe_list,exp_list))

with open(pheno_file) as f:
    pheno = [line.rstrip("\n") for line in f]

with open(pheno_list_file) as f:
    pheno_list = [line.rstrip("\n") for line in f]

####for test###############
#gene_sel = ["MAPK1","ADAM32","TP53","BRCA1"]
PlatAnnoProbExtr(PlatName="GPL570-55999",ExpDict=exp_dict,SampList=sample_list,OutputPath="/work/home/OM/shenyr/gene_gene/dataset/test/GSE8894_avexp.txt")



