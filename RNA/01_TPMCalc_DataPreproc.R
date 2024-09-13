library("readr") #数据处理
library(tximport)  #从各种RNA测序计数器中导入数据的R包，比如说salmon定量结果
library(readr) #读取文件


#以下3行代码按照自己的数据要求来进行参考基因组的加载
#library(AnnotationHub)  
library(ensembldb)  #用于与Ensembl数据库交云的R包
library(EnsDb.Hsapiens.v75)  #特定版本的Ensembl人类基因组数据库
#目前暂时使用ensemble的注释，hg38版本按照源码是直接ah之后选择第一个注释，但是hg19需要另外手动加载注释数据库v75，不能直接检索

#AnnotationHub 是一个更广泛的资源发现和访问工具，它不仅限于 Ensembl 数据。
#ensembldb 是一个专门用于处理 Ensembl 数据的包，它允许用户创建和管理本地数据库副本。
#EnsDb.Hsapiens.v75 是一个具体的数据库对象，它是通过 ensembldb 包创建的，并且包含了特定版本的 Ensembl 人类基因组注释数据


#都是为了安装RNASeqBits，提供了一系列工具和函数，用于帮助RNA测序分析流程中的各种任务
library("devtools")
#install_github("AmyOlex/RNASeqBits")
library(RNASeqBits) # https://github.com/AmyOlex/RNASeqBits

library(NMF) #NMF（非负矩阵分解）包提供函数来进行非负矩阵分解，这是一种数据降维技术，常用于基因表达数据分析中
library(limma) #limma 是一个用于基因表达微阵列和RNA测序数据分析的R包，它提供了一系列的统计方法，用于识别差异表达的基因和进行基因表达分析，后面在脚本2中改用了edgeR


# quan.sf文件放置的地方,目前全盘使用ensemble里的fa+gtf文件
setwd("/mnt/disk4/haitao/bysj_seu/geo_data/hic/script8/RNA/ens_fa_gtf/")

## 获取指定目录下的所有文件的路径,其实就是前面的setwd
files <-  list.files(".", pattern = "\\.sf$", full.names = TRUE)
#files <- file.path(".", list.files("."))
## 去除前后缀
names(files) <- sub("_quant.sf", "", sub("./", "", files))
#到时候可以根据需求只深入BT549 vs HMEC部分
#目前策略是读入数据的时候只读入HEMC vs BT549的数据，或者是全部读入数据之后后期提取出HMEC vs BT549的数据！！！

## Build a tx2gene object that will summarize all transcript expressions to the gene level
## 构建一个tx2gene对象，该对象将把所有转录本的表达量总结到基因水平,这里相比原代码就是修改了ensemble数据库版本，整体上这些数据都能接下去用
#这里tx2gene一共有两列，第一列是转录本名字tx_id——ENST，第二列是gene名字gene_id——ENSG，
#然后我们输入的quant.sf文件的第一列name也是转录本（但是是由fa文件决定的）——ENST
EnsDb_h <-EnsDb.Hsapiens.v75
df_human <- transcripts(EnsDb_h,return.type = "DataFrame") #从 EnsDb_h 中提取转录本信息，并将结果存储在 df_human 中，这是一个数据框形式的转录本注释信息
tx2gene_h <- df_human[,c("tx_id","gene_id")] #从转录本数据框中选取 tx_id（转录本ID）和 gene_id（基因ID）两列，用于构建转录本到基因的映射关系

## Import Salmon files and summarize to the gene level
## Note it will only recognize the transcript IDs that are in the txi list.  Thus, importing the same files that contain both
txi_h <- tximport(files, type = "salmon", tx2gene = tx2gene_h, ignoreTxVersion = TRUE, importer = function(x) read_tsv(x, col_types="cnnnn"))
#使用 tximport() 函数将转录本表达量汇总到基因水平，并创建一个 txi_h 对象来存储汇总后的表达量信息，tx2gene = tx2gene_h 将之前构建的转录本到基因的映射关系应用到表达量汇总过程中
#我们手头拥有的是转录本数据，然后通过这个映射将转录本的表达量汇总到了gene水平

#### Re-calculate the TPM values from the given counts and lengths
#### 重新计算给定计数和长度值的TPM（Transcripts Per Million）值
source("/mnt/disk4/haitao/bysj_seu/geo_data/hic/script8/RNA/RNA-seq_script/calc.tpm.fromSalmon.R") #一个R包，从github上下载到本地
#名为 calc.tpm.fromSalmon 的函数，用于从Salmon和tximport产生的RNA-seq数据的汇总数据中计算TPM（每百万转录本数）值
#这个函数的作用是将Salmon和tximport产生的RNA-seq数据的汇总数据转换为TPM表达值，这在基因表达分析中是常用的一种标准化方式
TPM_h <- calc.tpm.fromSalmon(txi_h, colnames(txi_h$counts)) #总之是获得了TPM值

#用round()函数将txi_h$counts转换为一个数据框，并将结果保存在变量counts_h中。这里使用round()函数是为了将计数值取整
##保留原始数据count,以及手头标准化之后的TPM数据
counts_h <- round(as.data.frame(txi_h$counts))
#names(counts_h) <- names(TPM_h)

################### Save raw TPM values for GEO ####################
#合并基因信息和RNA-Seq数据，将基因符号添加到TPM+count数据——因为原始的TPM以及count只有ENSG这种形式的gene id，需要整理通俗的gene id
## 从注释数据库中提取基因符号，并将其添加到TPM数据中，第一列
df_human_gene <- genes(EnsDb_h, return.type = "DataFrame")[,c("gene_id","symbol")]
#将基因信息与之前计算的TPM值和原始计数数据合，是TPM
TPM_h2_annot <- merge(df_human_gene, TPM_h, by.x="gene_id", by.y="row.names", all.x=FALSE, all.y=TRUE)
#将基因符号与原始计数数据进行合并，是counts
counts_h_annot <- merge(df_human_gene, counts_h, by.x="gene_id", by.y="row.names", all.x=FALSE, all.y=TRUE)

## 保存文件，作为原始文件，所以也就只是保留了带有gene注释的TPM+count数据
write.table(TPM_h2_annot, file="TPM_allgenes_HMECvsTNBC.tsv", quote=FALSE, sep="\t", row.names=FALSE)
#对应processed_data/TPM_human_AllGenes_2021_BulkRNASeq_10.28.21.tsv
write.table(counts_h_annot, file="Count_allgenes_HMECvsTNBC.tsv", quote=FALSE, sep="\t", row.names=FALSE)
#对应processed_data/Count_human_AllGenes_2021_BulkRNASeq_10.28.21.tsv


#下面的分析实际上并没有在后续代码分析中使用过，但是还是照样处理
##################### Process TPM values for Analyses ###################
#对TPM值进行进一步的处理，以便进行后续的分析
## 移除所有样本中表达为零的基因，确保在后续分析中只使用具有表达的基因，从count入手处理TPM，只有TPM是用于后续分析的
counts_h_noZero <- filter.removeZeroSumRows(counts_h)
TPM_h_noZero <- TPM_h[row.names(counts_h_noZero),]

## 从注释数据库中提取基因符号，并将其添加到原始数据中，第一列，但是前面已经做过了
#df_human_gene <- genes(EnsDb_h, return.type = "DataFrame")[,c("gene_id","symbol")]
#同理合并，只不过是特殊的去除0的count以及TPM数据，操作都是和前面同样的
counts_h_noZero_annot <- merge(df_human_gene, counts_h_noZero, by.x="gene_id", by.y="row.names", all.x=FALSE, all.y=TRUE)
TPM_h_noZero_annot <- merge(df_human_gene, TPM_h_noZero, by.x="gene_id", by.y="row.names", all.x=FALSE, all.y=TRUE)

#对处理后的TPM数据取对数（以2为底），并将基因符号（Gene Symbols）添加到数据中，以便更好地理解和分析数据
#还是同理，对上面的去0之后的TPM数据进行对数化，然后还是TPM合并，操作是同样的
TPMlog2_h_noZero <- log2(TPM_h_noZero+1)
TPMlog2_h_noZero_annot <- merge(df_human_gene, TPMlog2_h_noZero, by.x="gene_id", by.y="row.names", all.x=FALSE, all.y=TRUE)

##保存文件,首先是去除0的已表达的count和TPM
write.table(counts_h_noZero_annot, file="Count_rm0_HMECvsTNBC.tsv", quote=FALSE, sep="\t", row.names=FALSE)
#对应processed_data/Count_human_ZeroExpRemoved_2021_BulkRNASeq_10.28.21.tsv
write.table(TPM_h_noZero_annot, file="TPM__rm0_HMECvsTNBC.tsv", quote=FALSE, sep="\t", row.names=FALSE)
#对应processed_data/TPM_human_ZeroExpRemoved_2021_BulkRNASeq_10.28.21.tsv

##在前面基础上的log2转化的TPM
write.table(TPMlog2_h_noZero_annot, file="TPM_log2_rm0_HMECvsTNBC.tsv", quote=FALSE, sep="\t", row.names=FALSE)
#processed_data/log2TPM_human_ZeroExpRemoved_2021_BulkRNASeq_10.28.21.tsv

#继续对上面的TPM数据进行操作
## 对已经进行对数转换的TPM数据进行上四分位数归一化（Upper Quantile Normalization），以便更好地进行比较和分析
## 对只表达的gene的数据，计算上四分位数+归一化数据，总之是差异表达分析所需要的归一化分析
TPMlog2_h_noZero.quantileAll <- apply(TPMlog2_h_noZero, 2, function(x){quantile(x[x>0], 0.75)})
TPMlog2_h_noZero.norm <- as.data.frame(t(t(TPMlog2_h_noZero) / TPMlog2_h_noZero.quantileAll))
TPMlog2_h_noZero.norm_annot <- merge(df_human_gene, TPMlog2_h_noZero.norm, by.x="gene_id", by.y="row.names", all.x=FALSE, all.y=TRUE)

write.table(TPMlog2_h_noZero.norm_annot, file="TPM_QuantileNorm_log2_rm0_HMECvsTNBC.tsv", quote=FALSE, sep="\t", row.names=FALSE)
#对应processed_data/log2TPM_human_ZeroExpRemoved_QuantileNorm_2021_BulkRNASeq_10.28.21.tsv




