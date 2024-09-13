#1，构建索引index，从ucsc浏览器或者ENSEMBL下载参考基因组fa文件、基因注释gtf文件；假设为GRCh38.fa、GRCh38.gtf
#建议使用ensembl，最好是grch38也就是hg38

#2，获取cDNA序列，gffread从GTF/GFF3文件中提取序列
gffread GRCh38.gtf -g GRCh38.fa -w GRCh38.transcript.fa.tmp

# gffread生成的fasta文件同时包含基因名字和转录本名字
grep '>' GRCh38.transcript.fa.tmp | head

#3，去掉空格后面的字符串，保证cDNA文件中fasta序列的名字简洁，不然后续会出错
cut -f 1 -d ' ' GRCh38.transcript.fa.tmp >GRCh38.transcript.fa

#4，获取所有基因组序列的名字存储于decoy中
grep '^>' GRCh38.fa | cut -d ' ' -f 1 | sed 's/^>//g' >GRCh38.decoys.txt

#5，合并cDNA和基因组序列一起，注意：cDNA在前，基因组在后
cat GRCh38.transcript.fa GRCh38.fa >GRCh38_trans_genome.fa

#6，salmon index构建索引 （更慢，结果会更准）
salmon index -t GRCh38_trans_genome.fa -d GRCh38.decoys.txt -i GRCh38.salmon_sa_index

#salmon index：Salmon软件的索引构建命令。
#-t GRCh38_trans_genome.fa：指定合并后的cDNA和基因组序列文件。
#-d GRCh38.decoys.txt：指定decoy序列文件，用于提高比对的准确性。
#-i GRCh38.salmon_sa_index：指定输出的索引文件名。

#以上到index步骤所有部分的代码实际上就是参考https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/，对照中文做调整https://mp.weixin.qq.com/s/aDxnkK3jhyP6E_wfvza6jQ?search_click_id=11630710159653975537-1726118237844-0581657105

#7，对上游RNA-seq的fq文件进行质控QC清洗，fastqc/multiqc之后再trim_galore，如果html中报告error不能再改变就不用再QC了，获取clean fq
#建议参考ChIP、ATAC-seq等分析模块上游的QC脚本

#8，salmon quant定量化,下面的参数都是针对index+quant的mode而言的，如果使用STAR+salmon，该脚本就没有参考价值了
#对于双端测序PE，
fq1=$1
fq2=$2
output_folder=$3
salmon quant --gcBias -l A -1 ${fq1} -2 ${fq2}  -i GRCh38.salmon_sa_index路径 -g GRCh38.gtf路径 -o ${output_folder} -p 15
#-p是线程
#--gcBias，矫正GC偏差用，一般与上游multiqc报告结合使用，但是you can simply run Salmon with --gcBias in any case，所以都加上
#-l A允许Salmon自动推断测序文库类型
#-i index索引路径
#-g 提供转录本到gene的映射文件，也可以提供上面的gtf文件，只要提取其中的transcript_id以及gene_id两列即可

#如果是单端测序SE，
salmon quant --gcBias -l A -r ${fq1}   -i GRCh38.salmon_sa_index路径 -g GRCh38.gtf路径 -o ${output_folder} -p 20
