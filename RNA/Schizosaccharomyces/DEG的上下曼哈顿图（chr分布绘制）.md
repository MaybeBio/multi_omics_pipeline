记一次商务工作：

整体脚本文件记录参考：multi_omics_pipeline/RNA/Schizosaccharomyces（zip压缩）  
![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1732885446422-f1f81d4d-a156-4a1d-be09-2e1be2d5721b.png)

![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1732885471885-78eae7cd-6c5f-4804-8bc2-1c3b62ac45e6.png)

![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1732885493981-59bae777-637f-4706-b29e-30c9db282479.png)

![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1732885508339-68582af6-55ce-4d30-a001-c8adeb6f1a69.png)

![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1732885523492-14ad9928-0bb3-4eb1-a608-f0aae76a71b3.png)

![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1732885547247-cfcad8cc-41ac-405e-ba5c-a096d41d5d5b.png)

![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1732885597354-c72b8bd2-3d83-4a61-8ae5-4e99c23ed05b.png)

![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1732885637478-1cdb5979-9b13-4e7a-9f38-6bf58a6765a2.png)

总之就是一个酿酒酵母物种的一个DEG文件，然后要求按照绘制曼哈顿图那样，但是上调绘制在上侧（类似曼哈顿图），然后下调绘制在下侧；

需要参考基因组的注释文件（原始序列fa文件也可以，作为绘图的x轴参考标度），

以及DEG文件中的y轴（FC或者是logFC作为标度），

但是DEG文件中并没有基因的坐标数据，也就是回帖到染色体上的x轴；

所以需要从参考基因组的注释文件中自己提取。

该客人提供了一个公共网站，也就是包含该物种的参考基因组信息、注释信息所在的一个数据库：

[https://www.pombase.org/data/](https://www.pombase.org/data/)

![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1732948734405-25a0c41c-baa2-429b-9f6e-4dfcf4ec52a0.png)

因为该客人关注的DEG中大部分是非编码gene，所以CDS序列坐标是不能用了（因为CDS是蛋白质编码gene，所以坐标是收集不齐的）；

然后exon的坐标也不能用，因为exon的坐标实际上是离散的，也就是同一个gene id会有很多行，就是同一个gene会有很多个转录本；

其实在search中一搜gene id就能够知道这个gene到底有多少个转录本以及对应的详细的坐标信息：  
![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1732949039468-2f4a9067-6774-47f0-8e02-fdab4d781873.png)

比如说：  
![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1732949104499-94dd6be2-cc58-4350-bc57-679931184597.png)

![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1732949114323-70647094-5d3b-4ee0-b781-5dd4e551641b.png)

所以基本上只能用gff3中的注释信息（基本上参考基因组以及所有gene id对应的坐标信息都得从这里来获取）



1，首先是清洗DEG数据，毕竟这位客人对其下游数据的清理可谓是惨不忍睹：  
（1）首先是要读入到tibble中的各列数据格式清洗：  
str等都可以看，其中几乎所有数值列都是chr（我tmd），可以用as.numeric，当然直接excel中偷懒也差不多

![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1732949384493-7f85d43e-9e4a-4317-8156-1492e357dc6a.png)

![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1732949404797-1530ff05-7e40-4bff-82df-879f6e70c0c9.png)

![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1732949420738-d7443769-88ec-44b8-9a14-da800f688852.png)

其余的简单的数据操转换等，直接tidverse，可以参考我之前的博客

[https://blog.csdn.net/weixin_62528784/article/details/142744220](https://blog.csdn.net/weixin_62528784/article/details/142744220)

（2）对参考基因组数据的获取以及清洗：

gff3数据格式一览：  
![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1732949638148-1be62884-6bc6-4911-80eb-f43490e5d066.png)

我们主要关注的就是第3列是gene的那一行，其余接下来的几行其实都是对这个gene的描述，然后坐标其实就在第4/5列，第一列是chr，最后一列是对该gene（object）属性的描述，我们要的ID在第一个分号前面

提取坐标信息其实很简单，直接上awk，下面提供其中一种处理方式

```plain
#从gff3文件最后一列获取gene id
awk '$3 == "gene" {print $NF}' Schizosaccharomyces_pombe_all_chromosomes.gff3 | awk -F ";" '{print $1}' | awk -F "=" '{print $2}' > gene_id.bed
#找到第3列是gene的行，只打印最后1列，然后这一列我们用分号读入，只取第一列（也就是第一个分号前面的字段），然后再用等于号读入，只取第2个字段，那么获取的就是gene_id了
#指定按照什么分隔符读入是-F ""选项，指定按照什么分隔符输出是BEGIN {OFS="\t"}，此处只有1列所有不指定输出格式，输入倒是注意

#获取bed坐标
awk 'BEGIN {OFS="\t"}$3 == "gene" {print $1, $4, $5}' Schizosaccharomyces_pombe_all_chromosomes.gff3 > gene_coord.bed
#这里比较简单，其实就是按照制表符输出格式打印chr、start、end部分

#tsv拼接paste不指定分隔符的话默认是tsv拼接，当然也可以使用awk进行复杂的文件拼接等
paste gene_coord.bed gene_id.bed > ref_gene.tsv

#当然也可以使用awk进行一步处理等，此处略
```

（3）使用left_join取交集数据

2，正式绘图：  
该酵母物种只有3条染色体chr，分开绘制

（1）x轴标度：用gene bed中的start、end都可以，当然我也算了peak_summit（gene body center）

反正gene那么小，然后绘制起来整个chr长度缩放起来其实看起来都一样

（2）y轴标度使用logFC，因为反正要画上下曼哈顿，用FC幅度上不统一，很容易出现极端值，而且自己加个sign也是费事

（3）参考基因组的长度：  
还是使用samtools工具，其实就是获取chr_sizes文件

```plain
#主要是获取每条chr的长度
samtools faidx Schizosaccharomyces_pombe_all_chromosomes.fa 
#获取fai文件

awk '{print $1 "\t" $2}' Schizosaccharomyces_pombe_all_chromosomes.fa.fai > chr_size.bed
```

![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1732950593468-9f6bb08a-2b84-4ec9-b6d5-8bd5c3f29098.png)

其实绘制出来之后可以发现：chr1/2/3总长（**这里应该是总长，没必要将其中特殊区域提取出来再算剩余长度**）

还有一些特殊区域的长度，后面其实也没用上

（4）客人要求的特殊区域：

主要是端粒、着丝粒，以及mating_type region（上面提到的）  
然后具体坐标后面客人才提供（我tmd，我还在注释文件中吭哧吭哧想办法计算整合）

虽然数据上有出入，但那不是我的职责，所以按照提供的坐标去绘制这些区域就是了

![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1732950904991-63dd8014-da78-45ff-a446-3750ec5d17b4.png)

然后就是前后延伸5kb（其实就是侧翼序列flank），当然其实加5kb没有区别

3，草稿效果：  
![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1732885446422-f1f81d4d-a156-4a1d-be09-2e1be2d5721b.png?x-oss-process=image%2Fformat%2Cwebp)

![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1732951025531-c3202a7b-a897-4c5c-8cee-ad3edc68a1eb.png)

（1）其实我一开始就觉得都能分上下了就没必要给上下调的配色了，不过无所谓

后面出图的时候：有一种就是我没给上下调配色，只给gene类型分组配色（其实就是coding、non-coding分组）





![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1732951212554-859452e4-60c4-4b14-a775-dc501b491684.png)

（2）还有个显著性q-value，其实我一直都想用上，该客人压根门外汉不懂

所以我就直接硬上了显著性气泡图，

指标可以用倒数1/q，或者是-log10





（3）分组颜色：

就像1中说的那样，up以及down既然绘制了上下相对位置，就没有必要分组；

然后就是蛋白质编码gene与非编码gene之间的分组，这个我觉得还合理点，所以其实可以上下分轴的时候绘制一个gene类型分类；

当然如果硬要分成2x2分组分类的话，也可以在ggplot读入的时候使用unite合并变量

![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1733031442911-0f098957-939c-4da7-bde7-e26b70929449.png)

![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1733031453558-a246d588-9d81-4233-b4bc-8b2689bfd6cd.png)

（4）另外gene名字可以ggrepl添加一下：其实效果不好，可以不加

可以使用geom_text_repel或geom_label_repel

![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1733030869469-1567941b-bbe6-4f6d-ad74-f6f35d0a5f79.png)

（5）另外就是x轴上面的标度（scale_x_continuous(breaks)），

另外是否要用scale_x_log10，因为考虑到有gene分布比较零散，所以不建议使用scale_x_log10，

毕竟着丝粒以及端粒等特殊区域实际上更小，除非是跳跃展示区域，否则缩放没有意义

然后chr长度可以在矩形上限制，见（6）

![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1733031464825-72d0b0dd-16ad-435c-85f7-86e315d8e7d0.png)

（6）如何将chr等区域绘制到y=0上：

使用geom_rect+annotate在y=0直线上绘制矩形区域，以表征chr，

展示效果类似于核型图，其实有很多绘图包可以使用现成的核型绘图  
![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1733031361255-cae11c4b-5ca7-46f3-b335-c15af8477bb5.png)

注意该区域标度可以限制

（6）后期AI或者福熙修图：  
主要是x轴去除，将标度移动到y=0处，position坐标可以进一步使用（5）中的进行细化

![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1733031122667-3ab2020e-1c08-4dcc-bb23-832f2e851fc5.png)

其中chrx是对上小组配色，简单绘制绘图

chrx_facet是使用简单分面

chrx_41是不分组上下调，但是对基因类型coding与否进行分组配色

chrx_42是使用4分组

chrx_42_2是使用不同的配色组

chrx_42_label/text是在chrx_42基础上配色

基本草稿效果类似：  
![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1733031535269-da027631-7a1a-4ed8-966c-0b8888229125.png)

  












