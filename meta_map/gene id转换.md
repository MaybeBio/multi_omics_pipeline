涉及到gene的传统分析环节，<font style="color:rgb(0, 0, 0);">就是表达量矩阵数据处理，最频繁的就是传统的表达量芯片的差异分析和富集分析；</font>

<font style="color:rgb(0, 0, 0);">这些分析都是基于基因的，而基因有多种多样的id体系，而且不同的数据分析环节经常是需要进行id的转换！</font>

<font style="color:rgb(0, 0, 0);">#1，不同数据库的gene的id体系：  
</font><font style="color:rgb(0, 0, 0);">在基因组学和分子生物学研究中，基因的标识符是理解和交流基因信息的关键。关于基因的symbol、人类基因命名委员会（HGNC）、Entrez ID和Ensembl ID的解释：</font>![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1726103629917-0c2f8eeb-59cb-48ef-810a-969f47a73e23.png)

![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1726103642450-9fcdb00c-4044-447a-9243-18a3d77e907b.png)

理解这些标识符的重要性：

标准化：基因符号和ID提供了一种标准化的方式来引用基因，这对于科学出版物和数据库的一致性至关重要。

数据检索：这些标识符是检索基因相关信息的关键，如序列数据、表达模式、功能注释等。

跨数据库兼容性：不同的数据库可能使用不同的标识符系统，但它们通常可以相互转换，以便于数据的整合和分析。

在实际研究中，研究人员可能会根据需要使用这些不同的标识符来查找特定基因的信息，或者在比较不同研究结果时进行基因标识符的转换。



**<font style="color:rgb(0, 0, 0);">暂时没有统一的命名体系的基因</font>**

<font style="color:rgb(0, 0, 0);">基因命名和标识符系统是多样化的，反映了基因的不同类型和功能 ：</font>![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1726103762032-9af7d930-42b4-43f8-90b2-e7248325fcf7.png)![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1726103774757-0d31644e-faf7-4b07-ab83-d85003bb64df.png)









#2，如何转换gene id：

##（1）在线网页转换工具（比较权威的生物数据库）：  
<font style="color:rgba(0, 0, 0, 0.9);">①bioDBnet（Biological DataBase network，生物数据库网络）</font>

[https://biodbnet-abcc.ncifcrf.gov/db/db2db.php](https://biodbnet-abcc.ncifcrf.gov/db/db2db.php)

![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1726103916620-86a1d258-7ead-4abd-b00c-fe1384ddbfa4.png)

参考[https://mp.weixin.qq.com/s/VZ3gh-DlB0MqUfTpmRG80A](https://mp.weixin.qq.com/s/VZ3gh-DlB0MqUfTpmRG80A)

<font style="color:rgba(0, 0, 0, 0.9);">②Ensembl数据库是涵盖Human、Mouse、Zebrafish等模式动物，拟南芥、水稻等模式植物，真菌、原生动物、后生动物、细菌的基因组数据库。</font>

[https://asia.ensembl.org/index.html](https://asia.ensembl.org/index.html)

![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1726103952741-2cee3af9-3eef-4b0b-976a-e06ce7b59ebb.png)

参考[https://mp.weixin.qq.com/s/VZ3gh-DlB0MqUfTpmRG80A](https://mp.weixin.qq.com/s/VZ3gh-DlB0MqUfTpmRG80A)

<font style="color:rgb(63, 63, 63);">③DAVID（The Database for Annotation, Visualization and Integrated Discover）是一款整合了生物学数据和分析工具的生物信息数据库，为基因或蛋白提供系统的生物学注释信息。</font>

[https://david.ncifcrf.gov/home.jsp](https://david.ncifcrf.gov/home.jsp)

![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1726103990228-4a2bf643-590a-4685-89e9-1c1df6f81919.png)

参考[https://mp.weixin.qq.com/s/VZ3gh-DlB0MqUfTpmRG80A](https://mp.weixin.qq.com/s/VZ3gh-DlB0MqUfTpmRG80A)

<font style="color:rgb(63, 63, 63);">④g:Profiler中的g:convert可以实现各种基因，蛋白质，微阵列探针和许多其他类型的ID之间进行转换。涵盖60多个物种、至少40多种ID类型。</font>

[https://biit.cs.ut.ee/gprofiler/convert](https://biit.cs.ut.ee/gprofiler/convert)

![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1726104021458-cff11783-40c0-42c1-a7bb-442e587294d4.png)

参考[https://mp.weixin.qq.com/s/VZ3gh-DlB0MqUfTpmRG80A](https://mp.weixin.qq.com/s/VZ3gh-DlB0MqUfTpmRG80A)

<font style="color:rgb(63, 63, 63);">⑤KEGG（Kyoto Encyclopedia of Genes and Genomes）是一个整合了基因组、化学和系统功能信息的数据库。</font>

[https://www.kegg.jp/kegg/tool/conv_id.html](https://www.kegg.jp/kegg/tool/conv_id.html)

![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1726104054995-70c12fe3-5076-4318-a279-4107db8a38c0.png)

参考[https://mp.weixin.qq.com/s/VZ3gh-DlB0MqUfTpmRG80A](https://mp.weixin.qq.com/s/VZ3gh-DlB0MqUfTpmRG80A)

<font style="color:rgba(0, 0, 0, 0.9);">⑥Uniprot（Universal Protein）是信息最丰富、资源最广的蛋白质数据库。整合Swiss-Prot、 TrEMBL 和 PIR-PSD 三大数据库。</font>

[https://www.uniprot.org/id-mapping](https://www.uniprot.org/id-mapping)

![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1726104090127-500d344d-6b3a-4c93-8ec8-ca51477814be.png)

参考[https://mp.weixin.qq.com/s/VZ3gh-DlB0MqUfTpmRG80A](https://mp.weixin.qq.com/s/VZ3gh-DlB0MqUfTpmRG80A)

⑦metascape:  
[https://metascape.org/gp/index.html#/main/step1](https://metascape.org/gp/index.html#/main/step1)

##（2）使用R包进行gene id转换：

参考[https://mp.weixin.qq.com/s/2Z2qMnak7F2G_bo-9YbNxw](https://mp.weixin.qq.com/s/2Z2qMnak7F2G_bo-9YbNxw)

[https://mp.weixin.qq.com/s/Ba-KJc2_6PmEj-KtcQ9CGQ](https://mp.weixin.qq.com/s/Ba-KJc2_6PmEj-KtcQ9CGQ)

[https://blog.csdn.net/weixin_40739969/article/details/89354167](https://blog.csdn.net/weixin_40739969/article/details/89354167)

[https://www.biostars.org/p/22/](https://www.biostars.org/p/22/)  
①bioMart：除了上面的1个在线ensembl数据库，还有一个R包——biomaRt：  
参考[https://mp.weixin.qq.com/s/mb4rXd_sW0toPbGrJRs-Jg](https://mp.weixin.qq.com/s/mb4rXd_sW0toPbGrJRs-Jg)

②clusterProfiler::bitr：也就是bitr包：

参考[https://mp.weixin.qq.com/s/mb4rXd_sW0toPbGrJRs-Jg](https://mp.weixin.qq.com/s/mb4rXd_sW0toPbGrJRs-Jg)

③Genekitr：

参考[https://www.genekitr.fun/gene-id-conversion-1](https://www.genekitr.fun/gene-id-conversion-1)

④AnnotationDbi::mapIds：

参考[https://mp.weixin.qq.com/s/KsZGc4m8QXEZgKZbuSMhjg](https://mp.weixin.qq.com/s/KsZGc4m8QXEZgKZbuSMhjg)

⑤gprofiler：也有R包，和前面一样还有1个在线网站

##（3）自己写脚本：

参考[https://mp.weixin.qq.com/s/2Z2qMnak7F2G_bo-9YbNxw](https://mp.weixin.qq.com/s/2Z2qMnak7F2G_bo-9YbNxw)  
①use Web Crawler：

![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1726104771085-e70b2151-4269-43d2-b59d-0bde5af8029f.png)

参考[https://mp.weixin.qq.com/s/5u_7wJ1raSrltzum4RbSQA](https://mp.weixin.qq.com/s/5u_7wJ1raSrltzum4RbSQA)

②获取gtf/gff注释文件，然后自己写linux脚本（awk、sed、grep等3剑客）

参考[https://mp.weixin.qq.com/s/rKErPfTCFfcVJX_laGfpmg](https://mp.weixin.qq.com/s/rKErPfTCFfcVJX_laGfpmg)

![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1726104835909-4fcd3aed-e038-4728-a624-fbfd1f953788.png)

##（4）最常见的gene id转换：

ENSEMBL、SYMBOL、ENTREZID，这3种格式之间的转换

参考：

[https://mp.weixin.qq.com/s/mBNiTFKhGmLtMLWkrqLD8g](https://mp.weixin.qq.com/s/mBNiTFKhGmLtMLWkrqLD8g)

#3，id转换过程中会遇到的问题：

（1）转换数据不全：通常是转换之后id大幅损失（<font style="color:rgba(0, 0, 0, 0.9);">id匹配率低下问题）</font>

**建议优先用ensembl+bioMart**  
参考[https://mp.weixin.qq.com/s/kw2x2Tku53SiGeFAMoAwaw](https://mp.weixin.qq.com/s/kw2x2Tku53SiGeFAMoAwaw)

[http://www.bio-info-trainee.com/8032.html](http://www.bio-info-trainee.com/8032.html)

参考数据库、R包等都会影响  
<font style="color:rgba(0, 0, 0, 0.9);">org.Hs.eg.db一般会偏少，</font>

**应对策略：**

①使用<font style="color:rgba(0, 0, 0, 0.9);">biomart多点</font>

②使用<font style="color:rgba(0, 0, 0, 0.9);">genecode数据库的gtf注释，也就是从gtf上获取注释（但是注意哟啊使用匹配的gtf版本信息！！！）</font>

<font style="color:rgba(0, 0, 0, 0.9);">③TCGA的gene id转换问题可以去ucsc的xena浏览器里下载几乎完美匹配的id</font>

![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1726106066073-73a92767-702c-417d-af1c-a0b0a18d4d2d.png)

<font style="color:rgba(0, 0, 0, 0.9);">④如果用bitr少用org.Hs.eg.db数据库等</font>

<font style="color:rgba(0, 0, 0, 0.9);">⑤多用</font><font style="color:rgb(1, 1, 1);">Ensembl ID</font>

<font style="color:rgb(1, 1, 1);">⑥参考</font>[https://github.com/ixxmu/mp_duty/issues/1014](https://github.com/ixxmu/mp_duty/issues/1014)

[https://mp.weixin.qq.com/s/evF9kFNcjcJYYDB4vgTbow](https://mp.weixin.qq.com/s/evF9kFNcjcJYYDB4vgTbow)

R包影响比参考数据库影响比较小

![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1726106260781-fb1db12b-8300-4f80-8f3d-6dba44a7938e.png)

![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1726106267113-e8403ad6-a0ef-45cc-ab35-3f350f209bbe.png)

![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1726106305692-752fd08d-479d-4492-a632-c2c4053a151f.png)

**既然用ensembl多，那就用ensembl自家的biomaRT包，即3.3法推荐常用！！!**

<font style="color:rgb(1, 1, 1);">其他方法补充</font>[https://github.com/paulgeeleher/pRRophetic2/blob/master/pRRophetic/vignettes/prepareData.R](https://github.com/paulgeeleher/pRRophetic2/blob/master/pRRophetic/vignettes/prepareData.R)



![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1726105536271-fbe37ac2-4894-4b70-9043-d955db9025fa.png)

![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1726105305954-b3934408-7410-4266-b82a-d1b62287f4d3.png)

![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1726105566027-483d7272-24db-419e-8bf6-a40bfa359003.png)

![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1726105577973-694109b2-e847-4c41-a687-c44dcd8fc6f7.png)  
（2）要多少的匹配率才合适？

![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1726105553051-f349b066-d159-4b30-816b-c456899d5ec7.png)

参考[https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-229](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-229)

![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1726109217379-2ee05c41-1c22-4c00-b6e2-86f210237eba.png)

参考：[https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1779800/](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1779800/)

![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1726110156639-4a7ef633-94be-4bb5-89c3-5c1bf7131979.png)



这个问题到github issue、google、以及一些转换R包的issue、research gate、biostars上都没怎么找到，还找了一些生信英文期刊；

最后终于找到这么一点：

[https://www.researchgate.net/figure/Entrez-ID-to-gene-symbol-conversion-accuracy_tbl1_230831496](https://www.researchgate.net/figure/Entrez-ID-to-gene-symbol-conversion-accuracy_tbl1_230831496)

![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1726110855899-f6a932d6-ecea-41e1-89e7-531b4699737e.png)

![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1726110980622-7b500700-9856-476c-8c5f-2892f2f38a9b.png)

**从上面的table来看，个人认为最好是80%以上，当然越多越好！**

参考文献：[https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-229](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-13-229)

![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1726111196777-2639bbf6-d028-40df-ab2f-38ccdbc11b8d.png)

（3）一个基因有两个id等奇怪问题：

参考：[http://www.bio-info-trainee.com/8914.html](http://www.bio-info-trainee.com/8914.html)

![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1726106501751-8a72e25d-390e-4aa5-ad14-474af14b210b.png)

id的重复值问题：  
![](https://cdn.nlark.com/yuque/0/2024/png/33753661/1726106772879-aa0bd64d-db9b-4d3f-9df4-d692a258f286.png)

（4）尽量参考高质量的生信创作者：

生信技能树、生信菜鸟团、生信益站、生信宝典、生信星球+一些生物技术分析公司公众号，其他的个人博客少看、看多了误导性很强

