# 写在前面
因为我研究的物种比较小众，很多注释不完全，R包AnnotationHub中也没有对应信息，所以无法使用公共数据库进行kegg富集分析。所以自己尝试使用KAAS造一个自己的基因集，然后再进行使用Y叔的clusterProfiler进行富集分析。我觉得这样的好处是更和自己的物种相贴切，不会有一些pathway自己物种中没有但是公共库中存在的情况（当然，也有可能应该是有这个pathway的，kegg没有注释到）。
# 我的大体思路
1. 使用KEGG注释网站KAAS将自己的序列对比到KEGG数据库的中，得到基因与功能蛋白（K）的关系。
2. 使用KO数据库的mapping功能，将功能蛋白（K）与pathway（ko）对应上，并得到pathway的注释信息。
3. 得到gene与pathway（ko）的对应关系。
4. 利用clusterProfiler的enricher功能分析基因富集情况。
## 1. KAAS自动注释
具体使用请参考简书这篇文章
[如何使用KAAS进行KEGG注释](https://www.jianshu.com/p/289e91670f43)

![KAAS](https://upload-images.jianshu.io/upload_images/15318746-d15fe5f000b59d26.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

![KAAS程序](https://upload-images.jianshu.io/upload_images/15318746-5db9ea68b865ec73.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

![i基因集选择](https://upload-images.jianshu.io/upload_images/15318746-747e247a72ed76cf.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

当分析完成后你邮件会收到一个网址，打开网址得到类似这个网页

![结果](https://upload-images.jianshu.io/upload_images/15318746-29899e6739eff807.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

html里有download KO list选项，里边就是你的基因和功能蛋白（K）的关系，没有的可能是没有和数据库中的比对到。

![基因和功能蛋白（K）的关系](https://upload-images.jianshu.io/upload_images/15318746-09278a83f2bc1950.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

![ko_list](https://upload-images.jianshu.io/upload_images/15318746-1d8d9eadf80164b5.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

## 2. 功能蛋白K与pathway（ko）对应

将pathway（ko）提取出来，放入[https://www.genome.jp/kegg/ko.html](https://www.genome.jp/kegg/ko.html)中，将K number填入，单击map pathway。
![map_pathway](https://upload-images.jianshu.io/upload_images/15318746-fe031cb397106342.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

点击show matched object 可以获得该pathway（ko）下的功能蛋白K信息。

![该pathway（ko）下的功能蛋白K信息](https://upload-images.jianshu.io/upload_images/15318746-a06a81677469ab91.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

将此信息整理成K2ko的形式，即每个功能蛋白K对应的pathway（ko）和pathway（ko）的注释信息term2name。（Excel和R语言均可处理，我的方法比较笨，就不在这里讲了）。

![K2ko](https://upload-images.jianshu.io/upload_images/15318746-679e000d600b4544.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

![term2name](https://upload-images.jianshu.io/upload_images/15318746-1201499b5f70cf41.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

## 3. 得到gene与pathway（ko）的对应关系

再将上一步得到的gene和K关系的结果准备好，使用R语言merge函数整理得到pathway（ko）和gene的对应信息term2gene。


![K2gene](https://upload-images.jianshu.io/upload_images/15318746-582b5b1d6d1370fb.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

``````
gene2ko=merge(k2gene,K2ko,by="K")
write.table(gene2ko,"gene2ko.tab",row.names = F,sep = "\t")
``````

![term2gene](https://upload-images.jianshu.io/upload_images/15318746-132be6143bb0d0fa.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

## 4. 利用clusterProfiler的enricher功能分析基因富集情况。
之后利用clusterProfiler包进行一些常规分析。

```````````
library("clusterProfiler")
# 导入基因列表
gene <- read.csv("test_kegg_gene.txt",header = F,sep=",")
gene <- as.factor(gene$V1)
# 导入注释文件
term2gene <- read.csv("./kegg/ko2gene.csv",header=T,sep=",")
term2name <- read.csv("term2keggName.csv",header=F,sep=",")
# 富集分析
x <- enricher(gene,TERM2GENE=term2gene,TERM2NAME=term2name,pvalueCutoff = 0.05, 
               pAdjustMethod = "BH",qvalueCutoff = 0.2) head(x)
# 绘制条形图
barplot(x)
# 绘制气泡图
dotplot(x)
`````````````````````

![条形图](https://upload-images.jianshu.io/upload_images/15318746-225febd4aa85fb55.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

![气泡图](https://upload-images.jianshu.io/upload_images/15318746-5038e81cbb4074e6.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)
