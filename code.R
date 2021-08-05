
pdata <- read.table('clinical.txt', header=TRUE, row.names=1, sep='\t')
dmat <- read.table('ESCA_lncRNAs_counts.txt', header=TRUE, row.names=1, sep='\t', check.names=F)


# 处理多个探针对应同一个基因
library('hgu95av2.db')
library(CLL)
data(sCLLex)
sCLLex=sCLLex[,1:8] ## 样本太多，我就取前面8个
group_list=sCLLex$Disease
exprSet=exprs(sCLLex)
exprSet=as.data.frame(exprSet)
exprSet$probe_id=rownames(exprSet)
head(exprSet)
probe2symbol=toTable(hgu95av2SYMBOL)
dat=merge(exprSet,probe2symbol,by='probe_id')
results=t(sapply(split(dat,dat$symbol),function(x) colMeans(x[,1ncol(x)-1)])))
############################################################################################

ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
		(qx[6]-qx[1] > 50 && qx[2] > 0) ||
		(qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
	exprs(gset) <- log2(ex) }



dmat <- round(dmat)


group <- factor(pdata$sample_type.samples)
design <- model.matrix(~group)
rownames(design) <- colnames(dmat)
colnames(design) <- levels(group)


DGElist <- DGEList(counts=dmat, group=group)

keep_gene <- rowSums( cpm(DGElist) > 10 )
table(keep_gene)
DGElist <- DGElist[ keep_gene, , keep.lib.sizes = FALSE ]

# 标准化
y <- calcNormFactors(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)

# 差异分析 exactTest
et <- exactTest(y,pair=c("normal","tumor"))
topTags(et)

# 差异分析 quasi-likelihood F-tests
fit <- glmQLFit(ylnc, design)
delnc <- glmQLFTest(fit, coef=2)
DElnc <- topTags(delnc, n=Inf, sort.by = "logFC", adjust.method = "BH")$table




# DGElist <- calcNormFactors( DGElist )
# DGElist <- estimateGLMCommonDisp(DGElist, design)
# DGElist <- estimateGLMTrendedDisp(DGElist, design)
# DGElist <- estimateGLMTagwiseDisp(DGElist, design)


# fit <- glmFit(DGElist, design)
# results <- glmLRT(fit, coef=2) 
# nrDEG_edgeR <- topTags(results, n = nrow(DGElist))
# nrDEG_edgeR <- as.data.frame(nrDEG_edgeR)
# head(nrDEG_edgeR)
# 
# dge.sig <- nrDEG_edgeR[nrDEG_edgeR$FDR < 0.05 & abs(nrDEG_edgeR$logFC) >= 2,]
# 
# norm.counts <- cpm(DGElist, log=TRUE)
# norm.counts.sig <- norm.counts[rownames(dge.sig),]


################################################################
# limma 多组比较

design1 <- model.matrix(~0+factor(pheno$group1))
rownames(design1) <- rownames(pheno)
colnames(design1) <- levels(factor(pheno$group1))
contrast.matrix1 <- makeContrasts(OMI-Control, OFI-Control, levels = design1)

fit_mrna1 <- lmFit(expr.mrna, design1)
fit2_mrna1 <- contrasts.fit(fit_mrna1, contrast.matrix1)
fit2_mrna1 <- eBayes(fit2_mrna1)

head(fit2_mrna1$coefficients, 2)

nrDEG_OMI_CTL <- topTable(fit2_mrna1, coef=1, adjust='fdr', n=Inf)
nrDEG_OFI_CTL <- topTable(fit2_mrna1, coef=2, adjust='fdr', n=Inf)


#################################################################

dds <- DESeqDataSetFromMatrix(countData=dataset, colData=pdata, design= ~sample_type.samples)
dds$sample_type.samples <- relevel(dds$sample_type.samples, ref="Solid_Tissue_Normal")

dds <- DESeq(dds)
res <- results(dds)

resSig <-  subset(res, pvalue < 0.05 & abs(log2FoldChange) > 2)

vsdSig <- vsddf[rownames(resSig),]


anno_col <- pdata[,'sample_type.samples', drop=FALSE]
col = colorRampPalette(c("lightblue", "yellow", "orange", "red"),bias=3)(300)
pheatmap(vsdSig, annotation_col=anno_col, col=col, show_rownames=FALSE, show_colnames=FALSE)


myfun <- function(x) {
  cox.res <- coxph(Surv(pheno$OS.time, pheno$OS) ~ x)
  summcph <- summary(cox.res)
  p.value<-signif(summcph$wald["pvalue"], digits=2)
    beta<-signif(summcph$coef[1], digits=2);#coeficient beta
    HR <-signif(summcph$coef[2], digits=2);#exp(beta)
    HR.confint.lower <- signif(summcph$conf.int[,"lower .95"], 2)
    HR.confint.upper <- signif(summcph$conf.int[,"upper .95"],2)
    HR <- paste0(HR, " (", 
                 HR.confint.lower, "-", HR.confint.upper, ")")
    res<-c(beta, HR, p.value)
    names(res)<-c("beta", "HR (95% CI for HR)", 
          "p.value")
    return(res)
}

cox_res <- apply(norm.counts, 1, myfun)
cox.res <- as.data.frame(t(cox_res))
sig_res <- cox_res[cox_res < 0.05]


library(survival)

outcome <- c("survival.vector")


## The lines below should not need modification.

## Create list of models
list.of.models <- lapply(seq_along((predictors)), function(n) {

    left.hand.side  <- outcome
    right.hand.side <- apply(X = combn(predictors, n), MARGIN = 2, paste, collapse = " + ")

    paste(left.hand.side, right.hand.side, sep = "  ~  ")
})

## Convert to a vector
vector.of.models <- unlist(list.of.models)

vector.of.models <- c(list.of.models[[8]])

## Fit coxph to all models
list.of.fits <- lapply(vector.of.models, function(x) {

    formula    <- as.formula(x)
    fit        <- coxph(formula, data = dfcox)
    result.AIC <- extractAIC(fit)

    data.frame(num.predictors = result.AIC[1],
               AIC            = result.AIC[2],
               model          = x)
})

## Collapse to a data frame
result <- do.call(rbind, list.of.fits)


## Sort and print
library(doBy)
orderBy(~ AIC, result)


# Venn Plot ggvenn
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 4
  )


# stacked barplot with p-value
df <- data.frame(type=rep(c('Non-thyroidits', 'Thyroidits'), each=2),
	genotype=rep(c('Wt', 'Mut'), 2),
	num=c(345,4,59,4))

ggplot(data=df, aes(x=type, y=num, fill=genotype)) +
  geom_bar(stat="identity", width = 0.6) +
  labs(title='RET mutation frequency in thyroidits', x='', y='No. of samples') +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_signif(annotations = 'p-value=0.022', y_position = 365 ,xmin="Non-thyroidits", xmax="Thyroidits")

# boxplot 
p <- ggboxplot(tmp2, x='group', y='ILMN_1691714', color='group', add='jitter')
p + stat_compare_means(method='t.test')


##### 两组基因表达dotplot
ggplot(df_tfs_expr, aes(x=sample_type.samples, y=FCER2, fill=sample_type.samples)) + 
    geom_boxplot(notch = TRUE)+
	geom_dotplot(binaxis='y', stackdir='center',stackratio=1.5, dotsize=0.5) +
    stat_compare_means(method='t.test')

# 蜂群图加箱型图
beeswarm(mRNAsi ~ Tissue, data = pdata, pch = 16,
         method = 'swarm', pwcol = as.numeric(Tissue), xlab = '', 
         ylab = 'Follow-up time (months)', 
         labels = c('Normal', 'Tumor'))
boxplot(mRNAsi ~ Tissue, data = pdata, add = T,
        names = c("",""), col="#0000ff22")


# Visualize: Specify the comparisons you want
my_comparisons <- list( c("0.5", "1"), c("1", "2"), c("0.5", "2") )
ggboxplot(ToothGrowth, x = "dose", y = "len",
          color = "dose", palette = "jco")+ 
  stat_compare_means(comparisons = my_comparisons, label.y = c(18, 16, 17))+ # Add pairwise comparisons p-value
  stat_compare_means(label.y = 50)     # Add global p-value

#蜂群图，添加一些参数美化
p <- ggplot(dat, aes(group, RNF25, color = group)) +
geom_beeswarm(cex = 1.5, show.legend = FALSE)

# correlation plot
pdf('mir944-GPR31.cor.pdf')
ggplot(data=temp, aes(x=mir_944, y=GPR31))+
    geom_point(color="blue")+
    stat_smooth(method="lm")+
    ylim(1, 5)+
    labs(x='miR-944', y='GPR31')+
    ggtitle('Relationship between miR-944 and GPR31\nCor=-0.28(p-value=1.43e-7)')+
    theme(plot.title = element_text(hjust = 0.5))
dev.off()

ggplot(data=pheno, aes(x=hypoxia, y=mRNAsi))+
    geom_point(color="grey", size=3)+
    geom_smooth(method='lm',se=T) +
    stat_cor(method = "pearson", label.x = 2.4, label.y = 0.9, size=6)+ # ggpubr add coefficient
    labs(x='Hypoxia Z-Score', y='mRNAsi')+
    theme(plot.title = element_text(hjust = 0.5),
      panel.background = element_blank(),
      axis.line = element_line(color="black"),
      axis.text.x = element_text(size=12),
      axis.title.y = element_text(size=16),
      axis.title.x = element_text(size=16),
      axis.text.y = element_text(size=12),
      legend.title = element_text(size=14),
      legend.text = element_text(size=12)
         )


# 分组小提琴图
ggplot(plot.info, aes(x=variable, y=value, fill=condition)) + 
    geom_violin(trim = T,position = position_dodge(0.7)) +
     geom_boxplot(width = 0.2, position = position_dodge(0.7))+ # position 控制组间间距
    labs(y='Relative expression level', x='') +
    scale_fill_discrete(name='Condition') +
    theme(axis.text.x = element_text(size=12, angle = 90),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=16),
          legend.title = element_text(size=14),
          legend.text = element_text(size=12)
         )
ggsave('result/hubs_boxplot.pdf', width = 8, height = 6)

# 处理多个探针对应同一个基因
result <- aggregate(x=expr[,1:18], by=list(expr$symbol), FUN=mean, na.rm=T)


# volcano plot
pdf('DEanalysis/volcano.pdf')  
EnhancedVolcano(nrDEG,
    lab = '',
    x = 'logFC',
    y = 'P.Value',
    pCutoff = 0.05, FCcutoff = 1, colAlpha = 0.6)
dev.off()

tiff('DEG/MAplot.tiff', width = 6, height = 6, units = 'in', res=800)
plotMA(resLFC, ylim=c(-2,2))
dev.off()

# custom volcano plot
library(ggrepel)

df.cor$threshold = as.factor(ifelse(df.cor$pval < 0.05 & abs(df.cor$coef) >= 0.2, ifelse(df.cor$coef > 0 ,'Positive','Negtive'),'NoSignifi'))

ggplot(data = df.cor, aes(x = coef, y = -log10(pval), colour=threshold)) +
  geom_point(alpha=0.9, size=1.5) +
  scale_color_manual(values=c("black","red")) +
  xlim(-1, 1) +
  labs(x='Pearson correlation coefficient (Pearson test)', y='-log10(pvalue)', title='SDC2') +
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=14),
        axis.title.x = element_text(size=14),
          legend.title = element_text(size=14),
          legend.text = element_text(size=12),
        plot.title = element_text(hjust = 0.5, size=16)
         ) +
    geom_text_repel(
        data = subset(df.cor, df.cor$pval < 0.05 & abs(df.cor$coef) >= 0.2),
        aes(label = drug),
        size = 3,
        box.padding = unit(0, "lines"),
        point.padding = unit(0, "lines"), segment.color = "black", show.legend = FALSE )
ggsave('result/SDC2_IC50_corr.pdf')



anno1 <- select(hgu133a.db, keys=rownames(set1),columns=c('SYMBOL'),keytype='PROBEID')


# survival plot with custom annotation
ggsurv <- ggsurvplot(sfit, pval=F, data=data)

ggsurv$plot <- ggsurv$plot +
  ggplot2::annotate(
    "text",
    x = Inf, y = Inf,
    vjust = 1, hjust = 1,
    label = "n(High-risk) = 60 \n n(Low-risk) = 59 \n p < 0.001",
    size = 5
  )

ggsurv


##################################################################################
#############################  甲基化芯片处理流程  #################################
library(data.table)
library(ChAMP)


# load 450k beta matrix
a <- fread('rawdata//TCGA-ESCA.methylation450.tsv.gz', data.table=FALSE)
rownames(a)=a[,1]
a=a[,-1]

beta=as.matrix(a)
beta=impute.knn(beta)
betaData=beta$data
betaData=betaData+0.00001
a=betaData

myLoad=champ.filter(beta = a,pd = pdata)
save(myLoad,file = 'methy.step1-output.Rdata')

myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=5)
dim(myNorm) 
pD=myLoad$pd
save(myNorm,pD,file = 'step2-champ_myNorm.Rdata')

# DMP
group_list=myLoad$pd$group
myDMP <- champ.DMP(beta = myNorm,pheno=group_list)

# probe features
data(probe.features)
probe.features
head(probe.features)



##################################################################################
############################### CiberSort 结果可视化 ##############################

# 将宽数据转换成长数据
plot.info2 <- melt(ciber, id=c('Input.Sample', 'Tissue'), measure=names(ciber)[2:21])


# 堆积柱形图
ggplot(plot.info3, aes(x = sample, y =value, fill =variable))+
	geom_bar(position='stack',stat="identity") +
	xlab("Sample ID") + 
	ylab("Composition") + 
	labs(fill="Cell type") +  # 更改legend标题
	theme(axis.text.x = element_text(angle=90), 
	legend.text=element_text(size=8),
	legend.key.size=unit(0.5, 'cm')) +  # 设置legend方块和文字大小
	guides(fill=guide_legend(ncol=1)) #将legend设置成1列

# 分组小提琴图
ggplot(plot.info2, aes(x=variable, y=value, fill=conditioin)) + 
  # geom_boxplot(width=0.3, position=position_dodge(0.9)) +
      geom_violin(trim = T,position = position_dodge(0.7)) +
     geom_boxplot(width = 0.2, position = position_dodge(0.7))+ # position 控制组间间距
  xlab("") + 
  ylab("Relative expression") +
    labs(fill='Condition') +
  theme(axis.text.x = element_text(size=12, angle=90),
          axis.text.y = element_text(size=12),
         axis.title.y = element_text(size=14),
          legend.title = element_text(size=14),
          legend.text = element_text(size=12)
         ) +
  scale_fill_manual(values = c("#AED4D5","#F9CC88","#FD475D")) +
  stat_compare_means(aes(group = conditioin), 
    label = "p.signif", 
    method="t.test", hide.ns=T)


##################################################################################
#################### 分析GEO 中illumina beadchip non-normalized 数据###############

# GSE60436_non-normalized_data.txt.gz 为例

library(limma)
x <- read.ilmn('GSE60436_non-normalized_data.txt.gz',expr="SAMPLE ",probeid="ID_REF")
y <- neqc(x)

boxplot(y$E, outline=FALSE)

##################################################################################
#################### clusterProfiler and GOplot 可视化 ###########################

ego.all <- enrichGO(gene       = green.df$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05)

egox <- setReadable(ego.all, 'org.Hs.eg.db', 'ENTREZID')

GO <- egox[,c(1,2,3,9,7)]
GO$geneID <- str_replace_all(GO$geneID,"/",",") ### 修改geneID这一列的分隔符号
names(GO)=c("Category","ID","Term","Genes","adj_pval")

genedata = data.frame(ID=modProbes,logFC=rnorm(length(modProbes),mean = 0,sd=2))
circ <- circle_dat(GO,genedata)

reduced_circ <- reduce_overlap(circ, overlap = 0.9)
# labels值超过数据的上线会报错，调小
GOBubble(reduced_circ, labels = 5)


######################################################################################
############################# cox regression #########################################

covariates <- c("age", "sex",  "ph.karno", "ph.ecog", "wt.loss")
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(time, status)~', x)))
                        
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = lung)})
# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                          x <- summary(x)
                          p.value<-signif(x$wald["pvalue"], digits=2)
                          wald.test<-signif(x$wald["test"], digits=2)
                          beta<-signif(x$coef[1], digits=2);#coeficient beta
                          HR <-signif(x$coef[2], digits=2);#exp(beta)
                          HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                          HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                          HR <- paste0(HR, " (", 
                                       HR.confint.lower, "-", HR.confint.upper, ")")
                          res<-c(beta, HR, wald.test, p.value)
                          names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                        "p.value")
                          return(res)
                          #return(exp(cbind(coef(x),confint(x))))
                         })
res <- t(as.data.frame(univ_results, check.names = FALSE))
as.data.frame(res)


###############################################################################################
# forestplot

forestplot(data[1:3], data[4:6], graph.pos=3, 
  vertices=TRUE, zero=1,
  boxsize = 0.4, #设置点估计的方形大小
  lineheight = unit(8,'mm'),#设置图形中的行距
  colgap = unit(2,'mm'),#设置图形中的列间距
  lwd.zero = 2,#设置参考线的粗细
  lwd.ci = 2,#设置区间估计线的粗细
  col=fpColors(box='#458B00',summary="#8B008B",lines = 'black',zero = '#7AC5CD'),
  #使用fpColors()函数定义图形元素的颜色，从左至右分别对应点估计方形，汇总值，区间估计线，参考线
  lwd.xaxis=2,#设置X轴线的粗细
  lty.ci = "solid",
  )

#############################################################################################

# 3d PCA plot 
# https://r-graphics.org/recipe-miscgraph-3d-save
p <- pca3d(pca_res, group=group, show.ellipses=F, ellipse.ci=0.75, show.plane=T,
  show.group.labels=F, show.centroids=F, show.shadows=F, show.scale=T, col=colors,
  axe.titles=c('PC1(44.28%)','PC2(25.87%)','PC3(13.20%)'))
p

legend3d("topleft", p$groups, col=c('red','green','blue'), pch = p$pch, cex=1, inset=c(0.02,0.2))

# 保存图像
rgl.snapshot('3dplot.png', fmt = 'png')
rgl.postscript('3dplot.pdf', fmt = 'pdf')

rgl.postscript('3dplot.ps', fmt = 'ps')


#############################################################################################
#####################################  ROC curve ############################################
library(survivalROC)


cutoff = 365*3

rs.3 <- survivalROC(Stime=df_rs$OS.time, status=df_rs$OS, 
  predict.time=cutoff, method="KM", marker=df_rs$RS)


plot(rs.1$FP, rs.1$TP, ## x=FP,y=TP
     type="l",col="black", lwd=2,
     xlim=c(0,1), ylim=c(0,1),   
     xlab=paste( "False Positive Rate"), ##连接
     ylab="True Positive Rate",
     main="Time dependent ROC")
abline(0,1,col="black", lwd=1.6)

lines(rs.3$FP, rs.3$TP, type="l",lwd=2,col="blue",xlim=c(0,1), ylim=c(0,1))
lines(rs.5$FP, rs.5$TP, type="l",lwd=2,col="red",xlim=c(0,1), ylim=c(0,1))

legend(0.67,0.1,c(paste("AUC at 1 year ",round(rs.1$AUC,3)),
                 paste("AUC at 3 year ",round(rs.3$AUC,3)),
      paste("AUC at 5 year ",round(rs.5$AUC,3))),
                 x.intersp=1, y.intersp=0.8,
                 lty= 1 ,lwd= 2,col=c("black", "blue","red"),
                 bty = "n",# bty框的类型
                 seg.len=1,cex=0.9)# 


# 二分类样本ROC曲线
library("pROC")

auc1 <- roc(group~SLC25A40, data=data, smooth=FALSE)

plot(auc1, print.auc=TRUE, col="blue", main="ROC of SLC25A40",
  identity.lty=2, identity.lwd=1)


###################################################################################
############################## 预后模型散点图  #####################################
p <- ggplot(data, aes(x=index, y=RS, color=riskScore)) + 
    geom_point() +
  labs(x='Patients (increasing risk score)', y='Risk score') +
  geom_hline(aes(yintercept=1.3042), colour="#990000", linetype="dashed") +
  geom_vline(aes(xintercept=244), colour="#990000", linetype="dashed") +
  theme_bw() +
  theme(axis.title.x = element_text(size=16),
     axis.text.x = element_text(size=12),
         axis.title.y = element_text(size=16),
     axis.text.y = element_text(size=12),
          legend.title = element_text(size=14),
          legend.text = element_text(size=12)
         )
p

p <- ggplot(data, aes(x=index, y=time, color=BCR)) + 
    geom_point() +
  labs(x='Patients (increasing risk score)', y='Time to BCR (year)') +
  #geom_hline(aes(yintercept=1.3042), colour="#990000", linetype="dashed") +
  geom_vline(aes(xintercept=244), colour="#990000", linetype="dashed") +
  theme_bw() +
  theme(axis.title.x = element_text(size=16),
     axis.text.x = element_text(size=12),
         axis.title.y = element_text(size=16),
     axis.text.y = element_text(size=12),
          legend.title = element_text(size=14),
          legend.text = element_text(size=12)
         ) +
   scale_color_manual(values=c('#00BFC4', '#F8766D'))
p



## pheatmap customise group colors
colors <- c('#00DAE0', '#FF9289')
names(colors) <- c('Low', 'High')
mycolors <- list(Risk=colors)
pheatmap(data.plot, annotation_col = anno_col, show_colnames = FALSE, cluster_cols = FALSE, annotation_colors = mycolors, )



#####################################################################################
##############################  nomogram  列线图  ####################################
library(rms)

# 打包数据
dd<-datadist(data)

options(datadist="dd")

res.cox2 <- cph(Surv(OS.time, OS) ~ Age+Gender+Metastasis+Primary.site+RiskScore, data = data,
      x=TRUE, y=TRUE, surv=TRUE)


surv <- Survival(res.cox2) # 构建生存概率函数
# function(x) surv(365, x)#计算1年事件发生概率
# function(x) surv(1825, x)#计算5年事件发生概率

nom.cox <- nomogram(res.cox2,
                    fun=list(function(x) surv(365, x),
             function(x) surv(1095, x),
                             function(x) surv(1825, x)),
                    funlabel=c("1-year Survival Probability",
          "3-year Survival Probability", 
          "5-year Survival Probability"),
                    maxscale=100,
                    fun.at=c(0.01,seq(0.1,0.9,by=0.2),0.95,0.99))

#maxscale参数指定最高分数，一般设置为100或者10分
#fun.at设置生存率的刻度
# 绘制列线图
plot(nom.cox)


# 绘制calibration 图
f2 <- cph(Surv(OS.time, OS) ~ Metastasis+RiskScore, data = data,
    dist='lognormal', x=TRUE, y=TRUE, surv = T,na.action=na.delete,time.inc = 1825) 

#参数m=30表示每组30个样本进行重复计算, m约等于患者总数/3 （** 重要 **）
cal2<-calibrate(f2, cmethod="KM", method="boot",u=1825, m=30, B=100) 

#pdf("calibration_5y.pdf",width = 8,height = 8)
plot(cal2,
     lwd = 2,#error bar的粗细
     lty = 1,#error bar的类型，可以是0-6
     errbar.col = c("#2166AC"), #error bar的颜色
     xlim = c(0,1),ylim= c(0,1),
     main="Calibration Curve (3 year)",
     xlab = "Nomogram-prediced probability of 5-year survival",
  ylab = "Actual probability of 5-year survival",
     cex.lab=1.2, cex.axis=1, cex.main=1.2, cex.sub=0.6) #字的大小

lines(cal2[,c('mean.predicted',"KM")], 
      type = 'b', #连线的类型，可以是"p","b","o"
      lwd = 2, #连线的粗细
      pch = 16, #点的形状，可以是0-20
      col = c("red")) #连线的颜色
#mtext("")

box(lwd = 1) #边框粗细

abline(0,1,lty = 3, #对角线为虚线
       lwd = 2, #对角线的粗细
       col = c("#224444")#对角线的颜色
       )
