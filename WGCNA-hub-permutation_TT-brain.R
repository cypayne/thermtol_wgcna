### Adapted from Morgan et al 2020
### Permutations - are hub genes more likely to lie in 

## Random sample the full list of hub genes, each time counting the number of sampled
## genes that fall within a module. 
## Then determine whether the real number of hub genes within the module lies 
## outside of the 95% quantiles of the normal distribution.

library(dplyr)

setwd("~/Box/Schumer_lab_resources/Project_files/Thermal_tolerance_projects/Data/LTREB_qtl/final_rqtl-run_14I2021/LTREB-WGCNA_st07_onehot")

## Determine how many hub genes are transgressive (low, high, both) at 22C
## 15 total
geneset22c<-read.csv("~/Box/Schumer_lab_resources/Project_files/Thermal_tolerance_projects/Data/dge/brain_w-mito/transgressive-exp/TT-brain-22c-w-mito_DGE_lfc-shr_all.csv_with-annots.csv_with-F1info.csv",h=T,sep=",")
hub.genes_Kwithin<-read.csv('TT-brain-WGCNA-onehot_hubgenes_kWithin0.95_MM0.85.csv',header=T,sep=',')
hub.genes_KTotal<-read.csv('TT-brain-WGCNA-onehot_hubgenes_kTotal0.95_MM0.85.csv',header=T,sep=',')
hub.genes<-hub.genes_Kwithin$X
#hub.genes<-read.csv('TT-brain-WGCNA-onehot_hubgenes_kTotal0.95_MM0.85.csv',header=T,sep=',')
num_hub_genes = 600 # num hub genes using 0.95q cutoff

trans.both<-subset(geneset22c,subset=F1_transgress.high=="1"|F1_transgress.low=="1") #225
trans.both_hubgenes<-subset(subset(geneset22c,subset=F1_transgress.high=="1"|F1_transgress.low=="1"), Gene %in% hub.genes) #15
trans.low_hubgenes<-subset(subset(geneset22c,subset=F1_transgress.low=="1"), Gene %in% hub.genes) #12
trans.high_hubgenes<-subset(subset(geneset22c,subset=F1_transgress.high=="1"), Gene %in% hub.genes) #3

## Determine how many hub genes are transgressive (low, high, both) at 33C
## NONE
geneset33c<-read.csv("~/Box/Schumer_lab_resources/Project_files/Thermal_tolerance_projects/Data/dge/brain_w-mito/transgressive-exp/TT-brain-33c-w-mito_DGE_lfc-shr_all.csv_with-annots.csv_with-F1info.csv",h=T,sep=",")
trans.both<-subset(geneset33c,subset=F1_transgress.high=="1"|F1_transgress.low=="1") #83
trans.both_hubgenes<-subset(subset(geneset33c,subset=F1_transgress.high=="1"|F1_transgress.low=="1"), Gene %in% hub.genes) #0
trans.low_hubgenes<-subset(subset(geneset33c,subset=F1_transgress.low=="1"), Gene %in% hub.genes) #0
trans.high_hubgenes<-subset(subset(geneset33c,subset=F1_transgress.high=="1"), Gene %in% hub.genes) #0


## Permutations - expected transgressive hub genes
## 22C
# How many hub genes are expected to have transgressive f1 expression (both low and high)?
allGenes_withinTransgressive<-c()
for (i in 1:1000) {
  #randomly sample num_hub_genes rows from the GeneSet dataframe 
  random_genes<-sample_n(geneset22c, num_hub_genes)
  ran<-subset(random_genes,subset=F1_transgress.high=="1"|F1_transgress.low=="1")
  number_transgressive<-nrow(ran)
  allGenes_withinTransgressive<-append(allGenes_withinTransgressive,number_transgressive)
}
hist(allGenes_withinTransgressive)
quantile(allGenes_withinTransgressive,probs=c(0.9,0.95,0.99))
#Results: 22C:Kwithin=600 - 90% quantile - 10, 95% - 12, 99% - 14. 15 (the number of hub genes that are transgressive) is greater than the 95% quantile so this is significant.

# How many hub genes are expected to have low transgressive f1 expression?
allGenes_withLowTransExp<-c()
for (i in 1:1000) {
  #randomly sample num_hub_genes rows from the GeneSet dataframe
  random_genes<-sample_n(geneset22c, num_hub_genes)
  ran<-subset(random_genes,subset=F1_transgress.low=="1")
  number_transgressive<-nrow(ran)
  allGenes_withLowTransExp<-append(allGenes_withLowTransExp,number_transgressive)
}
hist(allGenes_withLowTransExp)
quantile(allGenes_withLowTransExp,probs=c(0.9,0.95,0.99))
#Results: 22C - 5% quantile - 1, 90% quantile - 6, 95% - 7, 99% - 8. 12 is greater than the 95% quantile so this is significant.

# How many hub genes are expected to have high transgressive f1 expression?
allGenes_withHighTransExp<-c()
for (i in 1:1000) {
  #randomly sample num_hub_genes rows from the GeneSet dataframe
  random_genes<-sample_n(geneset22c, num_hub_genes)
  ran<-subset(random_genes,subset=F1_transgress.high=="1")
  number_transgressive<-nrow(ran)
  allGenes_withHighTransExp<-append(allGenes_withHighTransExp,number_transgressive)
}
hist(allGenes_withHighTransExp)
quantile(allGenes_withHighTransExp,probs=c(0.05,0.9,0.95,0.99))
#Results: 22C - 5% quantile - 1, 90% quantile - 6, 95% quantile - 7, 99% - 8. 1 < 3 < 7 = not significant.


## 33C
# How many hub genes are expected to have transgressive f1 expression (both low and high)?
allGenes_withinTransgressive<-c()
for (i in 1:1000) {
  #randomly sample num_hub_genes rows from the GeneSet dataframe 
  random_genes<-sample_n(geneset33c, num_hub_genes)
  ran<-subset(random_genes,subset=F1_transgress.high=="1"|F1_transgress.low=="1")
  number_transgressive<-nrow(ran)
  allGenes_withinTransgressive<-append(allGenes_withinTransgressive,number_transgressive)
}
hist(allGenes_withinTransgressive)
quantile(allGenes_withinTransgressive,probs=c(0.01,0.05,0.9,0.95))
#Results: 22C:Kwithin=600 - 1% quantile - 0, 5% - 0, 90% - 5, 95% - 6. 0 (the number of hub genes that are transgressive at 33c) is not significant.


### WGCNA module permutations
# load all of the significant module names (sigMEs object)
load(file = "TT-brain-WGCNA_sigMEs.RData")
sigMEs

## Permutations - expected number of transgressive genes in a module
header="TT-brain-WGCNA-onehot_"
output_file=paste0(header,"num-trans-genes-in-ME_permutation-results.csv")
header_line = data.frame('module_name','num_MEgenes','num_MEtrans.both.22c','num_MEtrans.low.22c','num_MEtrans.high.22c','num_MEtrans.both.33c','num_MEtrans.low.33c','num_MEtrans.high.33c','22c_q0.01','22c_q0.05','22c_q0.1','22c_q0.9','22c_q0.95','22c_q0.99','33c_q0.01','33c_q0.05','33c_q0.1','33c_q0.9','33c_q0.95','33c_q0.99')
write.table(header_line,output_file,col.names=F,row.names=F,sep=",")
#modules_of_interest = c('MEpalevioletred3','MEpurple')
for( MEname in sigMEs ) { 
  MEname <- substring(MEname, 3)
  MEfile <- read.csv(paste0(header,MEname,"_genes.csv"),header=T)
  genes <- MEfile$x
  num_MEgenes <- length(genes)
  # number of transgressive genes in module
  MEtrans.both22c<-subset(subset(geneset22c,subset=F1_transgress.high=="1"|F1_transgress.low=="1"), Gene %in% genes)
  MEtrans.low22c<-subset(subset(geneset22c,subset=F1_transgress.low=="1"), Gene %in% genes)
  MEtrans.high22c<-subset(subset(geneset22c,subset=F1_transgress.high=="1"), Gene %in% genes)
  
  MEtrans.both33c<-subset(subset(geneset33c,subset=F1_transgress.high=="1"|F1_transgress.low=="1"), Gene %in% genes)
  MEtrans.low33c<-subset(subset(geneset33c,subset=F1_transgress.low=="1"), Gene %in% genes)
  MEtrans.high33c<-subset(subset(geneset33c,subset=F1_transgress.high=="1"), Gene %in% genes)
  
  row2out <- data.frame(MEname,num_MEgenes,dim(MEtrans.both22c)[1],dim(MEtrans.low22c)[1],dim(MEtrans.high22c)[1],dim(MEtrans.both33c)[1],dim(MEtrans.low33c)[1],dim(MEtrans.high33c)[1])
  
  # output the expression profiles for the transgressive genes
  write.csv(MEtrans.both22c,paste0(header,MEname,"_transF1-expression-22c.csv"))
  write.csv(MEtrans.both33c,paste0(header,MEname,"_transF1-expression-33c.csv"))
  
  # perms for 22c
  allGenes_withinTransgressive22c<-c()
  for (i in 1:1000) {
    #randomly sample the number of genes in module 
    random_genes<-sample_n(geneset22c, num_MEgenes)
    ran<-subset(random_genes,subset=F1_transgress.high=="1"|F1_transgress.low=="1")
    number_transgressive<-nrow(ran)
    allGenes_withinTransgressive22c<-append(allGenes_withinTransgressive22c,number_transgressive)
  }
  #hist(allGenes_withinTransgressive22c)
  quantiles22c<-quantile(allGenes_withinTransgressive22c,probs=c(0.01,0.05,0.1,0.9,0.95,0.99))
  row2out<-append(row2out,quantiles22c)
  
  # perms for 33c
  allGenes_withinTransgressive33c<-c()
  for (i in 1:1000) {
    #randomly sample the number of genes in module 
    random_genes<-sample_n(geneset33c, num_MEgenes)
    ran<-subset(random_genes,subset=F1_transgress.high=="1"|F1_transgress.low=="1")
    number_transgressive<-nrow(ran)
    allGenes_withinTransgressive33c<-append(allGenes_withinTransgressive33c,number_transgressive)
  }
  #hist(allGenes_withinTransgressive33c)
  quantiles33c<-quantile(allGenes_withinTransgressive33c,probs=c(0.01,0.05,0.1,0.9,0.95,0.99))
  row2out<-append(row2out,quantiles33c)
  
  # output results for each module to the output file
  write.table(row2out,output_file,col.names=F,row.names=F,sep=",",append=T)
  
}


## Permutations - expected number of hub genes in a module
# --> are there more highly connected genes in the thermtol modules?
# observed value: subsample the hub genes in the module
# simulations:
#   randomly sample genes from the entire gene set (22c -- 946 genes)
#   see how many of these randomly sampled genes fall inside the module
# use quantile function to determine the 0.9 and 0.95 percentile cutoffs, compare to observed value
hub.genes<-hub.genes_KTotal$X #705
header="TT-brain-WGCNA-onehot_"
output_file=paste0(header,"num-KTotal-hub-genes-in-ME_permutation-results.csv")
header_line = data.frame('module_name','num_MEgenes','num_MEhub.genes','q0.01','q0.05','q0.1','q0.9','q0.95','q0.99')
write.table(header_line,output_file,col.names=F,row.names=F,sep=",")
#modules_of_interest = c('MEpalevioletred3','MEpurple')
for( MEname in sigMEs ) { 
  MEname <- substring(MEname, 3)
  MEfile <- read.csv(paste0(header,MEname,"_genes.csv"),header=T)
  MEgenes <- MEfile$x
  num_MEgenes <- length(MEgenes)
  # number of hub genes in module
  MEhub.genes<-subset(subset(geneset22c, Gene %in% hub.genes), Gene %in% MEgenes)
  row2out <- data.frame(MEname,num_MEgenes,dim(MEhub.genes)[1])
  
  # output the expression profiles for the transgressive genes
#  write.csv(MEhub.genes,paste0(header,MEname,"_hub.gene-expression.csv"))
  
  # sample from either 22c or 33c, doesn't change hub gene status
  allGenes_withinHub<-c()
  for (i in 1:1000) {
    #randomly sample the number of genes in module
    random_genes<-sample_n(geneset22c, num_MEgenes)
    ran<-subset(random_genes,Gene %in% hub.genes)
    number_hub<-nrow(ran)
    allGenes_withinHub<-append(allGenes_withinHub,number_hub)
  }
  quantiles<-quantile(allGenes_withinHub,probs=c(0.01,0.05,0.1,0.9,0.95,0.99))
  row2out<-append(row2out,quantiles)
  
  # output results for each module to the output file
  write.table(row2out,output_file,col.names=F,row.names=F,sep=",",append=T)
  
}


## Permutations - how many transgressive hub genes do we expect?
## try both kWithin and kTotal
# hub.genes<-hub.genes_Kwithin$X
hub.genes<-hub.genes_KTotal$X #705
header="TT-brain-WGCNA-onehot_"
output_file=paste0(header,"num-trans-kTotal-hub-genes_permutation-results.csv")
header_line = data.frame('num_hub_genes','num_hub_trans.both.22c','num_hub_trans.low.22c','num_hub_trans.high.22c','num_hub_trans.both.33c','num_hub_trans.low.33c','num_hub_trans.high.33c','22c_q0.01','22c_q0.05','22c_q0.1','22c_q0.9','22c_q0.95','22c_q0.99','33c_q0.01','33c_q0.05','33c_q0.1','33c_q0.9','33c_q0.95','33c_q0.99')
write.table(header_line,output_file,col.names=F,row.names=F,sep=",")

# number of hub genes with transgression in F1s
hub.trans.both22c<-subset(subset(geneset22c,subset=F1_transgress.high=="1"|F1_transgress.low=="1"), Gene %in% hub.genes)
hub.trans.low22c<-subset(subset(geneset22c,subset=F1_transgress.low=="1"), Gene %in% hub.genes)
hub.trans.high22c<-subset(subset(geneset22c,subset=F1_transgress.high=="1"), Gene %in% hub.genes)
  
hub.trans.both33c<-subset(subset(geneset33c,subset=F1_transgress.high=="1"|F1_transgress.low=="1"), Gene %in% hub.genes)
hub.trans.low33c<-subset(subset(geneset33c,subset=F1_transgress.low=="1"), Gene %in% hub.genes)
hub.trans.high33c<-subset(subset(geneset33c,subset=F1_transgress.high=="1"), Gene %in% hub.genes)

row2out <- data.frame(num_hub_genes,dim(hub.trans.both22c)[1],dim(hub.trans.low22c)[1],dim(hub.trans.high22c)[1],dim(hub.trans.both33c)[1],dim(hub.trans.low33c)[1],dim(hub.trans.high33c)[1])
  
# output the expression profiles for the hub genes
write.csv(hub.trans.both22c,paste0(header,"_all-kTotal-hub-expression-22c.csv"))
write.csv(hub.trans.both33c,paste0(header,"_all-kTotal-hub-expression-33c.csv"))

# perms for 22c
allGenes_withinTransgressive22c<-c()
for (i in 1:1000) {
  #randomly sample the number of genes in module 
  random_genes<-sample_n(geneset22c, num_hub_genes)
  ran<-subset(random_genes,subset=F1_transgress.high=="1"|F1_transgress.low=="1")
  number_transgressive<-nrow(ran)
  allGenes_withinTransgressive22c<-append(allGenes_withinTransgressive22c,number_transgressive)
}
#hist(allGenes_withinTransgressive22c)
quantiles22c<-quantile(allGenes_withinTransgressive22c,probs=c(0.01,0.05,0.1,0.9,0.95,0.99))
row2out<-append(row2out,quantiles22c)
  
# perms for 33c
allGenes_withinTransgressive33c<-c()
for (i in 1:1000) {
  #randomly sample the number of genes in module 
  random_genes<-sample_n(geneset33c, num_hub_genes)
  ran<-subset(random_genes,subset=F1_transgress.high=="1"|F1_transgress.low=="1")
  number_transgressive<-nrow(ran)
  allGenes_withinTransgressive33c<-append(allGenes_withinTransgressive33c,number_transgressive)
}
#hist(allGenes_withinTransgressive33c)
quantiles33c<-quantile(allGenes_withinTransgressive33c,probs=c(0.01,0.05,0.1,0.9,0.95,0.99))
row2out<-append(row2out,quantiles33c)
  
# output results for each module to the output file
write.table(row2out,output_file,col.names=F,row.names=F,sep=",",append=T)
  