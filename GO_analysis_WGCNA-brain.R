## Gene Ontology enrichment analysis of WGCNA modules ##
# takes in deseq2 dge output for all genes to build universe
#                        and for sig genes to test enrichment

library("GOstats")
library("GSEABase")
library("biomaRt")

setwd("~/Box/Schumer_lab_resources/Project_files/Thermal_tolerance_projects/Data/LTREB_qtl/final_rqtl-run_14I2021/LTREB-WGCNA_st07_onehot/")
## TT brain transgressive
c22<-read.csv(file="~/Box/Schumer_lab_resources/Project_files/Thermal_tolerance_projects/Data/dge/brain_w-mito/transgressive-exp/TT-brain-22c-w-mito_DGE_lfc-shr_all.csv_with-annots.csv_with-F1info.csv",head=TRUE)
dgeres<-c22 # 19143 genes total, including mito
c33<-read.csv(file="~/Box/Schumer_lab_resources/Project_files/Thermal_tolerance_projects/Data/dge/brain_w-mito/transgressive-exp/TT-brain-33c-w-mito_DGE_lfc-shr_all.csv_with-annots.csv_with-F1info.csv",head=TRUE)
#dgeres<-c33

mart <- useMart(biomart = "ensembl", dataset = "xmaculatus_gene_ensembl", host="uswest.ensembl.org")
#attributes <- listAttributes(mart)
#attributes[1:50,]

# match go ids to ensembl gene ids
# c22 = 89701
results <- getBM(attributes = c("go_id","external_gene_name","ensembl_gene_id","kegg_enzyme"), filters=c("ensembl_gene_id"),values=dgeres$Gene, mart = mart)

# subset gene universe to only include genes with valid go ids and external gene names
gene_universe<-subset(results,nchar(results$go_id)>0 & nchar(results$external_gene_name) >0)
gene_universe$ensembl_id<-gene_universe[,3]
gene_universe[,3]<-as.numeric(as.factor(gene_universe[,3]))
gene_universe$Evidence<-rep("ISA",length(gene_universe[,3]))
colnames(gene_universe)<-c("frame.go_id","frame.gene_name","frame.gene_id","frame.KEGG","frame.ensembl","frame.Evidence")

goframeData <- data.frame(gene_universe$frame.go_id,gene_universe$frame.Evidence,gene_universe$frame.gene_id)

goFrame <- GOFrame(goframeData,organism="Xiphophorus")
goAllFrame <- GOAllFrame(goFrame)
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())

# c22 = 78787
universe <- goframeData$gene_universe.frame.gene_id

## GO enrichment

# TT brain transgressive
# to get immune genes not in misexp set: enrich for all !(misexpressed)
notmisexp.22c <- subset(c22, is.na(sig.F1_exp_profile))
dgesig<-notmisexp.22c

# 22c
head(dgeres)
trans.high.22c<-subset(c22, F1_transgress.high.111 == 1)
dgesig<-trans.high.22c
# trans.low<-subset(dgeres, F1_transgress.low.114 == 1)
# dgesig<-trans.low

# 33c
#head(dgeres)
trans.high.33c<-subset(c33, F1_transgress.high.57 == 1)
#dgesig <- trans.high.33c
#trans.low<-subset(dgeres, F1_transgress.low.26 == 1)
#dgesig<-trans.low

# trans.high at both 22c and 33c
both.high <- trans.high.22c[trans.high.22c$Gene %in% trans.high.33c$Gene,]
both.high 

#read in WGCNA module ids and match them to the gene_universe ids
header = "TT-brain-GO-wgcna_"
## choose modules of interest from those with at least one sig corr, groups of interesting modules below
#temp strongest corr: modules_of_interest = c('lightgreen','salmon','palevioletred3')
#temp only: modules_of_interest = c('black','darkorange','palevioletred3')
#notF1s: modules_of_interest = c('purple','red','magenta','floralwhite','skyblue','skyblue3','brown','midnightblue')
#notParents: modules_of_interest = c('salmon','royalblue','lightgreen')
#all: modules_of_interest = c('cyan')
#allbut1parent: modules_of_interest = c('paleturquoise','grey60')
#justF1s: modules_of_interest = c('brown4')

for(module in sigMEs) {
  module <- substring(module, 3)
  dgesig<-read.csv(file=paste0("./TT-brain-WGCNA-onehot_",module,"_genes.csv"),head=TRUE)
  genes_sig <- dgesig$x
  genes_match<-gene_universe[gene_universe$frame.ensembl %in% genes_sig,]
  genes_match_sig <- genes_match$frame.gene_id

  write.csv(genes_match,paste0(header,module,"_GO-gene-annots.csv"))

# set params for hyperGTest
# testDirection = "over" for overrepresented genes, = "under" for underrepresented
# interested in overrep genes
# three categories: 
#   cellular component (CC; where gene products are active)
#   molecular function (MF; the biological function of gene or gene product) 
#   biological process (BP; pathways or larger processes that multiple gene products involved in).
# Output: 
# ExpCount is the expected count and the Count is how many instances of that term were actually oberved 
# in your gene list while the Size is the number that could have been found in your gene list if every 
# instance had turned up. Values like the ExpCount and the Size are going to be affected by what is included 
# in the gene universe as well as by whether or not it was a conditional test.
  for(GO_type in c("CC","BP")) {
    params <- GSEAGOHyperGParams(name="Xiphophorus maculatus genes",
                             geneSetCollection=gsc,
                             geneIds = genes_match_sig,
                             universeGeneIds = universe,
                             ontology = GO_type,
                             pvalueCutoff = 0.95,
                             conditional = FALSE,
                             testDirection = "over")

    # get hypergeometric test results
    OverGtest <- hyperGTest(params)
    results_OverGtest<-summary(OverGtest)

    ## Output an additional column with the genes falling under each GO category
    if(GO_type == "BP") {
      # FOR BP
      res_GOBPID <- results_OverGtest$GOBPID 
      cats <- geneIdsByCategory(OverGtest) 
      # get list of lists: all geneids per GOBPID
      extract_geneids <- sapply(res_GOBPID, function(go) cats[[go]])
      # this gives dataframe of two columns: GOBPID, and the list of gene ids falling under it
      gobpid2geneid <- data.frame(GOBPID=names(extract_geneids), GO_GENEID_list=matrix(extract_geneids))
      gobpid2geneid[sapply(gobpid2geneid, is.list)] <- apply(gobpid2geneid[sapply(gobpid2geneid, is.list)], 1, function(x) paste(unlist(x), sep=", ", collapse=", "))
      # merge GO_GENEID_list with the rest of GO results
      go_output <- merge(results_OverGtest,gobpid2geneid, by.x="GOBPID", by.y="GOBPID", all=TRUE)
      write.table(go_output,paste0(header, module, '-GOresults-BP_w-geneids-pval0.95.tsv'), row.names = FALSE, quote=FALSE, sep="\t")
    } else if (GO_type == "CC") {
      # FOR CC
      res_GOCCID <- results_OverGtest$GOCCID # for CC
      cats <- geneIdsByCategory(OverGtest) 
      # get list of lists: all geneids per GOBPID
      extract_geneids <- sapply(res_GOCCID, function(go) cats[[go]])
      # this gives dataframe of two columns: GOBPID, and the list of gene ids falling under it
      goccid2geneid <- data.frame(GOCCID=names(extract_geneids), GO_GENEID_list=matrix(extract_geneids))
      goccid2geneid[sapply(goccid2geneid, is.list)] <- apply(goccid2geneid[sapply(goccid2geneid, is.list)], 1, function(x) paste(unlist(x), sep=", ", collapse=", "))
      # merge GO_GENEID_list with the rest of GO results
      go_output <- merge(results_OverGtest,goccid2geneid, by.x="GOCCID", by.y="GOCCID", all=TRUE)
      write.table(go_output,paste0(header, module, '-GOresults-CC_w-geneids-pval0.95.tsv'), row.names = FALSE, quote=FALSE, sep="\t")
    }
  }
}
# 33c
#write.csv(results_OverBP,'/Users/cypayne/Box/Schumer_lab_resources/Project_files/Thermal_tolerance_projects/Data/dge/brain_w-mito/transgressive-exp/TT-brain-33c-w-mito_DGE_lfc-shr_all-F1trans.high-GOresults-BP.csv')
#write.csv(results_OverBP,'/Users/cypayne/Box/Schumer_lab_resources/Project_files/Thermal_tolerance_projects/Data/dge/brain_w-mito/transgressive-exp/TT-brain-33c-w-mito_DGE_lfc-shr_all-F1trans.high-GOresults-CC.csv')
#write.csv(results_OverBP,'/Users/cypayne/Box/Schumer_lab_resources/Project_files/Thermal_tolerance_projects/Data/dge/brain_w-mito/transgressive-exp/TT-brain-33c-w-mito_DGE_lfc-shr_all-F1trans.high-GOresults-MF.csv')
#write.csv(results_OverBP,'/Users/cypayne/Box/Schumer_lab_resources/Project_files/Thermal_tolerance_projects/Data/dge/brain_w-mito/transgressive-exp/TT-brain-33c-w-mito_DGE_lfc-shr_all-F1trans.low-GOresults-BP.csv')
#write.csv(results_OverBP,'/Users/cypayne/Box/Schumer_lab_resources/Project_files/Thermal_tolerance_projects/Data/dge/brain_w-mito/transgressive-exp/TT-brain-33c-w-mito_DGE_lfc-shr_all-F1trans.low-GOresults-CC.csv')
#write.csv(results_OverBP,'/Users/cypayne/Box/Schumer_lab_resources/Project_files/Thermal_tolerance_projects/Data/dge/brain_w-mito/transgressive-exp/TT-brain-33c-w-mito_DGE_lfc-shr_all-F1trans.low-GOresults-MF.csv')


# # 22c
write.csv(results_OverBP,'/Users/cypayne/Box/Schumer_lab_resources/Project_files/Thermal_tolerance_projects/Data/dge/brain_w-mito/transgressive-exp/TT-brain-22c-w-mito_DGE_lfc-shr_all-F1trans.high-GOresults-BP-2.csv')
# #write.csv(results_OverBP,'/Users/cypayne/Box/Schumer_lab_resources/Project_files/Thermal_tolerance_projects/Data/dge/brain_w-mito/transgressive-exp/TT-brain-22c-w-mito_DGE_lfc-shr_all-F1trans.high-GOresults-CC.csv')
# write.csv(results_OverBP,'/Users/cypayne/Box/Schumer_lab_resources/Project_files/Thermal_tolerance_projects/Data/dge/brain_w-mito/transgressive-exp/TT-brain-22c-w-mito_DGE_lfc-shr_all-F1trans.high-GOresults-MF.csv')
# #write.csv(results_OverBP,'/Users/cypayne/Box/Schumer_lab_resources/Project_files/Thermal_tolerance_projects/Data/dge/brain_w-mito/transgressive-exp/TT-brain-22c-w-mito_DGE_lfc-shr_all-F1trans.low-GOresults-BP.csv')
# #write.csv(results_OverBP,'/Users/cypayne/Box/Schumer_lab_resources/Project_files/Thermal_tolerance_projects/Data/dge/brain_w-mito/transgressive-exp/TT-brain-22c-w-mito_DGE_lfc-shr_all-F1trans.low-GOresults-CC.csv')
# write.csv(results_OverBP,'/Users/cypayne/Box/Schumer_lab_resources/Project_files/Thermal_tolerance_projects/Data/dge/brain_w-mito/transgressive-exp/TT-brain-22c-w-mito_DGE_lfc-shr_all-F1trans.low-GOresults-MF.csv')


# match sig go ids to gene nums
#go2gene_matchup <- data.frame("GOBPID" = genes_match$frame.go_id, "GeneID" = genes_match$frame.ensembl)

#c33high_go2gene_matchup <- genes_match
#write.csv(c33high_go2gene_matchup,'/Users/cypayne/Box/Schumer_lab_resources/Project_files/Thermal_tolerance_projects/Data/dge/brain_w-mito/transgressive-exp/TT-brain-33c-w-mito_DGE_lfc-shr_all-F1trans.high-GOsiggenesmatch.csv')
# c33low_go2gene_matchup <- genes_match
# write.csv(c33low_go2gene_matchup,'/Users/cypayne/Box/Schumer_lab_resources/Project_files/Thermal_tolerance_projects/Data/dge/brain_w-mito/transgressive-exp/TT-brain-33c-w-mito_DGE_lfc-shr_all-F1trans.low-GOsiggenesmatch.csv')

#c22high_go2gene_matchup <- genes_match
#write.csv(c22high_go2gene_matchup,'/Users/cypayne/Box/Schumer_lab_resources/Project_files/Thermal_tolerance_projects/Data/dge/brain_w-mito/transgressive-exp/TT-brain-22c-w-mito_DGE_lfc-shr_all-F1trans.high-GOsiggenesmatch.csv')
c22low_go2gene_matchup <- genes_match
write.csv(c22low_go2gene_matchup,'/Users/cypayne/Box/Schumer_lab_resources/Project_files/Thermal_tolerance_projects/Data/dge/brain_w-mito/transgressive-exp/TT-brain-22c-w-mito_DGE_lfc-shr_all-F1trans.low-GOsiggenesmatch.csv')

# to look at genes that fall under a given go category
# in sig genes
# 48, 49, 108, 126, 130, 139
cats <- geneIdsByCategory(OverGtest)
c <- cats$`GO:0000422` # autophagy of mitochondrion:  atg7, rb1cc1
c <- cats$`GO:0061726` # mitochondrion disassembly:   atg7, rb1cc1
c <- cats$`GO:0006839` # mitochondrial transport:     SLC25A30, slc25a38b, chchd4a
c <- cats$`GO:0045041` # protein import into mitochondrial intermembrane space: chchd4a
c <- cats$`GO:0000423` # mitophagy: atg7
c <- cats$`GO:0070125` # mitochondrial translational elongation: gfm1
genes_match[genes_match$frame.gene_id %in% c,]

#library(plyr)
#gobpid2geneid <- plyr::ldply(extract_geneids, cbind) # gives two column df, one-to-one GO:geneid map
# now add xmac gene id to this df


# Look at dif in GO enrichment between pseudoalign refs
mal <- read.csv("/Users/cypayne/Desktop/Schumer_lab/ThermTol/dge_TT/TT-brain-mal33cV22c-padj0.1-GOresults.csv", header=TRUE)
bir <- read.csv("/Users/cypayne/Desktop/Schumer_lab/ThermTol/dge_TT/TT-brain-bir33cV22c-padj0.1-GOresults.csv", header=TRUE)

intersection <- intersect(bir$Term,mal$Term)
dim(bir) #161
dim(mal) #154
length(intersection) #84

bir_overlap <- subset(bir,Term %in% intersection)
bir_nolap <- subset(bir,!(Term %in% intersection)) #77

mal_overlap <- subset(mal,Term %in% intersection)
mal_nolap <- subset(mal,!(Term %in% intersection)) #70

write.csv(bir_overlap, "/Users/cypayne/Desktop/Schumer_lab/ThermTol/dge_TT/TT-brain-bir33cV22c-padj0.1-GOresults.bir-overlap.csv")
write.csv(mal_overlap, "/Users/cypayne/Desktop/Schumer_lab/ThermTol/dge_TT/TT-brain-mal33cV22c-padj0.1-GOresults.mal-overlap.csv")
write.csv(bir_nolap, "/Users/cypayne/Desktop/Schumer_lab/ThermTol/dge_TT/TT-brain-bir33cV22c-padj0.1-GOresults.bir-nolap.csv")
write.csv(mal_nolap, "/Users/cypayne/Desktop/Schumer_lab/ThermTol/dge_TT/TT-brain-mal33cV22c-padj0.1-GOresults.mal-nolap.csv")
