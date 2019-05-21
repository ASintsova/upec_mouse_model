
#library(topGO)
#library("org.EcK12.eg.db")
#library(graph)
#library("graph")
#library("Rgraphviz")
#library("SparseM")
#library("GOFunction")
library(ROntoTools)
kpg <- keggPathwayGraphs("eco", relPercThresh = 0.0, updateCache = TRUE, verbose = TRUE)
kpg <- setEdgeWeights(kpg, edgeTypeAttr = "subtype",
                      edgeWeightByType = list(activation = 1, inhibition = -1,
                                              expression = 1, repression = -1),
                      defaultWeight = 0)

core <- read.csv("/Users/annasintsova/git_repos/HUTI-RNAseq/results/differential_expression_analysis/best_strains_DEseq.csv", header = TRUE)
colnames(core)[1] <- "bnum"
core <- core[complete.cases(core),]
kpn <- keggPathwayNames("eco")
fc <- core$log2FoldChange[core$padj <= 0.01] # fold change
names(fc)<- paste0("eco:", core$bnum[core$padj <= 0.01])
pv <-core$padj[core$padj <= 0.01] # p-value
names(pv) <- paste0("eco:", core$bnum[core$padj <= 0.01])
kpg <- setNodeWeights(kpg, weights = alphaMLG(pv), defaultWeight = 1)
ref <- paste0("eco:",core$bnum)
peRes <- pe(x = fc, graphs = kpg, ref = ref,  nboot = 200, verbose = TRUE)

s <- Summary(peRes, pathNames = kpn,  totalAcc = FALSE, totalPert = FALSE,
             pAcc = FALSE, order.by = "pPert")



p <- peRes@pathways[["path:eco00650"]]
g <- layoutGraph(p@map, layoutType = "dot")
graphRenderInfo(g) <- list(fixedsize = FALSE)
edgeRenderInfo(g) <- peEdgeRenderInfo(p)
nodeRenderInfo(g) <- peNodeRenderInfo(p)
renderGraph(g)


#---------------------------------------------------------------------
library("topGO")
library("org.EcK12.eg.db")
results <- "/Users/annasintsova/git_repos/upec_mouse_model/results/"
# get the gene list and convert to entrez
goDat <- read.csv(file="/Users/annasintsova/git_repos/HUTI-RNAseq/private_results/differential_expression_analysis/gene_association.ecocyc.csv", row.names=1, stringsAsFactors=F) # 19361
myDat <- read.csv(file="/Users/annasintsova/git_repos/upec_mouse_model/results/HM43_LB_vs_HM43_UTI_all_residuals.csv", row.names=1, stringsAsFactors=F)
myDat <- read.csv(file="/Users/annasintsova/git_repos/upec_mouse_model/results/HM43_UR_vs_HM43_UTI_all_residuals.csv", row.names=1, stringsAsFactors=F)
myDat <- read.csv(file="/Users/annasintsova/git_repos/upec_mouse_model/results/HM43_mouse_vs_HM43_UTI_all_residuals.csv", row.names=1, stringsAsFactors=F)

#FCok
DE <- which( myDat$resid > 1.5*sd(myDat$resid))    # 212 - after removing those with no gene Symbol
#FDRok <- which(myDat$padj < 0.05)           
#DE <- intersect(FDRok,FCok) 

geneNam <- as.vector(myDat$Gene.Name) #3428
DEgeneNam <- as.vector(myDat$Gene.Name[DE]) #212

#geneNam1 <- unlist(strsplit(geneNam,","))  # 2653
#DEgeneNam1 <- unlist(strsplit(DEgeneNam,","))  

# Normalize the sets ----
geneNam2 <- geneNam[which(geneNam %in% goDat$Symbol)] # 2858
DEgeneNam2 <- DEgeneNam[which(DEgeneNam %in% goDat$Symbol)]  # 167
goDat2 <- goDat[which(goDat$Symbol %in% geneNam2),] # 17276,
EcoliGenes <- rep(0,length(geneNam2))
EcoliGenes[which(geneNam2 %in% DEgeneNam2  )] <- 1
names(EcoliGenes) <- geneNam2

DEcoli <- names(EcoliGenes)[which(EcoliGenes == 1)]
FacGenes <- as.factor(EcoliGenes)


# make GO to genes lists -----
goBP <- goDat2[which(goDat2$category == "P"),]
goMF <- goDat2[which(goDat2$category == "F"),]
goCC <- goDat2[which(goDat2$category == "C"),]

BPterms <- unique(goBP$GO)  # 1142
MFterms <- unique(goMF$GO)  # 1309
CCterms <- unique(goCC$GO)  # 124

BPlist <- list()
for(ii in 1:length(BPterms)) {
  term <- BPterms[ii]
  genes <- unique(goBP$Symbol[which(goBP$GO == term)])
  BPlist[[ii]] <-  genes
}
names(BPlist) <- BPterms

MFlist <- list()
for(ii in 1:length(MFterms)) {
  term <- MFterms[ii]
  genes <- unique(goMF$Symbol[which(goMF$GO == term)])
  MFlist[[ii]] <-  genes
}
names(MFlist) <- MFterms

CClist <- list()
for(ii in 1:length(CCterms)) {
  term <- CCterms[ii]
  genes <- unique(goCC$Symbol[which(goCC$GO == term)])
  CClist[[ii]] <-  genes
}
names(CClist) <- CCterms


topBP <- new("topGOdata", description = "Ecoli BP", ontology = "BP", allGenes = FacGenes,  nodeSize = 1, annot = annFUN.GO2genes, GO2genes = BPlist)
topMF <- new("topGOdata", description = "Ecoli MF", ontology = "MF", allGenes = FacGenes,  nodeSize = 1, annot = annFUN.GO2genes, GO2genes = MFlist)
topCC <- new("topGOdata", description = "Ecoli CC", ontology = "CC", allGenes = FacGenes,  nodeSize = 1, annot = annFUN.GO2genes, GO2genes = CClist)

BP.genes <- genesInTerm(topBP)
MF.genes <- genesInTerm(topMF)
CC.genes <- genesInTerm(topCC)


# BP -------------------
BP.Fisher.elim <- runTest(topBP, algorithm = "elim", statistic = "fisher")  
BP.Fisher.elim.Table <- GenTable(topBP, elimFisher=BP.Fisher.elim, topNodes=1000)
BP.Fisher.elim.Table <- BP.Fisher.elim.Table[which(BP.Fisher.elim.Table$elimFisher < 0.05),]
BP.Fisher.elim.Table <- BP.Fisher.elim.Table[, !(names(BP.Fisher.elim.Table) %in% c("Term"))]
if(nrow(BP.Fisher.elim.Table) > 0) { 
  goIDs <- BP.Fisher.elim.Table$GO.ID
  Term <- unlist(lapply(goIDs, function(x) (Term(x)[[1]])))
  genesInTerm <- sapply(goIDs, function(x) paste(BP.genes[[which(names(BP.genes) == x)]], collapse=","))
  DEGenesInTerm <- sapply(goIDs, function(x) paste(intersect(BP.genes[[which(names(BP.genes) == x)]], DEgeneNam), collapse=","))
  aDF <- data.frame(Term, genesInTerm,DEGenesInTerm)
  bDF <- cbind(BP.Fisher.elim.Table,aDF)
  write.csv(bDF, file=paste0(results, "down_mouse_outlier_pathway_analysis_BP.csv") )
}
# MF -------------------
MF.Fisher.elim <- runTest(topMF, algorithm = "elim", statistic = "fisher")  
MF.Fisher.elim.Table <- GenTable(topMF, elimFisher=MF.Fisher.elim, topNodes=1000)
MF.Fisher.elim.Table <- MF.Fisher.elim.Table[which(MF.Fisher.elim.Table$elimFisher < 0.05),]  # 46
MF.Fisher.elim.Table <- MF.Fisher.elim.Table[, !(names(MF.Fisher.elim.Table) %in% c("Term"))]

if(nrow(MF.Fisher.elim.Table) > 0) { 
  goIDs <- MF.Fisher.elim.Table$GO.ID
  Term <- unlist(lapply(goIDs, function(x) (Term(x)[[1]])))
  genesInTerm <- sapply(goIDs, function(x) paste(MF.genes[[which(names(MF.genes) == x)]], collapse=","))
  DEGenesInTerm <- sapply(goIDs, function(x) paste(intersect(MF.genes[[which(names(MF.genes) == x)]], DEgeneNam), collapse=","))
  aDF <- data.frame(Term, genesInTerm,DEGenesInTerm)
  bDF <- cbind(MF.Fisher.elim.Table,aDF)
  write.csv(bDF, paste0(results, "down_mouse_outlier_pathway_analysis_MF.csv") ) 
}

# CC -------------------

CC.Fisher.elim <- runTest(topCC, algorithm = "elim", statistic = "fisher")  
CC.Fisher.elim.Table <- GenTable(topCC, elimFisher=CC.Fisher.elim, topNodes=200)
CC.Fisher.elim.Table <- CC.Fisher.elim.Table[which(CC.Fisher.elim.Table$elimFisher < 0.05),]  # 2
CC.Fisher.elim.Table <- CC.Fisher.elim.Table[, !(names(CC.Fisher.elim.Table) %in% c("Term"))]

if(nrow(CC.Fisher.elim.Table) > 0) { 
  goIDs <- CC.Fisher.elim.Table$GO.ID
  Term <- unlist(lapply(goIDs, function(x) (Term(x)[[1]])))
  genesInTerm <- sapply(goIDs, function(x) paste(CC.genes[[which(names(CC.genes) == x)]], collapse=","))
  DEGenesInTerm <- sapply(goIDs, function(x) paste(intersect(CC.genes[[which(names(CC.genes) == x)]], DEgeneNam), collapse=","))
  aDF <- data.frame(Term, genesInTerm,DEGenesInTerm)
  bDF <- cbind(CC.Fisher.elim.Table,aDF)
  write.csv(bDF, paste0(results, "down_mouse_outlier_pathway_analysis_CC.csv")) 
}
# -----


