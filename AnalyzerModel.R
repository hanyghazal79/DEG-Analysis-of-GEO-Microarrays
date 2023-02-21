library(GEOquery)
library(dplyr)
library(tidyverse)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(FactoMineR)
library(factoextra)
library(bigutilsr)
library(readr)
library(stringr)
library(Biobase)
library(limma)
library(org.Hs.eg.db)
library(survminer)
library(survival)
library(methods)
library(annotate)
library(pd.hugene.1.0.st.v1)
library(desirability)
library(desiR)
AnalyzerModel<-setRefClass("AnalyzerModel",
                           fields = list(geo="character", id="character"),
                           methods = list(
                             getGEOdata = function(){
                               gse<-NULL
                               if (grepl("/", geo) == TRUE)
                                 gse<-getGEO(filename = geo, GSEMatrix = TRUE)
                               else
                                 gse<-getGEO(geo, destdir = ".")
                               
                               if(length(gse)>1)
                                 gse<-gse[[1]]
                               #==
                               sampleInfo <- pData(gse)
                               gex<-exprs(gse)
                               
                               anno<-fData(gse)
                               if(grepl("GPL570", id)){
                                 colnames(anno)[11]<-"GENE_SYMBOL"
                                 #anno<- dplyr::select(anno, ID, GENE_SYMBOL,ENTREZ_GENE_ID)
                               }
                               else if(grepl("GPL14877", id)){
                                 keysVector<-sapply(anno$SPOT_ID, function(x){
                                   x <- as.character(x)
                                 })
                                 symbolData<-AnnotationDbi::select(org.Hs.eg.db,
                                                                   columns = c("SYMBOL"),
                                                                   keys = keysVector,
                                                                   keytype = "ENTREZID")
                                 anno$GENE_SYMBOL<-symbolData$SYMBOL
                                 anno$ENTREZ_GENE_ID<-symbolData$ENTREZID
                                 
                                 #anno<- dplyr::select(anno, ID, GENE_SYMBOL,ENTREZ_GENE_ID)
                                 
                               }
                               else if(grepl("GPL10558", id)){
                                 colnames(anno)[10] <- "ENTREZ_GENE_ID"
                                 colnames(anno)[13] <- "GENE_SYMBOL"
                                 #anno <- dplyr::select(anno, ID, GENE_SYMBOL,ENTREZ_GENE_ID)
                                 
                               }
                               else if(grepl("GPL13534", id)){
                                 anno$GENE_SYMBOL <- anno$UCSC_RefGene_Name
                                 anno <- anno %>% mutate(ENTREZ_GENE_ID = rep("", n()))
                               }
                               else if(grepl("GPL6244", id)){
                                 anno$GENE_SYMBOL <- anno$gene_assignment
                                 anno <- anno %>% mutate(ENTREZ_GENE_ID = rep("", n()))
                                 
                               }
                               anno<- dplyr::select(anno, ID, GENE_SYMBOL,ENTREZ_GENE_ID)
                               annogex <- cbind(anno, gex)
                               
                               result <- list(gse = gse, 
                                              sampleInfo = sampleInfo, 
                                              anno = anno, 
                                              gex = gex, 
                                              annogex = annogex)
                               return(result)
                             },
                             getPCNAcor = function(annogex){
                               PCNA.cor <- NULL
                               if("PCNA" %in% annogex$GENE_SYMBOL){
                                 loc <- filter(annogex, GENE_SYMBOL == "PCNA")
                                 probe_id <- loc$ID
                                 df<-annogex[ , !colnames(annogex) %in% c("ID", "GENE_SYMBOL", "ENTREZ_GENE_ID") ]
                                 df <- df %>% mutate_if(is.character, as.double)
                                 df <- t(df)
                                 PCNA.cor <- cor(df, df[ , probe_id])
                                 PCNA.cor <- as.data.frame(PCNA.cor)
                                 colnames(PCNA.cor) <- "PCNA.cor"
                               }
                               
                             },
                             getPCA = function(gex){
                               pca <- prcomp(t(gex))
                               return(pca)
                             },
                             # 
                             # cleanSampleColumns = function(sampleInfo){
                             #   
                             # },
                             # makeContrastsFromDesign = function(sampleInfo){
                             #   
                             # },
                             fitHighExpGenes = function(gex, annogex, designMatrix, contrasts){
                               cutoff<-median(gex)
                               highGEX<- gex > cutoff
                               kept<- rowSums(highGEX)>2
                               keptgex<-gex[kept,]
                               aw<-arrayWeights(keptgex, designMatrix)
                               limmFit<-lmFit(keptgex, designMatrix, weights = aw)
                               contrastFit<-contrasts.fit(limmFit, contrasts)
                               eBayesfit<-eBayes(contrastFit)
                               eBayesfit$genes<-annogex[kept,]
                               result<-list(keptgex = keptgex, fit = eBayesfit)
                               return(result)
                             },
                             makeDecideTest = function(eBayesfit){
                               decideTestDF<-as.data.frame(table(decideTests(eBayesfit)))
                               decideTestSummary<-as.data.frame(summary(decideTests(eBayesfit)))
                               result<- list(result=decideTestDF, details=decideTestSummary)
                               return(result)
                             },
                             #getDiffExpResult = function(fit, coef = "integer", p_cutoff="numeric", fc_cutoff="numeric", topN="integer", keptgex, annogex){
                            getDiffExpResult = function(fit, coef = "integer", annogex, noteStr){
                                 
                               df<-annogex[ , !colnames(annogex) %in% c("ID", "GENE_SYMBOL", "ENTREZ_GENE_ID") ]
                               df <- df %>% mutate_if(is.character, as.double)
                               #
                               sd <- apply(df, 1, sd, na.rm = TRUE)
                               sd <- as.data.frame(sd)
                               #
                               df <- t(df)
                               PCNA.cor <- NULL
                               if("PCNA" %in% annogex$GENE_SYMBOL){
                                 loc <- filter(annogex, GENE_SYMBOL == "PCNA")
                                 probe_id <- loc$ID
                                 PCNA.cor <- cor(df, df[ , probe_id])
                                 PCNA.cor <- as.data.frame(PCNA.cor)
                                 colnames(PCNA.cor) <- "PCNA.cor"
                               }
                               # GETTING DEGs
                               deg <- topTable(fit, coef, number = Inf)
                               cols <- c("ID", 
                                         "GENE_SYMBOL", 
                                         "ENTREZ_GENE_ID", 
                                         "logFC", 
                                         "AveExpr", 
                                         "t",
                                         "P.Value",
                                         "adj.P.Val",
                                         "B")
                               deg_short <- deg[ , cols]
                               sd <- sd[rownames(deg_short), c("sd")]
                               PCNA.cor <- PCNA.cor[rownames(deg_short), c("PCNA.cor")]
                               deg_short$SD <- sd
                               deg_short$PCNA.cor <- PCNA.cor
                               write.csv(deg_short, file = paste(id, "_DEG_", coef, "_", noteStr, ".csv", sep = ""))
                               return(deg_short)
                             },
                            rankGenes = function(deg, coef = "numeric", contrast = "character"){
                              # > DESIRABILITY RANKING PLOTTING OF DESI LINES [desiR package]<
                              xcut1 <- median(deg[ , c("AveExpr")])
                              xcut2 <- (xcut1 + max(deg[ , c("AveExpr")]))/2
                              sdcut1 <- median(deg[ , c("SD")])
                              sdcut2 <- (sdcut1 + max(deg[ , c("SD")]))/2
                              fc_cut1 <- min(deg[ , c("logFC")]) / 2
                              fc_cut2 <- fc_cut1 / 2
                              fc_cut4 <- max(deg[ , c("logFC")]) / 2
                              fc_cut3 <- fc_cut4 / 2
                              
                              png(filename = paste(analyzer$id, "_desi_lines_", coef, "_", contrast, ".png", sep = ""), width = 1400, height = 900)
                              par(mfrow=c(2,3), mar=c(4, 4, 4, 4), las=1)
                              # CUTTING AveExpr
                              hist(abs(deg[ , c("AveExpr")]), breaks = 30, col = "red", border = "white", main = "", xlab ="Mean expression")
                              des.line(abs(deg[ , c("AveExpr")]), "d.high", des.args=c(cut1 = xcut1, cut2 = xcut2, scale = 0.5))
                              # CUTTING SD FOR GENE EXP
                              hist(deg[ , c("SD")], breaks = 30, col = "red", border = "white", main = "", xlab = "Within group SD")
                              des.line(deg[ , c("SD")], "d.high", des.args = c(cut1 = sdcut1, cut2 = sdcut2, scale = 0.5))
                              # CUTTING P.value
                              hist(deg[ , c("P.Value")], breaks = 50, col = "red", border = "white", main = "", xlab = "p-value")
                              des.line(deg[ , c("P.Value")], "d.low", des.args = c(cut1 = 0.01, cut2 = 0.05, scale = 0.5))
                              # CUTTING logFC
                              hist(deg[ , c("logFC")], breaks = 30, col = "red", border = "white", main = "", xlab = "Log FC")
                              des.line(deg[ , c("logFC")], "d.ends", des.args = c(cut1 = fc_cut1, cut2 = fc_cut2, cut3 = fc_cut3, cut4 = fc_cut4, scale = 0.5))
                              # CUTTING PCNA.cor
                              hist(abs(deg[ , c("PCNA.cor")]), breaks = 50, col = "red", border = "white", main = "", xlab = "Absolute correlation with PCNA")
                              des.line(abs(deg[ , c("PCNA.cor")]), "d.low", des.args = c(cut1 = 0.15, cut2 = 0.3, scale = 0.5))
                              dev.off()
                              # <SETTING DESIRABILITY FEATURES>
                              deg$d1 <- d.high(abs(deg[ , c("AveExpr")]), cut1 = xcut1, cut2 = xcut2, scale = 0.5)
                              deg$d2 <- d.high(deg[ , c("SD")], cut1 = sdcut1, cut2 = sdcut2, scale = 0.5)
                              deg$d3 <- d.low(deg[ , c("P.Value")], cut1 = 0.01, cut2 = 0.05, scale = 0.5)
                              deg$d4 <- d.ends(deg[ , c("logFC")], cut1 = fc_cut1, cut2 = fc_cut2, cut3 = fc_cut3, cut4 = fc_cut4, scale = 0.5)
                              deg$d5 <- d.low(abs(deg[ , c("PCNA.cor")]), cut1 = 0.15, cut2 = 0.3, scale = 0.5)
                              deg$D <- d.overall(deg$d1, deg$d2, deg$d3, deg$d4, deg$d5, weights = c(0.1, 0.1, 1, 1, 0.5))
                              # <PLOTTING OVERALL DESIRABILITY>
                              png(filename = paste(analyzer$id, "_overall_desirability_", coef, "_", contrast, ".png", sep = ""), width = 900, height = 500)
                              plot(rev(sort(deg$D)), type = "l", xlab = "Rank", ylab = "Overall Desirability")
                              dev.off()
                              # <>
                              deg_ranked <- arrange(deg, desc(D))
                              write.csv(deg_ranked, file = paste(id, "_RankedDEG_", coef, "_", contrast, ".csv", sep = ""))
                              return(deg_ranked)
                              
                            },
                            getTopRankedGenes = function(deg_ranked, d = "numeric", coef = "numeric", contrast = "character"){
                              
                              #png(filename = paste(analyzer$id, "_VolcanoPlot_DEG_",coef, "_", contrast, ".png", sep = ""), width = 1440, height = 960)
                              ggplot(deg_ranked, aes(x = logFC, y = B, col = D, label = ifelse(D >= d, GENE_SYMBOL, ""))) +
                                geom_point() +
                                geom_text_repel()
                              #dev.off()
                              ggsave(paste(analyzer$id, "_VolcanoPlot_DEG_",coef, "_", contrast, ".png", sep = ""), width = 15.9, height = 10)
                              top_DEGs <- filter(deg_ranked, D >= d)
                              geneMatrix = gex[top_DEGs$ID, ]
                              geneNames <- annogex[top_DEGs$ID, "GENE_SYMBOL"]
                              result <- list(top_DEGs = top_DEGs, geneMatrix = geneMatrix, geneNames = geneNames)
                              write.csv(top_DEGs, file = paste(id, "_topDEG_", coef, "_", contrast, ".csv", sep = ""))
                              
                              return(result)
                            }
                            ,
                             heatMapTopGeneMatrix = function(geneMatrix, geneNames, coef, contrast){
                               png(filename = paste(id,"_heatmap_Top_DEG_",coef,"_", contrast, ".png", sep = ""), width = 1440, height = 960)
                               pheatmap(geneMatrix, labels_row = geneNames, scale = "row")
                               dev.off()
                             },
                             annotateGenes = function(geneNames, noteStr){
                               anno <- AnnotationDbi::select(org.Hs.eg.db,
                                                             columns = c("ENSEMBL", "GO"),
                                                             keys = geneNames,
                                                             keytype = "SYMBOL")
                               write.csv(anno, file = paste(id, "_ENSEMBL_GO_anno_", noteStr, ".csv", sep = ""))
                               
                               return(anno)
                             }
                             
                           )
)





