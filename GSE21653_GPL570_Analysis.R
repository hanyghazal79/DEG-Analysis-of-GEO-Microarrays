#===| INITIALIZING ANALYZER MODEL REFERENCE |===
source("C:/Users/Dr.Hany/RstudioProjects/GenesDiffExpAnalysis/AnalyzerModel.R")
analyzer<-AnalyzerModel(geo = paste(getwd(),"/GSE21653_series_matrix.txt.gz", sep = ""),
                        id="GSE21653_GPL570")
#===| GETTING GEO DATA |===
data<-analyzer$getGEOdata()
gse<-data[["gse"]]
gex<-data[["gex"]]
anno<-data[["anno"]]
annogex<-data[["annogex"]]
sampleInfo<-data[["sampleInfo"]]
#===| SELECTION OF DESIGN FEATURES FROM SAMPLES |===
sampleInfo<-dplyr::select(sampleInfo, group1 = characteristics_ch1.10, group2 = characteristics_ch1.11, group3 = characteristics_ch1.12)
filter<-sampleInfo%>%filter(group3 == "")
filter$group3<-filter$group1
sampleInfo[rownames(filter),]$group3<-filter$group3
sampleInfo<-dplyr::select(sampleInfo, group = group3)
sampleInfo<-sampleInfo%>%mutate(group=gsub("molecular subtype: ", "", group),
                                group=gsub(" ", "", group))
unique(sampleInfo$group)
sampleInfo$group<-sapply(sampleInfo$group, function(x){
  if(x %in% c("ERBB2", "LuminalA", "LuminalB", "tissue:breastcancertumor")) x <- "NONTNBC"
  else if(x == "Basal") x <- "TNBC"
  else x <- "HEALTHY"
})
#===| CREATING A DESIGN MATRIX AND CONTRASTS |===
designMatrix <- model.matrix(~0+sampleInfo$group)
colnames(designMatrix) <- c("HEALTHY", "NONTNBC", "TNBC")
contrasts <- makeContrasts(TNBC-NONTNBC, TNBC-HEALTHY, NONTNBC-HEALTHY, levels = designMatrix)
#===| FITTING HIGHLY EXPRESSED GENES ABOVE MEDIAN |===
fitResult <- analyzer$fitHighExpGenes(gex, annogex, designMatrix, contrasts)
keptgex <- fitResult[["keptgex"]]
fit <- fitResult[["fit"]]
#===| DECIDE TESTING |===
decide <- analyzer$makeDecideTest(fit)
decideResult <- decide[["result"]]
decideDetails <- decide[["details"]]
##===| DIFFERENTIAL EXPRESSION RESULT FOR [TNBC-NONTNBC] CONTRAST|===
deg_TNBC_NONTNBC <- analyzer$getDiffExpResult(fit, coef = 1, annogex = annogex, "TNBC_NONTNBC")
deg_ranked <- analyzer$rankGenes(deg_TNBC_NONTNBC, coef = 1, contrast = "TNBC_NONTNBC")
deg_top_ranked <- analyzer$getTopRankedGenes(deg_ranked, d = 0.95, coef = 1, contrast = "TNBC_NONTNBC")
topGeneMatrix <- deg_top_ranked[["geneMatrix"]]
topGeneNames <- deg_top_ranked[["geneNames"]]
#===|<HEATMAP TOP RANKED GENE MATRIX>|===
analyzer$heatMapTopGeneMatrix(topGeneMatrix, topGeneNames, 1, "TNBC_NONTNBC")
anno_top_genes <- analyzer$annotateGenes(topGeneNames, "TNBC_NONTNBC")

#=======|<ANALYZING CONTRAST 2 [TNBC-HEALTHY]>|=======
#===| DIFFERENTIAL EXPRESSION RESULT FOR [TNBC-HEALTHY] CONTRAST|===
deg_TNBC_HEALTHY <- analyzer$getDiffExpResult(fit, coef = 2, annogex = annogex, "TNBC-HEALTHY")
deg_ranked_2 <- analyzer$rankGenes(deg_TNBC_HEALTHY, coef = 2, contrast = "TNBC-HEALTHY")
deg_top_ranked_2 <- analyzer$getTopRankedGenes(deg_ranked_2, d = 0.95, coef = 2, contrast = "TNBC-HEALTHY")
topGeneMatrix_2 <- deg_top_ranked_2[["geneMatrix"]]
topGeneNames_2 <- deg_top_ranked_2[["geneNames"]]
#===|<HEATMAP TOP RANKED GENE MATRIX>|===
analyzer$heatMapTopGeneMatrix(topGeneMatrix_2, topGeneNames_2, 2, "TNBC-HEALTHY")
anno_top_genes_2 <- analyzer$annotateGenes(topGeneNames_2, "TNBC-HEALTHY")
##===| DIFFERENTIAL EXPRESSION RESULT FOR [NONTNBC-HEALTHY] CONTRAST|===
deg_NONTNBC_HEALTHY <- analyzer$getDiffExpResult(fit, coef = 3, annogex = annogex, "NONTNBC-HEALTHY")
deg_ranked_3 <- analyzer$rankGenes(deg_NONTNBC_HEALTHY, coef = 3, contrast = "NONTNBC-HEALTHY")
deg_top_ranked_3 <- analyzer$getTopRankedGenes(deg_ranked_3, d = 0.95, coef = 3, contrast = "NONTNBC-HEALTHY")
topGeneMatrix_3 <- deg_top_ranked_3[["geneMatrix"]]
topGeneNames_3 <- deg_top_ranked_3[["geneNames"]]
#===|<HEATMAP TOP RANKED GENE MATRIX>|===
analyzer$heatMapTopGeneMatrix(topGeneMatrix_3, topGeneNames_3, 3, "NONTNBC-HEALTHY")
anno_top_genes <- analyzer$annotateGenes(topGeneNames_3, "NONTNBC-HEALTHY")
