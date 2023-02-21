#===| INITIALIZING ANALYZER MODEL REFERENCE |===
source("C:/Users/Dr.Hany/RstudioProjects/GenesDiffExpAnalysis/AnalyzerModel.R")
analyzer<-AnalyzerModel(geo = paste(getwd(),"/GSE36693_series_matrix.txt.gz", sep = ""),id="GSE36693_GPL10558")
#===| GETTING GEO DATA |===
data<-analyzer$getGEOdata()
gse<-data[["gse"]]
gex<-data[["gex"]]
anno<-data[["anno"]]
annogex<-data[["annogex"]]
sampleInfo<-data[["sampleInfo"]]
#===| SELECTION OF DESIGN FEATURES FROM SAMPLES |===
sampleInfo<-dplyr::select(sampleInfo, title)
sampleInfo$title<-sapply(sampleInfo$title, function(x){
  x <- str_split(x, " ")[[1]][[2]]
})
unique(sampleInfo$title)
sampleInfo$group<-sampleInfo$title
sampleInfo<-dplyr::select(sampleInfo, group)
sampleInfo$group<-sapply(sampleInfo$group, function(x){
  if(x == "Basal") x <- "TNBC"
  else if(x == "Basal-like") x <- "BASALLIKE"
  else x <- "NONTNBC"
})
#===| CREATING A DESIGN MATRIX AND CONTRASTS |===
designMatrix <- model.matrix(~0+sampleInfo$group)
colnames(designMatrix) <- c("BASALLIKE", "NONTNBC", "TNBC")
contrasts <- makeContrasts(TNBC-NONTNBC, TNBC-BASALLIKE, BASALLIKE-NONTNBC, levels = designMatrix)
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
deg_top_ranked <- analyzer$getTopRankedGenes(deg_ranked, d = 0.985, coef = 1, contrast = "TNBC_NONTNBC")
topGeneMatrix <- deg_top_ranked[["geneMatrix"]]
topGeneNames <- deg_top_ranked[["geneNames"]]
#===|<HEATMAP TOP RANKED GENE MATRIX>|===
analyzer$heatMapTopGeneMatrix(topGeneMatrix, topGeneNames, 1, "TNBC_NONTNBC")
anno_top_genes <- analyzer$annotateGenes(topGeneNames, "TNBC_NONTNBC")

#=======|<ANALYZING CONTRAST 2 [TNBC-BASALLIKE]>|=======
#===| DIFFERENTIAL EXPRESSION RESULT FOR [TNBC-BASALLIKE] CONTRAST|===
deg_TNBC_BASALLIKE <- analyzer$getDiffExpResult(fit, coef = 2, annogex = annogex, "TNBC-BASALLIKE")
deg_ranked_2 <- analyzer$rankGenes(deg_TNBC_BASALLIKE, coef = 2, contrast = "TNBC-BASALLIKE")
deg_top_ranked_2 <- analyzer$getTopRankedGenes(deg_ranked_2, d = 0.985, coef = 2, contrast = "TNBC-BASALLIKE")
topGeneMatrix_2 <- deg_top_ranked_2[["geneMatrix"]]
topGeneNames_2 <- deg_top_ranked_2[["geneNames"]]
#===|<HEATMAP TOP RANKED GENE MATRIX>|===
analyzer$heatMapTopGeneMatrix(topGeneMatrix_2, topGeneNames_2, 2, "TNBC-BASALLIKE")
anno_top_genes_2 <- analyzer$annotateGenes(topGeneNames_2, "TNBC-BASALLIKE")
##===| DIFFERENTIAL EXPRESSION RESULT FOR [BASALLIKE-NONTNBC] CONTRAST|===
deg_BASALLIKE_NONTNBC <- analyzer$getDiffExpResult(fit, coef = 3, annogex = annogex, "BASALLIKE-NONTNBC")
deg_ranked_3 <- analyzer$rankGenes(deg_BASALLIKE_NONTNBC, coef = 3, contrast = "BASALLIKE-NONTNBC")
deg_top_ranked_3 <- analyzer$getTopRankedGenes(deg_ranked_3, d = 0.985, coef = 3, contrast = "BASALLIKE-NONTNBC")
topGeneMatrix_3 <- deg_top_ranked_3[["geneMatrix"]]
topGeneNames_3 <- deg_top_ranked_3[["geneNames"]]
#===|<HEATMAP TOP RANKED GENE MATRIX>|===
analyzer$heatMapTopGeneMatrix(topGeneMatrix_3, topGeneNames_3, 3, "BASALLIKE-NONTNBC")
anno_top_genes <- analyzer$annotateGenes(topGeneNames_3, "BASALLIKE-NONTNBC")
