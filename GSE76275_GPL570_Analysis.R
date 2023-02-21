#===| INITIALIZING ANALYZER MODEL REFERENCE |===
source("C:/Users/Dr.Hany/RstudioProjects/GenesDiffExpAnalysis/AnalyzerModel.R")
analyzer<-AnalyzerModel(geo = paste(getwd(),"/GSE76275_series_matrix.txt.gz", sep = ""),
                        id="GSE76275_GPL570")
#===| GETTING GEO DATA |===
data<-analyzer$getGEOdata()
gse<-data[["gse"]]
gex<-data[["gex"]]
anno<-data[["anno"]]
annogex<-data[["annogex"]]
sampleInfo<-data[["sampleInfo"]]
# SUBSERIES GSE76124 = TNBC 
analyzerV2<-AnalyzerModel(geo = paste(getwd(),"/GSE76124_series_matrix.txt.gz", sep = ""),id="GSE76124_GPL570")
dataGSE76124<- analyzerV2$getGEOdata()
samplesTNBC<-dataGSE76124[["sampleInfo"]]
rownames(samplesTNBC)
# SUBSERIES GSE76274 = NONTNBC
analyzerV3<-AnalyzerModel(geo = paste(getwd(),"/GSE76274_series_matrix.txt.gz", sep = ""),id="GSE76274_GPL570")
dataGSE76274<- analyzerV2$getGEOdata()
samplesNONTNBC<-dataGSE76274[["sampleInfo"]]
# CONFIG SAMPLE INFO OF THE WHOLE DATASET

#===| SELECTION OF DESIGN FEATURES FROM SAMPLES |===
sampleInfo<-sampleInfo%>%mutate(group=rep("", n()))
sampleInfo<-dplyr::select(sampleInfo, group=group)
sampleInfo$group<-sapply(rownames(sampleInfo), function(x){
  if(x %in% rownames(samplesTNBC)) sampleInfo$group<-"TNBC"
  else sampleInfo$group<-"NONTNBC"
})
#===| CREATING A DESIGN MATRIX AND CONTRASTS |===
designMatrix <- model.matrix(~0+sampleInfo$group)
colnames(designMatrix) <- c("NONTNBC", "TNBC")
contrasts <- makeContrasts(TNBC-NONTNBC, levels = designMatrix)
#===| FITTING HIGHLY EXPRESSED GENES ABOVE MEDIAN |===
fitResult <- analyzer$fitHighExpGenes(gex, annogex, designMatrix, contrasts)
keptgex <- fitResult[["keptgex"]]
fit <- fitResult[["fit"]]
#===| DECIDE TESTING |===
decide <- analyzer$makeDecideTest(fit)
decideResult <- decide[["result"]]
decideDetails <- decide[["details"]]
#===| DIFFERENTIAL EXPRESSION RESULT FOR [TNBC-NONTNBC] CONTRAST|===
deg_TNBC_NONTNBC <- analyzer$getDiffExpResult(fit, coef = 1, annogex = annogex, "TNBC_NONTNBC")
deg_ranked <- analyzer$rankGenes(deg_TNBC_NONTNBC, coef = 1, contrast = "TNBC_NONTNBC")
deg_top_ranked <- analyzer$getTopRankedGenes(deg_ranked, d = 0.95, coef = 1, contrast = "TNBC_NONTNBC")
topGeneMatrix <- deg_top_ranked[["geneMatrix"]]
topGeneNames <- deg_top_ranked[["geneNames"]]
#===|<HEATMAP TOP RANKED GENE MATRIX>|===
analyzer$heatMapTopGeneMatrix(topGeneMatrix, topGeneNames, 1, "TNBC_NONTNBC")
anno_top_genes <- analyzer$annotateGenes(topGeneNames, "TNBC_NONTNBC")
