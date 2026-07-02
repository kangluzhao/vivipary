#########################################V.PhyloMaker2#########################################
library("V.PhyloMaker2")
exam <- read.table("species.txt", sep = "\t", header = T)
tree <- phylo.maker(sp.list = exam, tree = GBOTB.extended.TPL, nodes = nodes.info.1.TPL, scenarios = "S3")
write.tree(tree$scenario.3, "species.tree")

#########################################TissueEnrich#########################################
library(TissueEnrich)
expr_stage_mean <- read.csv("OGsumExpr.csv", check.names = F, header = T, row.names = 1)
expr_stage_mean <- as.matrix(expr_stage_mean)
Summarized_data <- SummarizedExperiment(assays = SimpleList(expr_stage_mean),
                           rowData = row.names(expr_stage_mean),
                           colData = colnames(expr_stage_mean))

GeneRetrieval_output <- teGeneRetrieval(Summarized_data,
                                        foldChangeThreshold = 5,
                                        maxNumberOfTissues = 7,
                                        expressedGeneThreshold = 1)

result <- as.data.frame(assay(GeneRetrieval_output))
write.table(result, "TissueEnrichResult.txt", quote = F, sep = '\t', row.names = F)

dryseed <- result %>% filter(Tissue == "dry seeds" & Group == "Tissue-Enhanced")
expr_stage_mean <- as.data.frame(expr_stage_mean) %>% rownames_to_column("Gene")
dryseedFPKM <- left_join(dryseed, expr_stage_mean, by = "Gene")
write.csv(dryseedFPKM, "dryseedFPKM.csv", row.names = F)

#########################################circlize#########################################
library(circlize)
library(ComplexHeatmap)

graphics.off() 
circos.clear()

df <- read.csv("dryseedFPKMheatmap.csv", row.names = 1, check.names = FALSE)

normalize_01 <- function(x) {
  rng <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
  if (rng == 0) return(rep(0, length(x)))
  return((x - min(x, na.rm = TRUE)) / rng)
}

df_norm <- apply(df, 2, normalize_01)
rownames(df_norm) <- rownames(df)

hc_genes <- hclust(dist(t(df_norm)), method = "ward.D2")
df_ready <- df_norm[, hc_genes$order]

pdf("Circular_Heatmap.pdf", width = 14, height = 14)

circos.par(
  start.degree = 0,            
  gap.degree = 2,              
  canvas.xlim = c(-1.6, 1.6),  
  canvas.ylim = c(-1.6, 1.6)  
)

col_fun <- colorRamp2(c(0, 0.5, 1), c("#2166ac", "white", "#b2182b"))

circos.heatmap(
  df_ready, 
  col = col_fun,
  cluster = TRUE,              
  dend.side = "inside",        
  dend.track.height = 0.2,    
  rownames.side = "outside", 
  rownames.cex = 0.6, 
  rownames.font = 1,  
  track.height = 0.2
)

lgd <- Legend(
  title = "Relative\nexpression", 
  col_fun = col_fun, 
  direction = "vertical"
)
draw(lgd, x = unit(0.9, "npc"), y = unit(0.15, "npc"))

dev.off()

