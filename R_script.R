#dir.create("/run/media/senthilkumar/New/GSE163877_bulkrna")
setwd("/run/media/senthilkumar/New/GSE163877_bulkrna")
#getwd()
packages <- c("GEOquery", "DESeq2", "EnhancedVolcano", "data.table", "edgeR", "limma", "ggplot2", "ggrepel", "ggfortify", "stats", "sva", "magrittr", "dplyr", "tidyverse")
lapply(packages, library, character.only = TRUE)
#gse <- getGEO("GSE163877", GSEMatrix = TRUE)
#exp_mat <- getGEOSuppFiles("GSE163877", fetch_files = T)
file = "/run/media/senthilkumar/New/GSE163877_bulkrna/GSE163877/GSE163877_VBB_Counts.txt.gz"
#gse <- getGEO("GSE163877", GSEMatrix = TRUE)
list.files()
count_matrix = read.csv(file = file, sep="\t")
rownames(count_matrix) <- count_matrix[[1]]
count_matrix <- count_matrix[,-1]
keep = rowSums2(count_matrix > 0) >= 3
filtered_matirx <- count_matrix[keep, ]
meta <- data.frame(sample = metadata$geo_accession, condition = metadata$characteristics_ch1.1)
rownames(meta) <- meta$sample
meta <- meta %>% mutate( condition = gsub("diagnosis: Alzheimer's disease", "AD", condition)) %>% mutate(condition = gsub("diagnosis: None", "None", condition))
colnames(count_matrix)
filtered_matirx %>% rename("GSM4989044" = "Patient1") %>% 
rename("GSM4989045" = "Patient2") %>% 
rename("GSM4989046" = "Patient3") %>%
rename("GSM4989047" = "Patient4") %>%
rename("GSM4989048" = "Patient5") %>%
rename("GSM4989049" = "Patient6") %>%
rename("GSM4989050" = "Patient7")
dge <- DGEList(counts = filtered_matirx, samples = meta)
dge <- calcNormFactors(dge, method = "TMM")
dge_v <- voom(dge, plot = T)
#saveRDS(dge_v, file = "dge_v.rds")
dge_v$targets[["condition"]]

comparison <- "AD-None"
design <- model.matrix(~0+dge_v$targets[['condition']])
colnames(design)
colnames(design) <- gsub(".*]]", "", colnames(design))
contrasts_matrix <- makeContrasts(contrasts = comparison, levels = design)
#remove(contrasts)
fit <- lmFit(dge_v, design = design)
fit2 <- contrasts.fit(fit = fit, contrasts = contrasts_matrix)
fit2 <- eBayes(fit2)
dge_list <- topTable(fit = fit2, p.value = 1, colnames(contrasts_matrix), n=Inf, sort.by = 'p')
#topTable(fit = fit2, p.value = 0.05, colnames(contrasts_matrix), n=Inf, sort.by = 'p')
#dge_list
#fwrite(dge_list, "DEG_P1_lfc0.tsv", sep = "\t", row.names = T) 
#sum(dge_list$adj.P.Val < 0.05)
#cat("Genes with raw p-value < 0.05: ", sum(dge_list$P.Value < 0.05), "\n")
#sum(dge_list$P.Value < 0.05)
library(EnhancedVolcano)
library(ggrepel)
EnhancedVolcano(dge_list,
                lab = rownames(dge_list),
                x = 'logFC',
                y = 'P.Value',                 # Uses raw P.Value since adj.P.Val has no signal
                pCutoff = 0.05,                # Custom raw p-value threshold line
                FCcutoff = 1,               # Custom log2 Fold Change threshold line (1.5-fold change)
                pointSize = 1.0,
                labSize = 4.0,
                title = 'Volcano Plot: AD vs None (GSE163877)',
                subtitle = 'Filtered by Raw P-Value < 0.05 and |log2FC| > 1.0',
                legendPosition = 'right',
                legendLabSize = 10,
                legendIconSize = 3.0,
                drawConnectors = TRUE,         # Draws lines connecting gene labels to their points
                widthConnectors = 0.5,
                colConnectors = 'black')

library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
#library(org.Mm.eg.db)
# library(org.Rn.eg.db)
library(ggplot2)
library(DOSE)
library(enrichplot)
p_thershold <- 0.05
fc_threshold <- 1.0


# enrichment --------------------------------------------------------------

dge_ordered <- dge_list[order(dge_list$logFC, decreasing = T), ]
head(dge_ordered)
logFC <- dge_ordered$logFC
r_names <- rownames(dge_ordered)
names(logFC) <- r_names

gene_enrichment_go <- gseGO(geneList = logFC, OrgDb = org.Hs.eg.db, pvalueCutoff = 0.05, ont = "ALL", verbose = T,keyType = "ENSEMBL")
gene_enrichment_go_df <- gene_enrichment_go@result
fwrite(gene_enrichment_go_df, "gene_enrichment_go_df.tsv", sep = "\t")
gsea_dotplot <- dotplot(gene_enrichment_go, showCategory = 20, orderBy="GeneRatio", label_format = 50)
?dotplot
ggsave("dotplot_enrich_go_gsea.png", gsea_dotplot, device = "png", units = "cm", width = 26, height = 18)
