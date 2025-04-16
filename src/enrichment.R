# NO LONGER USING BUT KEEPING IT FOR REFERENCE JUST IN CASE

#install.packages("BiocManager")
#BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "DOSE"))

library(clusterProfiler)
library(org.Hs.eg.db)    # human gene annotation
library(enrichplot)      # for plots
library(dplyr)          # needed for some plotting


# pathway enrichment analysis: see what pathways the significant genes are involved in and if any are overrepresented (suggests that pathway may be dysregulated for that cancer)
    # use KEGG (curated pathways) and GO (more broad/general annotations) databases 
    # if run only on mutexc, if results show: the pairs are enriched for the same pathway (suggests redundancy), enriched for DNA repair aka occur in DNA repair functions/pathways more often (potential synthetic lethality insights), shared across cancers (maybe key vulnerability)
    # if run only on co: if two genes are in complimentary pathways (maybe targets for combo therapy)
# maybe can't do it on co because only 3 cohorts have significant results and even for those they have very few pairs of genes
# going to do per-cancer pathway enrichment for mutually exclusive pairs for all cancers that have at least 5 significant pairs per cohort
# find literature sources later for justification


# function to run enrichment for a singular cancer type / cohort
enrich_mutexc <- function(cancer_type, df) {
  genes <- unique(toupper(c(c(df$gene1, df$gene2))))  # get unique gene names
  
  if (length(genes) < 10) {  # threshold since don't want to run analysis if too few genes
    message(paste("Skipping", cancer_type, "- fewer than 10 genes"))
    return(NULL)
  }
  
  gene_df <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)  # convert to entrez ids
  entrez_ids <- unique(gene_df$ENTREZID)

  kegg_res <- enrichKEGG(gene = entrez_ids, organism = "hsa",
                         pAdjustMethod = "BH", qvalueCutoff = 0.05)

  go_res <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, ont = "BP",
                     pAdjustMethod = "BH", qvalueCutoff = 0.05)

  if (nrow(kegg_res) > 0) {
    png(file.path("results/enrichment_plots", paste0(cancer_type, "_kegg.png")), width = 900, height = 700)
    print(dotplot(kegg_res, showCategory = 10, title = paste("KEGG:", cancer_type)))
    dev.off()
  }

  if (nrow(go_res) > 0) {
    png(file.path("results/enrichment_plots", paste0(cancer_type, "_go.png")), width = 900, height = 700)
    print(dotplot(go_res, showCategory = 10, title = paste("GO BP:", cancer_type)))
    dev.off()
  }

  write.csv(as.data.frame(kegg_res), file.path("results/enrichment_tables", paste0(cancer_type, "_kegg.csv")), row.names = FALSE)
  write.csv(as.data.frame(go_res), file.path("results/enrichment_tables", paste0(cancer_type, "_go.csv")), row.names = FALSE)

  return(list(kegg = kegg_res, go = go_res))
}


# create directories to store results
dir.create("results/enrichment_plots", recursive = TRUE, showWarnings = FALSE)
dir.create("results/enrichment_tables", recursive = TRUE, showWarnings = FALSE)

# prep data and run
mutexc_list <- readRDS("results/mutexc_list.rds")
enrichment_results <- list()
for (cancer_type in names(mutexc_list)) {
  cat("Running enrichment for", cancer_type, "...\n")
  df <- mutexc_list[[cancer_type]]
  enrichment_results[[cancer_type]] <- enrich_mutexc(cancer_type, df)
}

