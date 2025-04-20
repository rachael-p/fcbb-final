# synthetic lethality: genetic interaction where occurence of two genetic events results in cell death
    # relevant since if two genes are mutually exclusive, they are potential SL candidates - possible that cell might die if both are mutated
    # can explore this hypothesis by checking with synthetic lethality databases
# how can this add new insights: 
    # provides context for SL pairs in databases (not sure if this exists already or not): a lot of pairs are cancer-agnostic so can provide evidence for relevance in specific cancer
        # ex. “These two genes are mutually exclusive in breast cancer, at high statistical significance (FDR 0.01) — this matches an SL pair in SynLethDB.”
        # see if there are therapies / research around these high-confidence mutually exclusive SL pairs already
    # our mutually exclusive pairs may propose new or unstudied SL candidates (future directions to further analyze this) - since SL pairs are often mutually exclusive since don't want the cell to die if trying to develop tumor
    # if find an SL pair from mutexc results but that pair is also co-occuring in a different cohort then it provides tissue-specific biological insights and can guide therapeutic strategies if SL databases don't specify tissue



library(dplyr)

# section 1: see which mutexc pairs match known SL pairs  
# function to append source and confidence info from SL database to my df
find_sl_info <- function(g1, g2, sl_db) {
  match_idx <- which(
    (sl_db$gene1 == g1 & sl_db$gene2 == g2) |
    (sl_db$gene1 == g2 & sl_db$gene2 == g1)
  )
  
  if (length(match_idx) > 0) {
    return(c(sl_db$r.source[match_idx[1]], sl_db$r.statistic_score[match_idx[1]]))
  } else {
    return(c(NA, NA))
  }
}

mutexc_list <- readRDS("results/mutexc_list.rds")
mutexc_list <- mutexc_list[sapply(mutexc_list, nrow) > 0]  # get rid of empty types
mutexc_df <- bind_rows(  # combine all rows in list and add cancer_type column
  lapply(names(mutexc_list), function(cancer_type) {
    df <- mutexc_list[[cancer_type]]
    df$Cancer_Type <- cancer_type
    return(df)
  })
)
print(paste("Total mutexc pairs: ", nrow(mutexc_df)))

sl_db <- read.csv("data/Human_SL.csv", stringsAsFactors = FALSE)
colnames(sl_db)[c(1, 3)] <- c("gene1", "gene2")
sl_db <- sl_db %>%
  mutate(
    gene1 = toupper(trimws(gene1)),
    gene2 = toupper(trimws(gene2)),
  )
sl_db <- sl_db[sl_db$r.statistic_score >= 0.85, ]  # only keeps high confidence SL pairs

mutexc_df$in_SL <- mapply(function(g1, g2) {
  any(
    (sl_db$gene1 == g1 & sl_db$gene2 == g2) |
    (sl_db$gene1 == g2 & sl_db$gene2 == g1)
  )
}, mutexc_df$gene1, mutexc_df$gene2)

SL_hits <- mutexc_df %>% filter(in_SL)
matched_info <- mapply(find_sl_info,
                       SL_hits$gene1,
                       SL_hits$gene2,
                       MoreArgs = list(sl_db = sl_db))
SL_hits$SL_source <- matched_info[1, ]
SL_hits$SL_score <- as.numeric(matched_info[2, ])
SL_hits$in_SL <- NULL

print(paste("Number of SL hits:", nrow(SL_hits)))
write.csv(SL_hits, "results/SL_hits_high_conf.csv", row.names = FALSE)



# section 2: trying to identify potential SL target
gene_counts <- table(c(mutexc_df$gene1, mutexc_df$gene2))  # get recurrence counts - if gene appears more freq then may be of greater interest
mutexc_df_filt <- mutexc_df %>%
  filter(in_SL == FALSE) %>%
  mutate(
    gene1_freq = gene_counts[gene1],
    gene2_freq = gene_counts[gene2],
    pair_recurrence = gene1_freq + gene2_freq
  )
mutexc_ranked <- mutexc_df_filt %>%  # rank it based on statistical significance and recurrence 
  arrange(q.value) %>%
  mutate(
    rank_score = (1 / q.value) + pair_recurrence
  ) %>%
  arrange(desc(rank_score))
top_10 <- mutexc_ranked %>% slice_head(n = 10)
write.csv(top_10, "results/potential_SL_targets.csv", row.names = FALSE)
