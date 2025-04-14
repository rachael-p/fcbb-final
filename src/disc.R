#options(repos=c(getOption("repos"), "https://ccb.nki.nl/software/discover/repos/r"))
#install.packages("discover")

library(discover)
library(readr)
library(dplyr)
library(tidyr)

# Read the file
file_path <- "data/mutations/ACC_mc3_gene_level.txt"
data <- read_tsv(file_path)

# Standardize gene names (first column)
gene_names <- toupper(data[[1]])
valid_rows <- !is.na(gene_names) & gene_names != ""
data <- data[valid_rows, ]
gene_names <- gene_names[valid_rows]

# Add gene names as a new column and remove old one
data_clean <- data[ , -1]
data_clean$GENE <- gene_names

# Remove duplicated gene names before setting rownames
data_clean <- data_clean %>%
  distinct(GENE, .keep_all = TRUE) %>%
  relocate(GENE)

# Convert to data.frame and set row names
data_clean_df <- as.data.frame(data_clean)
rownames(data_clean_df) <- data_clean_df$GENE
data_clean_df$GENE <- NULL  # Drop gene column now that it's row names

# Rename columns with cancer type prefix
cancer_type <- strsplit(basename(file_path), "_")[[1]][1]
colnames(data_clean_df) <- paste0(cancer_type, "_", colnames(data_clean_df))

# Filter rows that are not all 0 or all 1, drop NAs, force binary
binary_matrix <- data_clean_df %>%
  filter(rowSums(.) > 0 & rowSums(.) < ncol(.)) %>%
  drop_na() %>%
  mutate(across(everything(), ~ as.integer(pmin(1, .))))

events <- discover.matrix(binary_matrix)
subset <- rowSums(binary_matrix) > 2
result.mutex <- pairwise.discover.test(events[subset, ])
print(result.mutex)
print(as.data.frame(result.mutex))
