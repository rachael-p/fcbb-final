# install.packages("Rediscover")
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("maftools") 
#install.packages("bigmemory")

library(Rediscover)
library(readr)
library(bigmemory)

data <- read_delim("data/mutations/ACC_mc3_gene_level.txt", delim = "\t", show_col_types = FALSE)
numeric_part <- data[ , -1]
numeric_clean <- as.data.frame(lapply(numeric_part, function(col) as.numeric(as.character(col))))
is_binary <- function(row) {
  all(!is.na(row) & row %in% c(0, 1))
}
binary_only <- numeric_clean[apply(numeric_clean, 1, is_binary), ]
matrix_data <- as.big.matrix(binary_only)
sample_ids <- data$sample
PMA <- getPM(matrix_data)
mymutex = getMutex(A=matrix_data, PM=PMA)