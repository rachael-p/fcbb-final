#options(repos=c(getOption("repos"), "https://ccb.nki.nl/software/discover/repos/r"))
#install.packages("discover")

library(discover)
library(readr)
library(tidyr)
library(dplyr)

# import data
directory <- "../data/mutations"
file_list <- list.files(directory, pattern = "\\.txt$", full.names = TRUE)
mutexc_list <- list()
co_list <- list()
driver_ref <- read_csv("../results/mdg_per_cohort.csv", show_col_types = FALSE)
driver_ref$Mutated_Driver_List <- strsplit(driver_ref$Mutated_Driver_Genes, ",")
driver_ref$Mutated_Driver_List <- lapply(driver_ref$Mutated_Driver_List, function(x) toupper(trimws(x)))

for (file_path in file_list) {
    cancer_type <- strsplit(basename(file_path), "_")[[1]][1]
    data <- read_tsv(file_path, show_col_types = FALSE)

    # pre-process and clean data
    data[[1]] <- toupper(data[[1]])
    valid_rows <- !is.na(data[[1]]) & data[[1]] != ""  # logical vector
    data <- data[valid_rows, ]  # only keep rows with valid gene names
    names(data)[1] <- "GENE"
    data <- data %>%  # combine rows if have same gene twice in same file by ORing them
            group_by(GENE) %>%
            summarise(across(everything(), ~ as.integer(any(. == 1))), .groups = "drop")

    binary_matrix <- as.data.frame(data)
    rownames(binary_matrix) <- binary_matrix$GENE  # turn gene col into rownames instead
    binary_matrix$GENE <- NULL
    colnames(binary_matrix) <- paste0(cancer_type, "_", colnames(binary_matrix))  # prepend cancer type to sample names
    
    row_totals <- rowSums(binary_matrix)
    binary_matrix <- binary_matrix[row_totals > 0 & row_totals < ncol(binary_matrix), ]  # removes rows with all 0s or 1s (not useful for exc/co)
    binary_matrix <- binary_matrix %>%
        drop_na() %>%
        mutate(across(everything(), ~ as.integer(pmin(1, .))))  # makes sures all values are binary (0 or 1)


    # run tests
    driver_row <- driver_ref[driver_ref$Cohort == cancer_type, ]
    if (driver_row$Num_Mutated_Driver_Genes > 0) {
        driver_genes <- driver_row$Mutated_Driver_List[[1]]
        driver_subset <- intersect(driver_genes, rownames(binary_matrix))  # only include drivers in matrix

        if (length(driver_subset) > 1) {
            disc_matrix <- discover.matrix(binary_matrix)
            cat(sprintf("Running %d pairwise comparisons for %s \n", choose(length(driver_subset), 2), cancer_type))
            mutexc_result <- pairwise.discover.test(disc_matrix[driver_subset, ], alternative = "less")
            mutexc_list[[cancer_type]] <- as.data.frame(mutexc_result)  # only stores significant results, default FDR is 1%
            co_result <- pairwise.discover.test(disc_matrix[driver_subset, ], alternative = "greater")
            co_list[[cancer_type]] <- as.data.frame(co_result)
        } else {
            print(paste("Skipping:", cancer_type, "- only 1 driver gene found"))
        }
    } else {
        print(paste("Skipping:", cancer_type, "- no driver genes"))
    }
}

# save results to compressed binary files so they can be loaded into other scripts 
saveRDS(mutexc_list, "../results/mutexc_list.rds")
saveRDS(co_list, "../results/co_list.rds")

# save results to csv files (only creates a csv if there are significant pair results - otherwise cohort name goes added to txt file)
file.create("../results/no_mutexc.txt")
dir.create("../results/mutexc", recursive = TRUE, showWarnings = FALSE)
for (cancer_type in names(mutexc_list)) {
  if (nrow(mutexc_list[[cancer_type]]) > 0) {
    write.csv(mutexc_list[[cancer_type]],
        file = file.path("../results/mutexc", paste0(cancer_type, "_mutexc.csv")),
        row.names = FALSE)
  } else {
    write(cancer_type, file = "../results/no_mutexc.txt", append = TRUE)
  }
}

file.create("../results/no_co.txt")
dir.create("../results/co", recursive = TRUE, showWarnings = FALSE)
for (cancer_type in names(co_list)) {
  if (nrow(co_list[[cancer_type]]) > 0) {
    write.csv(co_list[[cancer_type]],
        file = file.path("../results/co", paste0(cancer_type, "_co.csv")),
        row.names = FALSE)
  } else {
    write(cancer_type, file = "../results/no_co.txt", append = TRUE)
  }
}

