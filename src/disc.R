#options(repos=c(getOption("repos"), "https://ccb.nki.nl/software/discover/repos/r"))
#install.packages("discover")

library(discover)
library(readr)
library(dplyr)
library(tidyr)

# set up data reading
driver_ref <- read_csv("results/mdg_per_cohort.csv")
directory <- "data/mutations"
file_list <- list.files(directory, pattern = "\\.txt$", full.names = TRUE)
result_list <- list()


for (i in seq_along(file_list)) {
    file_path <- file_list[i]
    cancer_type <- strsplit(basename(file_path), "_")[[1]][1]

    data <- read_tsv(file_path, show_col_types = FALSE)

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

    # Get driver gene list for this cancer type
    driver_row <- driver_ref[driver_ref$Cohort == cancer_type, ]

    # Parse comma-separated string of gene names
    if (nrow(driver_row) > 0 && !is.na(driver_row$Mutated_Driver_Genes)) {
        driver_genes <- strsplit(driver_row$Mutated_Driver_Genes, ",")[[1]] %>%
                        toupper() %>%
                        trimws()  # ensure upper case and no spaces

        # Subset binary matrix to only those driver genes (and make sure they exist in the matrix)
        driver_subset <- intersect(driver_genes, rownames(binary_matrix))

        if (length(driver_subset) > 1) {
            events <- discover.matrix(binary_matrix)
            num_genes <- length(driver_subset)
            num_comparisons <- choose(num_genes, 2)
            cat(sprintf("Running %d pairwise comparisons for %s...\n", num_comparisons, cancer_type))
            result.mutex <- pairwise.discover.test(events[driver_subset, ])
            result_list[[cancer_type]] <- as.data.frame(result.mutex)
            print(paste("Finished:", cancer_type, "with", length(driver_subset), "driver genes"))
        } else {
            print(paste("Skipping:", cancer_type, "- only 1 driver gene found"))
        }
    } else {
        print(paste("Skipping:", cancer_type, "- no driver genes listed"))
    }
}


print(result_list[["ACC"]])


# TODO: 
# save results to a file instead of printing
# ask to translate from this to base R 
# look at the results to process them - not sure if this is testing for exc or co 
    # what is the statistical significance 
    # how to pick which pairs we want to look at (record number) - maybe try groupwise or see if some have significance like pathways?
    # look at what the discover paper itself did 
    # also for the paper see if other studies have used discover vs. traditional tests
    # how can it help achieve goal of informing selection of therapeutic targets like synthetic lethal strategies
    # maybe go a bit deeper
# then find biological hypotheses
# also start looking at survival data 
    # how co or exc presence influence patient survival since now that we have a list of significant ones we can check if patients have them (but are the samples the same)
    # store the significant ones in a csv 