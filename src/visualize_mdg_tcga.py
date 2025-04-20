# Visualize mutated driver genes (mdg) in TCGA cohorts

# Get required libraries
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Define helper functions

def get_driver_genes(xlsx_file):
    """Return a list of driver genes from a .xlsx file"""
    driver_genes = pd.read_excel(xlsx_file, sheet_name = "Table S1", skiprows = 3)
    driver_genes = driver_genes.dropna(subset=["Gene"])
    driver_genes = driver_genes["Gene"].str.upper().dropna().unique()
    return driver_genes

def get_mutation_counts_frequency(mutations_directory, driver_genes):
    """
    Return a 1D array with number of mutations in each cohort (mutation file)
    and a DataFrame with mutation frequency for each driver gene per cohort 
    """
    mutation_counts = {}
    mutation_freqs = {}
    # Look in each cohort (mutation file)
    for filename in os.listdir(mutations_directory):
        if filename.endswith(".txt"):
            cohort = filename.split("_")[0]
            filepath = os.path.join(mutations_directory, filename)
            data = pd.read_csv(filepath, sep='\t', index_col=0)
            data.index = data.index.str.upper()
            # Filter the mutation data to just driver genes
            driver_gene_data = data.loc[data.index.isin(driver_genes)]
            # Update mutation counts 
            mutation_counts[cohort] = (driver_gene_data.select_dtypes(include='number').sum(axis=1) > 0).sum()
            # Update mutation frequency
            mutation_freqs[cohort] = (driver_gene_data.select_dtypes(include='number').sum(axis=1) / data.shape[1]).to_dict()
    # Convert results to DataFrames
    mutation_counts_df = pd.Series(mutation_counts).sort_index()
    mutation_freqs_df = pd.DataFrame(mutation_freqs).fillna(0)
    return mutation_counts_df, mutation_freqs_df    

def create_counts_barplot(mutation_counts_df):
    """Create a barplot showing number of driver gene mutations per cohort"""
    plt.figure(figsize=(10, 5))
    mutation_counts_df.plot(kind='bar')
    plt.title("Mutated Driver Gene Count per Cohort")
    plt.ylabel("Number of Mutated Driver Genes")
    plt.xlabel("Cohort")
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig("../results/mdg_counts.png", dpi=300)
    plt.close()
    
def create_freqs_heatmap(mutation_freqs_df):
    """Create a heatmap showing frequence of driver gene mutations per cohort"""
    # Filter to include only genes with freq ≥ 0.2 in at least one cohort
    filtered_df = mutation_freqs_df[mutation_freqs_df.max(axis=1) >= 0.2]
    plt.figure(figsize=(12, 8))
    sns.heatmap(filtered_df, cmap="Reds", linewidths=0.5, linecolor='gray')
    plt.title("Mutation Frequency of Driver Genes (0.2 in at least 1 cohort)")
    plt.ylabel("Driver Genes")
    plt.xlabel("Cohort")
    plt.tight_layout()
    plt.savefig("../results/mdg_freqs.png", dpi=300)
    plt.close()
    
def create_mdg_list_per_cohort(mutation_freqs_df):
    """Save a list of mutated driver genes (mdg) per cohort"""
    rows = []
    for cohort in mutation_freqs_df.columns:
        mutated_genes = mutation_freqs_df[cohort][mutation_freqs_df[cohort] > 0].index.tolist()
        rows.append({
            "Cohort": cohort,
            "Num_Mutated_Driver_Genes": len(mutated_genes),
            "Mutated_Driver_Genes": ", ".join(mutated_genes)
        })
    df = pd.DataFrame(rows)
    df.to_csv("../results/mdg_per_cohort.csv", index=False)

def create_top_mdg(mutation_freqs_df, num_top, freq_status):
    """
    Save the number of cohorts and average mutation frequency for 
    most num_top mutated driver genes
    """
    mdg = mutation_freqs_df[mutation_freqs_df.sum(axis=1) > 0]
    sorting_factor = "avg_freq" if freq_status else "num_cohorts"
    top_mdg = (
        mdg
        .apply(lambda row: pd.Series({
            "num_cohorts": (row > 0).sum(),
            "avg_freq": row[row > 0].mean()
        }), axis=1)
        .sort_values(by=sorting_factor, ascending=False)
        .head(num_top)
    )
    top_mdg.index.name = "mdg"
    top_mdg.to_csv(f"../results/mdg_top_{num_top}_{sorting_factor}.csv")

# Get driver genes from Bailey et al.'s reference list
driver_genes = get_driver_genes("../data/driver_genes.xlsx")
# Get mutation counts and frequency for visualization
mutation_counts_df, mutation_freqs_df = get_mutation_counts_frequency("../data/mutations/", driver_genes)
# Create the driver gene mutations counts and frequency visualizations
create_counts_barplot(mutation_counts_df)
create_freqs_heatmap(mutation_freqs_df)
# Store mutated driver genes per cohort and top 5 mutated driver gene
create_mdg_list_per_cohort(mutation_freqs_df)
create_top_mdg(mutation_freqs_df, 5, 0)
create_top_mdg(mutation_freqs_df, 5, 1)