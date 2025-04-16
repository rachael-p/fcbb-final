import os
import pandas as pd
import matplotlib.pyplot as plt

# Defining paths
co_dir = "../results/co"
mutexc_dir = "../results/mutexc"

def plot_co_mutexc(co_dir, mutexc_dir, include_UCEC_and_zeros):
    # Defining variables
    co_counts = {}
    mutexc_counts = {}

    # Get cohort names from filenames
    for file in os.listdir(co_dir):
        if file.endswith(".csv"):
            cohort = file.split("_")[0]
            co_df = pd.read_csv(os.path.join(co_dir, file))
            mutexc_df = pd.read_csv(os.path.join(mutexc_dir, f"{cohort}_mutexc.csv"))
            num_mutexc = len(mutexc_df)
            if not include_UCEC_and_zeros:
                if cohort != "UCEC" and num_mutexc != 0: # Ignore UCEC and 0 mutually exclusive mutations
                    co_counts[cohort] = len(co_df)
                    mutexc_counts[cohort] = num_mutexc
            else:
                co_counts[cohort] = len(co_df)
                mutexc_counts[cohort] = num_mutexc
            
    # Create DataFrame for plotting
    df = pd.DataFrame({
        "Cohort": list(co_counts.keys()),
        "Co-occurring": [co_counts[c] for c in co_counts],
        "Mutually Exclusive": [-mutexc_counts[c] for c in mutexc_counts]  # negative for diverging plot
    }).sort_values("Cohort")
    
    # Adjust title based on function parameters
    if include_UCEC_and_zeros:
        title = "Counts of Co-occurring vs Mutually Exclusive Driver Gene Pairs per Cohort"
        output_name = "co_mutexc_all"
    else:
        title = "Non-Zero Counts of Co-occurring vs Mutually Exclusive Driver Gene Pairs per Cohort"
        output_name = "co_mutexc_non_zero"

    # Plot
    plt.figure(figsize=(14, 6))
    bar_width = 0.4
    x = range(len(df))

    plt.bar(x, df["Co-occurring"], width=bar_width, label="Co-occurring", color="skyblue")
    plt.bar(x, df["Mutually Exclusive"], width=bar_width, label="Mutually Exclusive", color="lightcoral")

    plt.axhline(0, color='black', linewidth=0.8)
    plt.xticks(x, df["Cohort"], rotation=45, ha='right')
    plt.ylabel("Number of Gene Pairs")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"../results/{output_name}.png", dpi=300)
    
def count_gene_pairs(folder):
    gene_counts = {}
    for fname in os.listdir(folder):
        if fname.endswith(".csv"):
            cohort = fname.split("_")[0]
            df = pd.read_csv(os.path.join(folder, fname))
            for gene in pd.concat([df["gene1"], df["gene2"]]):
                if gene not in gene_counts:
                    gene_counts[gene] = {}
                gene_counts[gene][cohort] = gene_counts[gene].get(cohort, 0) + 1
    return pd.DataFrame.from_dict(gene_counts, orient='index').fillna(0).astype(int)

def plot_mutation_stacks(co_dir, mutexc_dir, is_top_10):
    # Count occurrences
    co_df = count_gene_pairs(co_dir)
    mutexc_df = count_gene_pairs(mutexc_dir)
    
    if is_top_10:
        # Count occurrences
        co_df = count_gene_pairs(co_dir)
        mutexc_df = count_gene_pairs(mutexc_dir)

        # Compute total involvement
        co_totals = co_df.sum(axis=1)
        mutexc_totals = mutexc_df.sum(axis=1)

        # Get top 10 genes from each
        top_co_genes = co_totals.sort_values(ascending=False).head(10).index
        top_mutexc_genes = mutexc_totals.sort_values(ascending=False).head(10).index

        plot_genes = sorted(set(top_co_genes).union(set(top_mutexc_genes)))
        label = "top_10"
    else:
        plot_genes = set(co_df.index).union(mutexc_df.index)
        label = "all"
    # Keep only genes that appear at least once in either
    co_df = co_df.reindex(plot_genes, fill_value=0)
    mutexc_df = mutexc_df.reindex(plot_genes, fill_value=0)

    # Plotting
    # cohorts = sorted(set(co_df.columns).union(mutexc_df.columns))
    colors = plt.cm.tab20.colors  # color palette

    plt.figure(figsize=(14, 10))
    x = range(len(plot_genes))
    gene_names = list(plot_genes)
    
    # CO-OCCURRING: Stack bars upward
    bottoms = [0] * len(plot_genes)
    co_cohorts = co_df.columns
    for i, co_cohort in enumerate(co_cohorts):
        heights = co_df[co_cohort].tolist()
        plt.bar(x, heights, bottom=bottoms, label=f"{co_cohort} (co)", color=colors[i % len(colors)])
        bottoms = [b + h for b, h in zip(bottoms, heights)]
        
    # MUTEXC: Stack bars downward
    bottoms = [0] * len(plot_genes)
    mutexc_cohorts = mutexc_df.columns
    for i, mutexc_cohort in enumerate(mutexc_cohorts):
        heights = (-mutexc_df[mutexc_cohort]).tolist()
        plt.bar(x, heights, bottom=bottoms, label=f"{mutexc_cohort} (mutexc)", color=colors[i % len(colors)], alpha=0.4)
        bottoms = [b + h for b, h in zip(bottoms, heights)]

    plt.xticks(x, gene_names, rotation=90, fontsize=8)
    plt.ylabel("Gene Involvement Count")
    plt.title("Driver Gene Involvement in Co-occurring and Mutually Exclusive Mutations (per Cohort)")
    plt.axhline(0, color="black")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(f"../results/{label}_mutation_stacks.png", dpi=300)
    
# Get counts of the whole picture, then go deeper
# plot_co_mutexc(co_dir, mutexc_dir, 1)
# plot_co_mutexc(co_dir, mutexc_dir, 0)

plot_mutation_stacks(co_dir, mutexc_dir, 1)
# plot_mutation_stacks(co_dir, mutexc_dir, 0)