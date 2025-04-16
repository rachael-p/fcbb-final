import os
import pandas as pd
import matplotlib.pyplot as plt

def plot_co_mutexc(include_UCEC_and_zeros):
    # Paths to result folders
    co_dir = "../results/co"
    mutexc_dir = "../results/mutexc"

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
    
# Get counts of the whole picture, then go deeper
plot_co_mutexc(1)
plot_co_mutexc(0)