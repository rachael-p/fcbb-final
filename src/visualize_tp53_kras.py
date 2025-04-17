import pandas as pd
import os
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import seaborn as sns

# A function to create a network diagram for a top gene based on 
# the number of cohorts it has mutually exclusive mutations
# with other genes
def visualize_top_genes(mutexc_dir, top_gene):
    top_gene_interactions = {}
    # Load all mutexc files and extract top gene interactions
    for file in os.listdir(mutexc_dir):
        if file.endswith(".csv"):
            cohort = file.split("_")[0]
            df = pd.read_csv(os.path.join(mutexc_dir, file))

            # Standardize gene names
            df["gene1"] = df["gene1"].str.upper()
            df["gene2"] = df["gene2"].str.upper()

            # Filter for top gene interactions
            top_gene_rows = df[(df["gene1"] == top_gene) | (df["gene2"] == top_gene)]
            for _, row in top_gene_rows.iterrows():
                partner = row["gene2"] if row["gene1"] == top_gene else row["gene1"]
                if partner not in top_gene_interactions:
                    top_gene_interactions[partner] = set()
                top_gene_interactions[partner].add(cohort)

    # Convert to DataFrame for plotting
    edges = [(partner, len(cohorts)) for partner, cohorts in top_gene_interactions.items()]
    edges_df = pd.DataFrame(edges, columns=["Gene", "NumCohorts"])
    edges_df = edges_df.sort_values("NumCohorts", ascending=False)

    # Create network
    G = nx.Graph()
    G.add_node(top_gene)
    for gene, count in edges:
        G.add_node(gene)
        G.add_edge(top_gene, gene, weight=count)

    # Assign colors to nodes based on number of cohorts
    # The top gene gets its own neutral color
    node_colors = []
    cmap = plt.cm.Blues
    norm = mcolors.Normalize(vmin=1, vmax=edges_df["NumCohorts"].max())  # normalize to range
    for node in G.nodes():
        if node == top_gene:
            node_colors.append("#cccccc")  # gray for top gene
        else:
            count = len(top_gene_interactions.get(node, []))
            node_colors.append(cmap(norm(count)))

    fig, ax = plt.subplots(figsize=(12, 10))
    pos = nx.spring_layout(G, seed=42)
    edges_weights = [G[u][v]['weight'] for u, v in G.edges()]

    # Draw nodes and edges first (without labels)
    nx.draw_networkx_nodes(G, pos,
                        node_color=node_colors,
                        node_size=1000,
                        ax=ax)
    nx.draw_networkx_edges(G, pos,
                        width=edges_weights,
                        edge_color='gray',
                        ax=ax)

    # Determine label colors based on node brightness
    label_colors = {}
    for node, color in zip(G.nodes(), node_colors):
        if isinstance(color, str):  # e.g., "#cccccc"
            rgb = mcolors.to_rgb(color)
        else:
            rgb = color  # already in RGB tuple format

        brightness = 0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2]
        label_colors[node] = "white" if brightness < 0.5 else "black"

    # Manually draw labels with node-specific font colors
    for node in G.nodes():
        x, y = pos[node]
        ax.text(x, y, node,
                fontsize=8,
                ha='center',
                va='center',
                color=label_colors[node])

    # Add colorbar correctly
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.7)
    cbar.set_label('Number of Cohorts (Mutually Exclusive)', fontsize=12)

    plt.title(f"{top_gene} Mutual Exclusivity Network (Node Color = Cohort Count)", fontsize=14)
    plt.tight_layout()
    plt.savefig(f"../results/{top_gene}_mutexc_network_colored.png", dpi=300)
    plt.close()
    
mutexc_dir = "../results/mutexc"
visualize_top_genes(mutexc_dir, "TP53")
visualize_top_genes(mutexc_dir, "KRAS")