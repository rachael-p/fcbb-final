import pandas as pd
import os
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors

# Load all mutexc files and extract TP53 interactions
mutexc_dir = "../results/mutexc"
tp53_interactions = {}

for file in os.listdir(mutexc_dir):
    if file.endswith(".csv"):
        cohort = file.split("_")[0]
        df = pd.read_csv(os.path.join(mutexc_dir, file))

        # Standardize gene names
        df["gene1"] = df["gene1"].str.upper()
        df["gene2"] = df["gene2"].str.upper()

        # Filter for TP53 interactions
        tp53_rows = df[(df["gene1"] == "TP53") | (df["gene2"] == "TP53")]
        for _, row in tp53_rows.iterrows():
            partner = row["gene2"] if row["gene1"] == "TP53" else row["gene1"]
            if partner not in tp53_interactions:
                tp53_interactions[partner] = set()
            tp53_interactions[partner].add(cohort)

# Convert to DataFrame for plotting
edges = [(partner, len(cohorts)) for partner, cohorts in tp53_interactions.items()]
edges_df = pd.DataFrame(edges, columns=["Gene", "NumCohorts"])
edges_df = edges_df.sort_values("NumCohorts", ascending=False)

# Create network
G = nx.Graph()
G.add_node("TP53")
for gene, count in edges:
    G.add_node(gene)
    G.add_edge("TP53", gene, weight=count)

# Assign colors to nodes based on number of cohorts
# TP53 gets its own neutral color
node_colors = []
cmap = cm.viridis
norm = mcolors.Normalize(vmin=1, vmax=edges_df["NumCohorts"].max())  # normalize to range
for node in G.nodes():
    if node == "TP53":
        node_colors.append("#cccccc")  # gray for TP53
    else:
        count = len(tp53_interactions.get(node, []))
        node_colors.append(cmap(norm(count)))

# Draw the network
plt.figure(figsize=(12, 10))
pos = nx.spring_layout(G, seed=42)
edges_weights = [G[u][v]['weight'] for u, v in G.edges()]

nx.draw(G, pos,
        with_labels=True,
        node_size=800,
        font_size=10,
        width=edges_weights,
        edge_color='gray',
        node_color=node_colors)

# Add colorbar
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = plt.colorbar(sm, shrink=0.7)
cbar.set_label('Number of Cohorts (Mutually Exclusive)', fontsize=12)

plt.title("TP53 Mutual Exclusivity Network (Node Color = Cohort Count)", fontsize=14)
plt.tight_layout()
plt.savefig("../results/tp53_mutexc_network_colored.png", dpi=300)
plt.close()

# Show top 5 mutexc partners
top5 = edges_df.head(5).copy()
top5["Partners"] = top5["Gene"]
top5 = top5.drop(columns="Gene")
print("\nTop 5 TP53 Mutually Exclusive Partners:")
print(top5[["Partners", "NumCohorts"]])
