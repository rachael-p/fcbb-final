# 🔬 TCGA Driver Gene Mutation Analysis

This project explores the mutation landscape of driver genes across The Cancer Genome Atlas (TCGA) cancer types. It identifies mutation frequency patterns, highlights co-occurring and mutually exclusive mutations, and visualizes gene-gene interactions, particularly focusing on driver genes like TP53 and KRAS.

---

## 📚 Table of Contents

- [Project Structure](#project-structure)
- [Data Sources](#data-sources)
- [Code Overview](#️code-overview)
  - [visualize_mdg_tcga.py](#visualize_mdg_tcgap)
  - [visualize_co_mutexc.py](#visualize_co_mutexcpy)
  - [visualize_tp53_kras.py](#visualize_tp53_kraspy)
- [How to Run](#how-to-run)
- [Biological Insights](#biological-insights)
- [Dependencies](#dependencies)
- [Contact](#contact)

---

## Project Structure

```bash
.
├── data/
│   ├── driver_genes.xlsx            # Reference driver gene list from Bailey et al. 2018
│   └── mutations/                   # TCGA mutation data (gene-by-sample matrices, one .txt file per cohort)
├── results/
│   ├── co/                          # CSV files of co-occurring gene pairs by cohort
│   ├── mutexc/                      # CSV files of mutually exclusive gene pairs by cohort
│   └── [output images and CSVs]     # All generated plots and tables
├── src/
│   ├── visualize_co_mutexc.py       # Co-occurrence / mutual exclusivity counts and gene-bar plots
│   ├── visualize_mdg_tcga.py        # Mutation frequency and driver gene summaries
│   └── visualize_tp53_kras.py       # Network visualizations for TP53 and KRAS
```

---

## Data Sources

- **Driver Genes Reference List** (Bailey et al., 2018):  
  https://pmc.ncbi.nlm.nih.gov/articles/instance/6029450/bin/NIHMS948705-supplement-8.xlsx

- **TCGA Mutation Data** (USC Xena, gene-level binary matrix):  
  https://xenabrowser.net/datapages/  
  Example: *TCGA.LAML.sampleMap/mutation_wustl_gene.gz*

---

## Code Overview

### `visualize_mdg_tcga.py`

- **Input:** `driver_genes.xlsx` and mutation files from `data/mutations/`
- **Output:**
  - `mdg_counts.png`: Barplot of the number of mutated driver genes (mdg) per cohort
  - `mdg_freqs.png`: Heatmap of driver gene mutation frequencies (genes with ≥20% mutations in any cohort)
  - `mdg_per_cohort.csv`: Table of mutated driver genes per cohort
  - `mdg_top_5_avg_freq.csv`: Top 5 driver genes by average mutation frequency
  - `mdg_top_5_num_cohorts.csv`: Top 5 driver genes by number of cohorts mutated in

---

### `visualize_co_mutexc.py`

- **Input:** CSV files in `results/co/` and `results/mutexc/`
- **Output:**
  - `co_mutexc_all.png`: Barplot of co-occurring vs. mutually exclusive pairs (all cohorts)
  - `co_mutexc_non_zero.png`: Barplot of co-occurring vs. mutually exclusive pairs with zero-count and excluding UCEC, because it makes it harder to see the visual trends among the other cohorts
  - `mutation_stacks_all.png`: Stacked barplot of all genes by type of interaction and cohort
  - `mutation_stacks_top_10.png`: Stacked barplot of top 10 most involved genes in mutually exclusive mutations

---

### `visualize_tp53_kras.py`

- **Input:** Mutually exclusive pair files in `results/mutexc/`
- **Output:**
  - `mutexc_network_colored_TP53.png`: TP53-centered network of mutual exclusivity
  - `mutexc_network_colored_KRAS.png`: KRAS-centered network of mutual exclusivity
  - Nodes are sized and colored based on cohort counts, and text contrast is auto-adjusted.

---

## How to Run

> All scripts assume you are running from the `src/` directory and paths are relative.

1. **Visualize mutation frequency and driver gene coverage:**
```bash
python visualize_mdg_tcga.py
```

2. **Visualize co-occurrence vs mutual exclusivity across cohorts:**
```bash
python visualize_co_mutexc.py
```

3. **Generate mutual exclusivity networks for TP53 and KRAS:**
```bash
python visualize_tp53_kras.py
```

---

## Biological Insights

This project enables:
- Identifying highly recurrent driver mutations across cancers
- Comparing mutation exclusivity vs co-occurrence patterns
- Exploring interaction networks of key driver genes (e.g., TP53, KRAS)
- Proposing hypotheses for potential functional interactions or pathway redundancies

---

## Dependencies

- Python 3.x
- `pandas`, `matplotlib`, `seaborn`, `networkx`, `openpyxl`

Install with:
```bash
pip install pandas matplotlib seaborn networkx openpyxl
```

---

## Contact
- Rawan Elshobaky: relshob1@jhu.edu
- Rachael Pei: rpei2@jhu.edu
- Charissa Luk: cluk4@jhu.edu

---

Source for driver gene list: https://pmc.ncbi.nlm.nih.gov/articles/instance/6029450/bin/NIHMS948705-supplement-8.xlsx

Source for TCGA mutation data (TCGA hub in USC Xena): https://xenabrowser.net/datapages/

Github for DISCOVER: https://github.com/NKI-CCB/DISCOVER

Data for synthetic lethal pairs: https://synlethdb.sist.shanghaitech.edu.cn/#/download 