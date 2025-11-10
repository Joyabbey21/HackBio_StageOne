"""Surprise_task(stage1)"""
#Gene Expression (Heatmap & Volcano Plot)
#Part A – Gene Expression Analysis

# Import libraries
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Step 1: Load the dataset from the GitHub link
url = "https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/hbr_uhr_top_deg_normalized_counts.csv"
df_gene = pd.read_csv(url)

df_gene.head(13)

# Step 2: Set gene names as index
df_heatmap = df_gene.set_index(df_gene.columns[0])  # The first column is gene names

# Step 3: Create clustered heatmap
plt.figure(figsize=(10, 8))
sns.clustermap(
    df_heatmap,
    cmap="Blues",
    standard_scale=1,
    linewidths=0.5,
    figsize=(10, 6)
)

#b. Volcano Plot
import numpy as np

# Step 1: Load the dataset
url_link = "https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/hbr_uhr_deg_chr22_with_significance.csv"

df_deg = pd.read_csv(url_link)

df_deg.head()

# Step 2: Clean and filter data
# Ensure no infinite or zero values for padj (adjusted p-value)
df = df_deg.replace([np.inf, -np.inf], np.nan)
df = df_deg.dropna(subset=["log2FoldChange", "PAdj"])
df = df[df["PAdj"] > 0]  # avoid padj = 0 to compute −log10

df_deg.head(20)

df_deg.tail(20)

# Step 3: Define thresholds for significance
log2fc_threshold = 1
padj_threshold = 0.05

# Step 4: Categorize each gene
def categorize(row):
    if row["PAdj"] < padj_threshold and row["log2FoldChange"] > log2fc_threshold:
        return "Upregulated"
    elif row["PAdj"] < padj_threshold and row["log2FoldChange"] < -log2fc_threshold:
        return "Downregulated"
    else:
        return "Not Significant"

df["Regulation"] = df.apply(categorize, axis=1)

# Step 5: Map colours
colour_map = {
    "Upregulated": "green",
    "Downregulated": "blue",
    "Not Significant": "orange"
}
df["Colour"] = df["Regulation"].map(colour_map)

# Step 6: Compute y-axis value (−log10 adjusted p-value)
df["negLog10Padj"] = -np.log10(df["PAdj"])

#Plot log2FoldChange vs log10(Padj) from the DEG results.
#Step 7: Create the plot

plt.figure(figsize=(10, 7))
plt.scatter(
    df["log2FoldChange"],
    df["negLog10Padj"],
    c=df["Colour"],
    alpha=0.7,
    edgecolor='none'
)

#Step 8: Add dashed lines for fold-change thresholds and significance threshold
plt.axvline(x=log2fc_threshold, color="black", linestyle="--", linewidth=1)
plt.axvline(x=-log2fc_threshold, color="black", linestyle="--", linewidth=1)
plt.axhline(y=-np.log10(padj_threshold), color="black", linestyle="--", linewidth=1)

#Step 9: Labels and title
plt.title("Volcano Plot of log2FoldChange vs log10(Padj)", fontsize=14)
plt.xlabel("log2(Fold Change)")
plt.ylabel("-log10(Adjusted P-value)")

plt.legend(
    handles=[
        plt.Line2D([0], [0], marker='o', color='w', label='Upregulated', markerfacecolor='green', markersize=8),
        plt.Line2D([0], [0], marker='o', color='w', label='Downregulated', markerfacecolor='blue', markersize=8),
        plt.Line2D([0], [0], marker='o', color='w', label='Not Significant', markerfacecolor='orange', markersize=8)
    ],
    title="Gene Regulation",
    loc='upper right'
)

plt.tight_layout()
plt.show()

#Part B
#Breast Cancer Diagnostic Data (Correlation & Scatter/Density Plots)

url_link2 = "https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/data-3.csv"

df_brst = pd.read_csv(url_link2)

df_brst.head()

df_brst.tail()

df_brst.info()

#c.Scatter Plot (radius vs texture)
plt.scatter(df_brst['radius_mean'], df_brst['texture_mean'],
            c=df_brst['diagnosis'].map({'M': 'blue', 'B': 'orange'}),
            alpha=0.7, edgecolors='k')

#Add labels and title
plt.xlabel("Radius Mean")
plt.ylabel("Texture Mean")
plt.title("Scatter Plot of Texture Mean vs Radius Mean by Diagnosis")

#Add legend
for diag, color in colors.items():
    plt.scatter([], [], c=color, label=f'{diag} ({ "Malignant" if diag == "M" else "Benign" })')
plt.legend(title="Diagnosis")

#Show plot
plt.show()


#d. Correlation Heatmap
# Clean column names (replace spaces with underscores for convenience)
df.columns = df.columns.str.strip().str.replace(' ', '_')

# Select the six key features
features = ['radius_mean', 'texture_mean', 'perimeter_mean', 'area_mean', 'smoothness_mean', 'compactness_mean']

# Compute correlation matrix
corr_matrix = df_brst[features].corr()

# Display correlation matrix
print("Correlation Matrix:\n")
print(corr_matrix)

# Plot heatmap
plt.figure(figsize=(8, 6))
sns.heatmap(corr_matrix, annot=True, cmap="Blues", fmt=".1f", linewidths=0.5)
plt.title("Correlation Heatmap of Key Features", fontsize=14)
plt.show()

#Create color mapping for diagnosis
colors = {'M': 'blue', 'B': 'orange'}

#e. Scatter Plot (smoothness vs compactness)
#Create scatter plot

plt.figure(figsize=(8, 6))
plt.scatter(df_brst['smoothness_mean'], df_brst['compactness_mean'],
            c=df_brst['diagnosis'].map(colors), alpha=0.7, edgecolor='k')

#Add gridlines, labels and title
plt.grid(True, linestyle='--', alpha=0.6)
plt.xlabel("Smoothness Mean", fontsize=12)
plt.ylabel("Compactness Mean", fontsize=12)
plt.title("Scatter Plot of Compactness vs Smoothness (Colored by Diagnosis)", fontsize=14)


#Add legend
for diag, color in colors.items():
    plt.scatter([], [], c=color, label=f'{diag} ({ "Malignant" if diag == "M" else "Benign" })')
plt.legend(title="Diagnosis")
plt.tight_layout()
plt.show()

#f. Density Plot (area distribution)
#Create the KDE plot
plt.figure(figsize=(8, 6))
sns.kdeplot(data=df_brst[df_brst['diagnosis'] == 'M'], x='area_mean', label='Malignant (M)', color='blue', fill=True, alpha=0.5)
sns.kdeplot(data=df_brst[df_brst['diagnosis'] == 'B'], x='area_mean', label='Benign (B)', color='orange', fill=True, alpha=0.5)

plt.xlabel("Area Mean", fontsize=12)
plt.ylabel("Density", fontsize=12)
plt.title("Density Plot of Area Mean by Diagnosis", fontsize=14)
plt.legend(title="Diagnosis")

# Show gridlines for better readability
plt.grid(True, linestyle='--', alpha=0.6)

plt.tight_layout()
plt.show()

#Write a python function for translating DNA to protein
def dna_to_protein(dna_seq):
    """Translate a DNA sequence into a protein sequence."""
    # Standard codon table
    codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }

    protein = ""
    # Translate in codons of 3 bases
    for i in range(0, len(dna_seq) - 2, 3):
        codon = dna_seq[i:i+3].upper()
        amino_acid = codon_table.get(codon, 'X')  # 'X' for unknown codons
        protein += amino_acid
    return protein

dna_sequence = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
print("Protein sequence:", dna_to_protein(dna_sequence))


#Write a python function for calculating the hamming distance between your slack username (e.g josoga) and twitter/X (joseph) handle (synthesize one if you don’t have one). Feel free to pad it with extra words if they are not of the same length.
def hamming_distance(str1, str2):
    """Calculate Hamming distance between two strings."""
    # Pad shorter string with underscores
    max_len = max(len(str1), len(str2))
    str1 = str1.ljust(max_len, '_')
    str2 = str2.ljust(max_len, '_')

    # Compare character by character
    distance = sum(ch1 != ch2 for ch1, ch2 in zip(str1, str2))
    return distance

# Example usage:
slack_username = "ABIODUN JOY"
twitter_handle = "Joy_abbey_"
print(f"Hamming distance between '{slack_username}' and '{twitter_handle}':",
      hamming_distance(slack_username, twitter_handle))