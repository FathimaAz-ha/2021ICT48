# Body-site specificity in the human microbiome: can unsupervised clustering recover body-site identity?

**Student ID:** 2021ICT48

## 1. Abstract

The human body hosts diverse microbial communities that vary significantly by anatomical site. This project investigates whether 16S rRNA profiles are distinct enough to identify body sites using unsupervised machine learning. We analysed the Human Microbiome Project dataset (3,285 samples, 15 body sites) using Bray-Curtis dissimilarity, Principal Coordinates Analysis (PCoA), and Ward hierarchical clustering. PCoA revealed clear separation of gut and vaginal microbiomes, while oral and skin samples overlapped. Clustering into five groups achieved a mean purity of 66.6%, heavily driven by the class imbalance of oral samples. The results confirm that 16S data can broadly differentiate body regions but struggles with fine-scale anatomical subsites and imbalanced datasets.

## 2. Introduction

The human microbiome plays a critical role in health and disease, influencing metabolism, immunity, and protection against pathogens. Different body sites—such as the gut, skin, mouth, and vagina—provide distinct ecological niches, resulting in varied microbial compositions. Understanding whether these communities are distinct enough to be identified without prior labels is essential for clinical diagnostics, forensic microbiology, and source-tracking in outbreaks.

This project aims to answer: **Can unsupervised clustering of 16S rRNA data recover human body-site identity?** We hypothesise that due to differing environmental pressures (e.g., oxygen levels, pH, moisture), broad body regions will form distinct clusters. However, closely related sites (e.g., within the oral cavity) may be harder to distinguish due to shared salivary environments.

## 3. Data

The data used is the `HMPv13` dataset from the Human Microbiome Project (HMP), accessed via the `MicrobeDS` R package. The dataset contains 3,285 samples across 15 anatomical sites and 5,830 Operational Taxonomic Units (OTUs).

The original dataset uses UBERON anatomy ontology labels (e.g., `UBERON:feces`). For this analysis, labels were cleaned and samples were grouped into five broad regions: \* **Oral:** 1,818 samples (8 subsites including tongue, tooth, palate) \* **Skin:** 779 samples (2 subsites: elbow, ear) \* **Vaginal:** 291 samples (3 subsites) \* **Gut (Feces):** 210 samples \* **Nasal:** 187 samples

This data is publicly available through the HMP and is openly licensed for research use.

## 4. Methods

All analysis was performed in R (v4.6.2) using the `phyloseq` and `tidyverse` packages.

**Preprocessing:** Rare OTUs (those with fewer than 10 total reads across all samples) were filtered out to reduce noise, leaving 3,130 OTUs. Counts were then converted to relative abundance.

**Distance and Ordination:** Bray-Curtis dissimilarity was used to quantify the distance between samples. This metric is well-suited for ecological abundance data because it handles the absence of species well and does not assume linear relationships. Principal Coordinates Analysis (PCoA) was used for 2D visualisation of the distance matrix.

**Clustering:** Unsupervised hierarchical clustering was performed using Ward.D2 linkage, which minimises within-cluster variance. The resulting dendrogram was cut into 5 clusters (matching the 5 broad regions). Cluster purity was assessed by calculating the proportion of the dominant region within each cluster.

## 5. Results

### 5.1 Exploratory Analysis (PCoA)

The PCoA plot (Figure 1) revealed that the first two axes explained 15.3% and 8.2% of the variance, respectively. Despite the high dimensionality of microbiome data, biological structure was visually apparent: \* **Vaginal** samples formed a tight, distinct cluster on the bottom left. \* **Gut** samples formed a separate cluster at the bottom. \* **Skin** samples spread along the right side. \* **Oral** and **Nasal** samples overlapped heavily in the upper-left quadrant.

*(Insert Figure 1 here: `figures/pcoa_region.png`)* *Figure 1: PCoA of HMP 16S profiles using Bray-Curtis dissimilarity.*

### 5.2 Unsupervised Clustering

Hierarchical clustering into 5 groups yielded a mean purity of **66.6%**. The cluster composition plot (Figure 2) showed that three out of five clusters were almost entirely composed of Oral samples (99.1%, 97.4%, and 98.6%). Cluster 3 was predominantly Skin (73.1%). Cluster 1 was the most mixed, dominated by Skin (35.4%) but containing the majority of Gut and Vaginal samples.

*(Insert Figure 2 here: `figures/cluster_composition_5.png`)* *Figure 2: Composition of 5 unsupervised clusters.*

## 6. Discussion

The results support the hypothesis that broad body regions have distinct microbial signatures, but highlight significant methodological challenges.

The high overlap between Oral and Nasal samples is biologically plausible, as both are part of the upper respiratory tract exposed to similar environmental microbes. Conversely, the Gut and Vaginal microbiomes were highly distinct, which is expected given their vastly different physiological environments (anaerobic vs. acidic).

The 66.6% purity for 5 clusters is largely driven by **class imbalance**. Because Oral samples constitute 55% of the dataset, the clustering algorithm split them into three separate clusters to minimise variance. This starved the smaller groups (Gut, Vaginal) of their own clusters, causing them to be absorbed into the heterogeneous Cluster 1. This demonstrates a known limitation of standard hierarchical clustering when applied to imbalanced ecological data.

At a finer resolution of 15 clusters (matching the 15 specific body sites), mean purity dropped to 31.4%. This reflects the difficulty of distinguishing between fine-scale oral subsites (e.g., tongue, tooth, palate) using 16S data alone, as they share the same salivary environment.

**Limitations:** The study relies on 16S rRNA sequencing, which provides taxonomic resolution only to the genus level, potentially missing species-level differences. Furthermore, the PCoA explains a low percentage of total variance, suggesting that factors like host genetics, diet, and temporal variation are not captured in this analysis.

## 7. Conclusion and Future Work

This study demonstrates that unsupervised machine learning can successfully identify broad human body regions from 16S microbiome data, achieving 66.6% purity. However, class imbalance and fine-scale anatomical similarities limit its accuracy.

Future work could apply supervised machine learning (e.g., Random Forest) to improve classification accuracy, or use metagenomic sequencing for higher taxonomic resolution. Additionally, techniques like SMOTE (Synthetic Minority Over-sampling Technique) could be used to address the class imbalance prior to clustering.

## 8. Reproducibility

To reproduce this analysis: 1. Clone this repository. 2. Open the `.Rproj` file in RStudio. 3. Install the required packages: `MicrobeDS`, `phyloseq`, `tidyverse`. 4. Knit the file `code/Analysis.Rmd` to generate the full HTML report and all figures.

## 9. References

- Battaglia, T. et al. (2015). MicrobeDS: A curated data package for microbiome analysis.
- McMurdie, P.J. and Holmes, S. (2013). phyloseq: An R package for reproducible interactive analysis and graphics of microbiome census data. PLoS ONE.
- Human Microbiome Project Consortium (2012). Structure, function and diversity of the healthy human microbiome. Nature.
