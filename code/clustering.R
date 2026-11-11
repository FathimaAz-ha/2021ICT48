# ============================================================
# 03_clustering.R
# Unsupervised hierarchical clustering and purity assessment
# ============================================================

library(MicrobeDS)
library(phyloseq)
library(tidyverse)

# Re-load and prep data (same as Day 2)
data(HMPv13)
ps <- HMPv13
meta <- as.data.frame(sample_data(ps))
meta$body_site_clean <- gsub("^UBERON:", "", meta$body_site)
meta$body_site_clean <- gsub("_", " ", meta$body_site_clean)
meta$body_region <- case_when(
  grepl("skin",                    meta$body_site_clean) ~ "Skin",
  grepl("vagina|introitus|fornix", meta$body_site_clean) ~ "Vaginal",
  grepl("feces",                   meta$body_site_clean) ~ "Gut",
  grepl("nares",                   meta$body_site_clean) ~ "Nasal",
  TRUE                                                   ~ "Oral"
)
sample_data(ps) <- meta

keep <- taxa_sums(otu_table(ps)) > 10
ps_f <- prune_taxa(keep, ps)
ps_rel <- transform_sample_counts(ps_f, function(x) x / sum(x))

# Compute distance again
bray <- phyloseq::distance(ps_rel, method = "bray")

# ---- 1. Hierarchical clustering ----
# Ward.D2 minimizes variance within clusters (good for ecological data)
hc <- hclust(bray, method = "ward.D2")

# ---- 2. Cut tree into 5 clusters (matching the 5 broad regions) ----
clusters_5 <- cutree(hc, k = 5)

# Add cluster assignments to the metadata
meta$cluster_5 <- factor(clusters_5)

# ---- 3. Compute Purity (Region level) ----
# Purity = for each cluster, what % belongs to the dominant region?
purity_5 <- meta %>%
  count(cluster_5, body_region) %>%
  group_by(cluster_5) %>%
  mutate(prop = n / sum(n)) %>%
  slice_max(prop, n = 1) %>%
  ungroup()

# Overall purity = average of the max proportions weighted by cluster size
overall_purity_5 <- sum(purity_5$n * purity_5$prop) / nrow(meta)
cat("Mean Purity (5 clusters, broad region):", round(overall_purity_5 * 100, 1), "%\n")

# Print the cluster composition table
print(purity_5)

# ---- 4. Plot the composition of the 5 clusters ----
p_bar <- meta %>%
  count(cluster_5, body_region) %>%
  ggplot(aes(x = cluster_5, y = n, fill = body_region)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_minimal() +
  labs(title = "Composition of 5 Unsupervised Clusters",
       x = "Cluster", y = "Proportion",
       fill = "Broad Region")

ggsave("figures/cluster_composition_5.png", p_bar, width = 6, height = 4, dpi = 300)
print(p_bar)

# ---- 5. Repeat at finer resolution (15 clusters matching 15 body sites) ----
clusters_15 <- cutree(hc, k = 15)
meta$cluster_15 <- factor(clusters_15)

purity_15 <- meta %>%
  count(cluster_15, body_site_clean) %>%
  group_by(cluster_15) %>%
  mutate(prop = n / sum(n)) %>%
  slice_max(prop, n = 1) %>%
  ungroup()

overall_purity_15 <- sum(purity_15$n * purity_15$prop) / nrow(meta)
cat("Mean Purity (15 clusters, fine body site):", round(overall_purity_15 * 100, 1), "%\n")