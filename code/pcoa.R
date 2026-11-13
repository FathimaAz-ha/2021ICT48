# ============================================================
# 02_pcoa.R
# Clean body-site labels, filter rare OTUs, compute
# Bray-Curtis dissimilarity, and plot PCoA.
# ============================================================

library(MicrobeDS)
library(phyloseq)
library(tidyverse)

data(HMPv13)
ps <- HMPv13

# ---- 1. Clean UBERON labels and add a broad region column ----
meta <- as.data.frame(sample_data(ps))
meta$body_site_clean <- gsub("^UBERON:", "", meta$body_site)
meta$body_site_clean <- gsub("_", " ", meta$body_site_clean)

meta$body_region <- case_when(
  grepl("skin",                      meta$body_site_clean) ~ "Skin",
  grepl("vagina|introitus|fornix",   meta$body_site_clean) ~ "Vaginal",
  grepl("feces",                     meta$body_site_clean) ~ "Gut",
  grepl("nares",                     meta$body_site_clean) ~ "Nasal",
  TRUE                                                     ~ "Oral"
)

sample_data(ps) <- meta

# sanity check: 5 regions expected
table(meta$body_region)

# ---- 2. Filter rare OTUs (keep taxa with total count > 10) ----
keep <- taxa_sums(otu_table(ps)) > 10
ps_f <- prune_taxa(keep, ps)
cat("OTUs before filtering:", ntaxa(ps),  "\n")
cat("OTUs after  filtering:", ntaxa(ps_f), "\n")

# ---- 3. Convert counts to relative abundance ----
ps_rel <- transform_sample_counts(ps_f, function(x) x / sum(x))

# ---- 4. Bray-Curtis dissimilarity (may take 1-5 min) ----
bray <- phyloseq::distance(ps_rel, method = "bray")

# ---- 5. PCoA ordination ----
ord <- ordinate(ps_rel, method = "PCoA", distance = bray)

# ---- 6. Variance explained by PC1-PC3 ----
eig <- ord$values$Eigenvalues
var_exp <- round(100 * eig / sum(eig[eig > 0]), 1)
cat("Variance explained by PC1-PC3 (%):", var_exp[1:3], "\n")

# ---- 7. Plot: broad region ----
p1 <- plot_ordination(ps_rel, ord, color = "body_region") +
  geom_point(size = 1.5, alpha = 0.7) +
  theme_minimal() +
  labs(title    = "PCoA of HMP 16S profiles (Bray-Curtis)",
       subtitle = "Coloured by broad body region",
       color    = "Region")

ggsave("figures/pcoa_region.png", p1, width = 7, height = 5, dpi = 300)
print(p1)

# ---- 8. Plot: fine body sites ----
p2 <- plot_ordination(ps_rel, ord, color = "body_site_clean") +
  geom_point(size = 1.2, alpha = 0.7) +
  theme_minimal() +
  labs(title    = "PCoA of HMP 16S profiles (Bray-Curtis)",
       subtitle = "Coloured by body site",
       color    = "Body site")

ggsave("figures/pcoa_site.png", p2, width = 8, height = 5, dpi = 300)
print(p2)