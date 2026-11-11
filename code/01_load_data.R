# Load the packages
library(MicrobeDS)
library(phyloseq)
library(tidyverse)

# Load the HMPv13 dataset (may take 10-30 seconds)
data(HMPv13)

# Look at the object
HMPv13

# How many samples and how many OTUs (bacterial groups)?
nsamples(HMPv13)
ntaxa(HMPv13)

# What variables are in the sample metadata?
sample_variables(HMPv13)

meta <- as.data.frame(sample_data(HMPv13))
table(meta$body_site)