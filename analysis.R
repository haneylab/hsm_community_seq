# Andrew Wilson
# 14 November 2019
# Song, Y. et al

#########
# Setup #
#########

# Load packages
library(tidyverse)
library(vegan)
library(phyloseq)
library(DESeq2)
library(cowplot)

# Set working directory
setwd("./compbio/hsm_community_seq/")
# Load data
hsm_otu_table <-
  read.table("otu_table.tsv", header = TRUE, sep = "\t")
hsm_16s_metadata <- read.csv("metadata.csv")
hsm_16s_taxonomy <-
  read.table("taxonomy.tsv", header = TRUE, sep = "\t")
condition_of_interest <- c("BS", "HSM", "Col0", "W.Col0")

## Fix Taxonomy file
hsm_16s_taxonomy$Taxon <-
  str_replace_all(hsm_16s_taxonomy$Taxon,  "[a-z]__", "")
hsm_16s_taxonomy <-
  separate(hsm_16s_taxonomy,
           col = Taxon,
           into = c("k", "p", "c", "o", "f", "g", "s")) %>%
  replace_na(list(
    k = "",
    p = "",
    c = "",
    o = "",
    f = "",
    g = "",
    s = ""
  ))

#############
# Diversity #
#############

# Rarefy OTUs for diversity measurements
otu_names <- as.character(hsm_otu_table$OTU_ID)
hsm_otu_table <- data.frame(t(hsm_otu_table[, -1]))
colnames(hsm_otu_table) <- otu_names
hsm_otus_rarified <- data.frame(rrarefy(hsm_otu_table, 11095))

# Separate samples by soil type
s17_samples <-
  as.character(hsm_16s_metadata$sample_name[hsm_16s_metadata$soil == "s17" &
                                              hsm_16s_metadata$condition %in% condition_of_interest])

s17_rarified_otus <-
  filter(hsm_otus_rarified,
         rownames(hsm_otus_rarified) %in% s17_samples)
hsm_otus_rarified <- s17_rarified_otus

s17_unrarified_otus <-
  filter(hsm_otu_table, rownames(hsm_otu_table) %in% s17_samples)

# Alpha Diversity measurements
alpha_diversity <- data.frame(
  "sample_name" = s17_samples,
  "inverse_simpsons" = diversity(hsm_otus_rarified, "invsimpson", MARGIN = 1),
  "shannon" = diversity(hsm_otus_rarified, MARGIN = 1),
  "chao1" = estimateR(hsm_otus_rarified)[2,]
) %>%
  gather(key = "metric", value = "index",-1) %>%
  inner_join(hsm_16s_metadata)


alpha_diversity_plot <-
  ggplot(alpha_diversity, aes(condition, index)) +
  facet_wrap(. ~ metric, scales = "free") +
  geom_boxplot() +
  theme_cowplot()

alpha_diversity_plot

# MDS of Bray-Curtis dissimilarity
## Run MDS on S17 samples
s17_mds <- data.frame(metaMDS(s17_rarified_otus)$points)
s17_mds$sample_name <- s17_samples
s17_mds <- inner_join(hsm_16s_metadata, s17_mds)
## Plot S17 MDS
s17_mds_plot <-
  ggplot(s17_mds, aes(MDS1, MDS2, colour = condition)) +
  geom_point(size = 2.5) +
  coord_fixed() +
  theme_cowplot()
s17_mds_plot

## Run PERMANOVA for each MDS plot
permanova_s17 <-
  adonis(s17_rarified_otus ~ condition, data = s17_mds)
capture.output(permanova_s17, file = "permanova.txt")


#############
# Abundance #
#############

# Import into Phyloseq
ps_taxonomy <- tax_table(as.matrix(hsm_16s_taxonomy))
ps_metadata <- sample_data(hsm_16s_metadata)
ps_otu_table <- otu_table(hsm_otu_table, taxa_are_rows = FALSE)
