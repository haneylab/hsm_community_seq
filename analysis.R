# Andrew Wilson
# 14 November 2019
# Song, Y. et al

#########
# Setup #
#########

# Load packages
library(tidyverse)
library(vegan)
library(ape)
library(cowplot)

# Set working directory
setwd("./compbio/hsm_community_seq/")
# Load data
hsm_otu_table <-
  read.table("otu_table.tsv", header = TRUE, sep = "\t")
hsm_16s_metadata <- read.csv("metadata.csv")
hsm_16s_taxonomy <-
  read.table("taxonomy.tsv", header = TRUE, sep = "\t")

condition_of_interest <- c("BS", "HSM", "Col0")

## Generate vector of ranks for phylum-level comparisons
taxa <- as.character(hsm_16s_taxonomy$Taxon)
plot_ranks <- sapply(taxa, function(x) {
  taxon_list <- unlist(strsplit(x, "; "))
  if (length(taxon_list) == 1 | taxon_list[2] == "p__") {
    return(paste("Unclassified", str_replace(taxon_list[1], "k__", "")))
  } else if (taxon_list[2] == "p__Proteobacteria" &
             length(taxon_list) >= 3) {
    return(str_replace(taxon_list[3], "c__", ""))
  } else {
    return(str_replace(taxon_list[2], "p__", ""))
  }
})

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
hsm_16s_taxonomy$plot_ranks <- plot_ranks
endosymbionts <- filter(hsm_16s_taxonomy, c == "Chloroplast" | f == "mitochondria")

hsm_otu_table <- filter(hsm_otu_table, !(OTU_ID %in% endosymbionts$Feature.ID))
rownames(hsm_otu_table) <- hsm_otu_table$OTU_ID
## Create master dataframe with all data
hsm_otus_with_tax <- hsm_otu_table[,-1]
hsm_otus_with_tax <-
  as.data.frame(vegan::decostand(hsm_otus_with_tax, method = "total", 2))
hsm_otus_with_tax$Feature.ID <- rownames(hsm_otus_with_tax)
hsm_otus_with_tax <-
  inner_join(hsm_otus_with_tax, hsm_16s_taxonomy) %>%
  pivot_longer(cols = 1:76,
               names_to = "sample_name",
               values_to = "abundance") %>%
  inner_join(hsm_16s_metadata)



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

ggsave(plot = alpha_diversity_plot, "alpha_diversity.pdf")

# NMDS/PCoA of Bray-Curtis dissimilarity
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
ggsave(plot = s17_mds_plot, "nmds.pdf")

## Run PCoA
s17_dist <- vegdist(s17_rarified_otus, "bray")
s17_pcoa <- pcoa(s17_dist)
s17_pcoa_df <- as.data.frame(s17_pcoa$vectors[, 1:2])
s17_pcoa_df$sample_name <- s17_samples
s17_pcoa_df <- inner_join(hsm_16s_metadata, s17_pcoa_df)
axis_names <- s17_pcoa$values$Relative_eig[1:2]

## Plot PCoA
s17_pcoa_plot <-
  ggplot(s17_pcoa_df, aes(Axis.1, Axis.2, colour = condition)) +
  geom_point(size = 2.5) +
  coord_fixed() +
  xlab(paste0("PC1", " (", round(axis_names[1] * 100, digits = 2), "%)")) +
  ylab(paste0("PC2", " (", round(axis_names[2] * 100, digits = 2), "%)")) +
  theme_cowplot()

ggsave(plot = s17_pcoa_plot, "pcoa.pdf")

## Run PERMANOVA for each MDS plot
permanova_s17 <-
  adonis(s17_rarified_otus ~ condition, data = s17_pcoa_df)
capture.output(permanova_s17, file = "permanova_pcoa.txt")


# Phylum-Level Bar Chart
phylum_chart_df <- hsm_otus_with_tax %>%
  group_by(sample_name, plot_ranks, soil, condition)%>%
  summarize(abundance = sum(abundance))%>%
  filter(condition %in% condition_of_interest & soil == "s17")

phylum_chart_df$taxon <- sapply(1:987, function(x){
  if(phylum_chart_df[x,5] >= 0.02){
    unlist(phylum_chart_df[x, 2])
  } else {
    return("Other")
  }
})

phylum_chart_df$taxon <- factor(phylum_chart_df$taxon,levels = c(
  "Acidobacteria", "Actinobacteria", "Bacteroidetes","Cyanobacteria", "Firmicutes",
  "Gemmatimonadetes", "Planctomycetes", "Alphaproteobacteria", "Betaproteobacteria",
  "Deltaproteobacteria", "Gammaproteobacteria", "Verrucomicrobia", "Other"))

phylum_chart_df <- phylum_chart_df %>%
  group_by(sample_name, taxon, soil, condition)%>%
  summarize(abundance = sum(abundance))

phylum_chart <- ggplot(phylum_chart_df, aes(x = sample_name, y = abundance, fill = taxon))+
  facet_grid(.~condition, scales = "free_x", space = "free")+
  geom_bar(stat = "identity")+
  scale_y_continuous(name = "Relative Abundance")+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggsave(plot = phylum_chart, "phylum_bar_chart.pdf")

## Code for multiple t tests goes here

#############
# Abundance #
#############


plot_family <-
  function(family,
           conditions = condition_of_interest,
           soil = "s17") {
    subset_df <-
      filter(hsm_otus_with_tax,
             f == family &
               condition %in% conditions & soil == "s17") %>%
      group_by(sample_name, condition) %>%
      summarize(ag_abundance = sum(abundance))
    capture.output(
      t.test(ag_abundance ~ condition, data = subset_df),
      file = paste0("output/", family, ".ttest.txt")
    )
    write_csv(subset_df, path = paste0("output/", family, "_abundances.csv"))
  }
