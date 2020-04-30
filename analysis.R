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
library(phyloseq)
library(DESeq2)

# Set working directory
setwd("~/compbio/hsm_community_seq/")
# Load data
hsm_otu_table <-
  read.table("data/otu_table.tsv", header = TRUE, sep = "\t")
hsm_16s_metadata <- read.csv("data/metadata.csv")
hsm_16s_taxonomy <-
  read.table("data/taxonomy.tsv", header = TRUE, sep = "\t")

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
hsm_otus_with_tax <- as.data.frame(vegan::decostand(hsm_otus_with_tax, method = "total", 2))
hsm_otus_with_tax$Feature.ID <- rownames(hsm_otus_with_tax)
hsm_otus_with_tax <- inner_join(hsm_otus_with_tax, hsm_16s_taxonomy) %>%
  pivot_longer(cols = 1:76, names_to = "sample_name", values_to = "abundance") %>%
  inner_join(hsm_16s_metadata)



#############
# Diversity #
#############
# Prepare Data

## Rarefy OTUs for diversity measurements
otu_names <- as.character(hsm_otu_table$OTU_ID)
hsm_otu_table <- data.frame(t(hsm_otu_table[, -1]))
colnames(hsm_otu_table) <- otu_names
hsm_otus_rarified <- data.frame(rrarefy(hsm_otu_table, 11095))

## Separate samples by soil type
s17_samples <-
  as.character(hsm_16s_metadata$sample_name[hsm_16s_metadata$soil == "s17" &
                                              hsm_16s_metadata$condition %in% condition_of_interest])

s17_rarified_otus <-
  filter(hsm_otus_rarified,
         rownames(hsm_otus_rarified) %in% s17_samples)
hsm_otus_rarified <- s17_rarified_otus

s17_unrarified_otus <-
  filter(hsm_otu_table, rownames(hsm_otu_table) %in% s17_samples)

# Alpha diversity 
## Calculate indices in vegan
alpha_diversity <- data.frame(
  "sample_name" = s17_samples,
  "inverse_simpsons" = diversity(hsm_otus_rarified, "invsimpson", MARGIN = 1),
  "shannon" = diversity(hsm_otus_rarified, MARGIN = 1),
  "chao1" = estimateR(hsm_otus_rarified)[2,]
) %>%
  inner_join(hsm_16s_metadata)

## Plot alpha diversity
alpha_diversity_plot <-
  gather(alpha_diversity, key = "index", value = "value", 2:4) %>%
  ggplot(aes(condition, value)) + 
  facet_wrap(. ~ index, scales = "free") +
  geom_boxplot() +
  theme_cowplot()

ggsave(plot = alpha_diversity_plot, "output/alpha_diversity.pdf")

## Alpha diversity stats
chao_pvals <- broom::tidy(TukeyHSD(aov(chao1 ~ condition, data = alpha_diversity)))
shannon_pvals <- broom::tidy(TukeyHSD(aov(shannon ~ condition, data = alpha_diversity)))
invsimp_pvals <- broom::tidy(TukeyHSD(aov(inverse_simpsons ~ condition, data = alpha_diversity)))
write_csv(data.frame("comparison" = chao_pvals$comparison,
                     "chao" = chao_pvals$adj.p.value,
                     "shannon" = shannon_pvals$adj.p.value,
                     "inv_simpsons" = invsimp_pvals$adj.p.value),
          path = "output/alpha_diversity_pvals.csv")

# Beta diversity
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

ggsave(plot = s17_pcoa_plot, "output/pcoa.pdf")

## Run PERMANOVA for each MDS plot
permanova_s17 <-
  adonis(s17_rarified_otus ~ condition, data = s17_pcoa_df)
capture.output(permanova_s17, file = "output/permanova_pcoa.txt")


## Phylum-Level Bar Chart
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

ggsave(plot = phylum_chart, "output/phylum_bar_chart.pdf")

##########
# DESeq2 #
##########

#Taxonomy Table
hsm_tax_ps <- as.matrix(hsm_16s_taxonomy[-1])
rownames(hsm_tax_ps) <- hsm_16s_taxonomy$Feature.ID
hsm_tax_ps <- tax_table(hsm_tax_ps)

#OTU Table
hsm_otus_ps <- otu_table(hsm_otu_table, taxa_are_rows = FALSE)

#Metadata
rownames(hsm_16s_metadata) <- hsm_16s_metadata$sample_name
hsm_metadata_ps <- sample_data(hsm_16s_metadata)

# Merge phyloseq object
hsm_ps <- phyloseq(tax_table=hsm_tax_ps,
                   sample_data=hsm_metadata_ps,
                   otu_table=hsm_otus_ps)
# Agglomerate at family level
hsm_ps <- tax_glom(hsm_ps, taxrank = "f")%>%
  phyloseq::subset_samples(soil == "s17")

# Run DESeq2
hsm_dds <- phyloseq_to_deseq2(hsm_ps, ~condition)
hsm_dds$condition <- relevel(hsm_dds$condition, ref = "Col0")
hsm_dds <- DESeq(hsm_dds)
deseq_results <- results(hsm_dds, name = "condition_HSM_vs_Col0", alpha = 0.05)

# DESeq Results formatting
hsm_da_families <- as.data.frame(deseq_results)
hsm_da_families$Feature.ID <- rownames(hsm_da_families)
hsm_da_families <- inner_join(hsm_da_families, hsm_16s_taxonomy)%>%
  transform(f = reorder(f, log2FoldChange))%>%
  filter(padj < 0.1)

hsm_da_families$plot_ranks <- factor(hsm_da_families$plot_ranks, levels = c(
  "Acidobacteria", "Actinobacteria", "Bacteroidetes","Cyanobacteria", "Firmicutes",
  "Gemmatimonadetes", "Planctomycetes", "Alphaproteobacteria", "Betaproteobacteria",
  "Deltaproteobacteria", "Gammaproteobacteria", "Verrucomicrobia", "Chloroflexi"))

# Plot DA taxa bubble plot
da_taxa_plot <- ggplot(hsm_da_families, aes(x = log2FoldChange, y = f, colour = plot_ranks, size = -log10(padj)))+
  geom_point()+
  geom_vline(xintercept = 0)+
  theme(axis.title.y=element_blank())

da_taxa_plot

ggsave(plot = da_taxa_plot, "output/differential_taxa.pdf")

# Write DESeq Results to csv
write_csv(hsm_da_families, "output/deseq_results.csv")
