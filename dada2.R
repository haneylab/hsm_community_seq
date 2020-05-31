library(dada2)

reads <- "/Users/haneylab/compbio/190219_CommunityAnalysis/trimmed_reads/"
list.files(reads)

f_seqs <- sort(list.files(reads, pattern = "_trimmed_R1.fastq", full.names = TRUE))
r_seqs <- sort(list.files(reads, pattern = "_trimmed_R2.fastq", full.names = TRUE))

sample.names <- sapply(strsplit(basename(f_seqs), "_"), `[`, 1)

#Quality profiles
plotQualityProfile(f_seqs[1:2])
plotQualityProfile(r_seqs[1:2])

f_filter <- file.path("filtered_reads", paste0(sample.names, "F_filtered.fastq.gz"))
r_filter <- file.path("filtered_reads", paste0(sample.names, "R_filtered.fastq.gz"))

names(f_filter) <- sample.names
names(r_filter) <- sample.names

out <- filterAndTrim(f_seqs, f_filter, r_seqs, r_filter, truncLen = c(200, 180),
                     maxN = 0, maxEE = c(4,2), truncQ = 2, rm.phix = TRUE, compress = TRUE,
                     multithread = TRUE)

#Learn error rates
errF <- learnErrors(f_filter, multithread = TRUE)
errR <- learnErrors(r_filter, multithread = TRUE)

plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)

#Dereplication/merging

mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sample in sample.names) {
  cat("Processing:", sample, "\n")
    derepF <- derepFastq(f_filter[[sample]])
    ddf <- dada(derepF, err = errF, multithread = TRUE)
    derepR <- derepFastq(r_filter[[sample]])
    ddr <- dada(derepR, err = errR, multithread = TRUE)
    merger <- mergePairs(ddf, derepF, ddr, derepR)
    mergers[[sample]] <- merger
}
 rm(derepF); rm(derepR)



st <- makeSequenceTable(mergers)
seqTab <- removeBimeraDenovo(st, method =  "consensus", multithread = TRUE)
saveRDS(seqTab, "output/sequence_table.rds")

#Assign taxonomy
tax <- assignTaxonomy(seqTab, "tax_db/silva_nr_v132_train_set.fa", multithread=TRUE)
saveRDS(tax, "output/taxonomy.rds")

getN <- function(x) sum(getUniques(x))
track <- cbind(out,sapply(mergers, getN), rowSums(seqTab))
track

