require(ggplot2)
require(reshape2)
require(dplyr)
require(MCL)
blast_results = read.csv2("~/GitHub/yabacus/blastn.tab", header = F, stringsAsFactors = F, sep="\t", col.names = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"))
# remove self hits
#blast_results <- blast_results[blast_results$qseqid != blast_results$sseqid, ]
blast_results$bitscore <- as.numeric(blast_results$bitscore)
blast_results$pident <- as.numeric(blast_results$pident)
blast_results$length <- as.numeric(blast_results$length)

names = c("pident", "length", "mismatch", "gapopen", "qstart", "evalue", "bitscore")
for (i in 2:10){
  blast_results[, paste0("k_", i)] <-  paste(kmeans(blast_results[, names], i)$cluster)
}
ggplot(blast_results, aes(x=bitscore, y=length, color=k_5)) + geom_point()

head(blast_results)
# filter those
max_len <- max(blast_results$length)
blast_results <- blast_results %>% 
#  filter(pident > 99) %>%
#  filter(length > max_len*.95) %>%
  group_by(qseqid) %>%
  mutate(norm_score = bitscore/max(bitscore)) %>%
  top_n(n=1, wt=norm_score) %>%
  as.data.frame()
str(blast_results)
blast_nxn <- dcast(blast_results, formula = qseqid ~ sseqid, value.var = "norm_score")
rownames(blast_nxn) <- blast_nxn$qseqid
str(blast_nxn)
blast_nxn$qseqid <- NULL
blast_nxn[is.na(blast_nxn)] <- 0
blast_matrix <- as.matrix(blast_nxn)
diag(blast_matrix) <- 0
mcl(blast_nxn, addLoops = T)
blast_dist <-dist(blast_nxn, diag = T)
blast_tall_dist <- melt(as.matrix(blast_dist), value.name = "val", varnames = c("qseqid", "sseqid"))
mcl(as.matrix(blast_dist), addLoops = T)
str(blast_tall_dist)
ggplot(blast_tall_dist, aes(x=sseqid, y=qseqid, alpha=val) ) + geom_point()
##########################################33

# graph results
