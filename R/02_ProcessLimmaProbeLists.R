library(org.Hs.eg.db)
library(annotate)
library(readr)
library(magrittr)
detach("package:dplyr", unload=TRUE)
library(dplyr)

results_nbm <- read_csv("results/meynert.neocortex.FALSE/limma/basal nucleus of meynert.allen_HBA.geneSummary.csv")
results_band <- read_csv("results/band.neocortex.FALSE/limma/diagonal band.allen_HBA.geneSummary.csv")

results_nbm %>% arrange(desc(tval)) %>% select(gene_symbol, tval)
results_band %<>% arrange(desc(tval)) %>% select(gene_symbol, tval)

write.table(results_nbm[1:2078,], file="results/NBM_top10.csv", sep="\t", row.names=F, col.names=F, quote=FALSE)
write.table(results_band[1:2078,], file="results/DB_top10.csv", sep="\t", row.names=F, col.names=F, quote=FALSE)

write.table(results_nbm[1:2078,], file="results/NBM_top10.csv", sep="\t", row.names=F, col.names=F, quote=FALSE)
write.table(results_band[1:2078,], file="results/DB_top10.csv", sep="\t", row.names=F, col.names=F, quote=FALSE)

write.table(results_nbm, file="results/NBM_ranked.csv", sep="\t", row.names=F, col.names=F, quote=FALSE)
write.table(results_band, file="results/DB_ranked.csv", sep="\t", row.names=F, col.names=F, quote=FALSE)

###

#adult_target <- "claustrum"
#fetal_target <- "claustrum"

#inclusionPattern <- "claustrum"

adult_target <- "meynert"
inclusionPattern <- "meynert"

inclusionPattern <- "band"
neoCortexOnly <- FALSE

baseFolder <- paste0(inclusionPattern, ".neocortex.", neoCortexOnly)

(adultFilename <- list.files(paste0("./results/",baseFolder,"/limma/"), pattern = paste0("*",adult_target,".allen_HBA.geneSummary.csv"), full.names = T)[1])

#load fetal and adult gene-wise data
adult <- read_csv(adultFilename, guess_max = 20000)

adultProbes <- read_csv(gsub("geneSummary.","", adultFilename))

#check counts
length(unique(adultProbes$gene_symbol))

adultProbes %<>% mutate(adj.p.value = p.adjust(p.value, method="fdr")) %>% arrange(p.value) %>% 
  mutate(sigUpRegulated = adj.p.value < 0.05 & t > 0) %>% print()

paste("Significantly enriched probes for adult", adult_target, ":", nrow(adultProbes %>% filter(sigUpRegulated)))

#get the probe level corrected max p-value
threshold <- max((adultProbes %>% filter(sigUpRegulated))$p.value)
adult %<>% arrange(-pValueWithDirection) %>% mutate(is_probe_corrected_significant = p.value <= threshold & direction > 0) %>% print()

sigAdultGenes <- unique(adult %>% filter(is_probe_corrected_significant) %>% .$gene_symbol)

cat("Significantly enriched genes for adult", adult_target, ":",  length(sigAdultGenes))

adult %<>% mutate(rank = rank(-pValueWithDirection))

geneToName <- adultProbes %>% select(gene_symbol, entrez_id) %>% distinct() %>% filter(!is.na(entrez_id))

#add gene names
geneToName %<>% mutate(name = unlist(lookUp(as.character(entrez_id), "org.Hs.eg", "GENENAME"))) %>% distinct() %>% print()
geneToName %<>% filter(!is.na(name))

adult <- left_join(adult, geneToName) %>% select(gene_symbol, name, entrez_id, everything())


#top 20
get_top_20 <- function(probe_table, gene_table) {
  #symbol, name, number of significant probes, p-value, fetal rank
  sigProbeCount <- probe_table %>% filter(sigUpRegulated) %>% group_by(gene_symbol) %>% summarize(sigProbes = dplyr::n())
  top20 <- gene_table %>% arrange(rank) %>% head(20)
  top20 %<>% inner_join(sigProbeCount)
  top20 %<>% mutate(p.value = signif(p.value, digits=3)) %>% 
    select(`Gene Symbol` = gene_symbol, `Name` = name, `Significant probes` = sigProbes, `p-value` = p.value) 
  top20
}

adult20 <- get_top_20(adultProbes, adult)

write_csv(adult, paste0(adultFilename, ".addedStats.csv") )
write_csv(adult20, paste0(adultFilename, ".top20.csv") )




