library(AnnotationDbi)
library(annotate)
library(org.Hs.eg.db)
library(readxl)
detach("package:dplyr", unload=TRUE)
library(dplyr)
library(magrittr)
library(limma)
library(readr)

results_nbm <- read_csv("results/meynert.neocortex.FALSE/limma/basal nucleus of meynert.allen_HBA.geneSummary.csv")
results_band <- read_csv("results/band.neocortex.FALSE/limma/diagonal band.allen_HBA.geneSummary.csv")

results_nbm %<>% arrange(desc(tval)) %>% mutate(qval = p.adjust(p.value, method="fdr"))
results_band %<>% arrange(desc(tval)) %>% mutate(qval = p.adjust(p.value, method="fdr"))

db_scz <- read_tsv("results/ToppGene/DB/ToppGene_DB_SCZ_overlap.csv")

db_scz <- left_join(db_scz, results_band, by=c("Gene Symbol" = "gene_symbol")) %>% arrange(desc(tval))
db_scz %<>% mutate(logpval= -log10(p.value)*direction) %>% arrange(desc(logpval))

write_tsv(db_scz, "results/ToppGene/DB/ToppGene_DB_SCZ_overlap_withstats.csv")

nbm_scz <- read_tsv("results/ToppGene/NBM/ToppGene_NBM_SCZ_overlap.csv")
nbm_scz <- left_join(nbm_scz, results_nbm, by=c("Gene Symbol" = "gene_symbol")) %>% arrange(desc(tval))
nbm_scz %<>% mutate(logpval= -log10(p.value)*direction) %>% arrange(desc(logpval))

write_tsv(nbm_scz, "results/ToppGene/NBM/ToppGene_NBM_SCZ_overlap_withstats.csv")

intersect(nbm_scz$`Gene Symbol`, cholgenes$gene_symbol)
intersect(db_scz$`Gene Symbol`, cholgenes$gene_symbol)
########################################################################

db_topp_disease <- read_tsv("results/ToppGene/DB/ToppGene_output_DB.txt") %>% filter(Category=="Disease") %>% head(20)
db_topp_drug <- read_tsv("results/ToppGene/DB/ToppGene_output_DB.txt") %>% filter(Category=="Drug") %>% head(20)

nbm_topp_disease <- read_tsv("results/ToppGene/NBM/ToppGene_output_NBM.txt") %>% filter(Category=="Disease") %>% head(20)
nbm_topp_drug <- read_tsv("results/ToppGene/NBM/ToppGene_output_NBM.txt") %>% filter(Category=="Drug") %>% head(20)

db_topp_disease %<>% mutate(colour=ifelse(Name=="Schizophrenia", "1", "0"))
nbm_topp_disease %<>% mutate(colour=ifelse(Name=="Schizophrenia", "1", "0"))

db$Structure <- "DB"
nbm$Structure <- "NBM"

combined <- rbind.fill(db, nbm)

plot1<- ggplot(db_topp_disease %>% head(14), aes(x=reorder(Name, -log10(`q-value FDR B&H`)), y=-log10(`q-value FDR B&H`), fill=colour)) + geom_bar(stat='identity', alpha=0.6) + coord_flip() + theme_Publication() + xlab("") + ylab("") + theme(legend.position = "none", axis.text.y=element_text(size=10))+ scale_fill_lancet() + scale_y_continuous(position="right")

plot2<- ggplot(nbm_topp_disease, aes(x=reorder(Name, -log10(`q-value FDR B&H`)), y=-log10(`q-value FDR B&H`), fill=colour)) + geom_bar(stat='identity', alpha=0.6) + coord_flip() + theme_Publication() + xlab("") + ylab("") + theme(legend.position = "none", axis.text.y=element_text(size=10))+ scale_fill_lancet() + scale_y_continuous(position="right")

fig_topp <- ggarrange(plot1, plot2, ncol=1, nrow=2)
ggsave(filename="Fig_Topp.png", plot=fig_topp, scale=1, width=5, height=8, units="in", dpi=300)

write_tsv(db_topp_drug, "results/ToppGene/DB/ToppGene_DB_drug_top20.csv")
