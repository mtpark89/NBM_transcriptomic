library(AnnotationDbi)
library(annotate)
library(org.Hs.eg.db)
library(readxl)
detach("package:dplyr", unload=TRUE)
library(dplyr)
library(magrittr)
library(limma)
library(readr)

neoCortexOnly <- TRUE
sourceExpression <- "allen_HBA"

if (sourceExpression == "allen_HBA") {
  folderPattern <- "normalized_microarray.*"
  sampleFilename <- "SampleAnnot.csv"
  probeFilename <- "Probes.csv"
  expressionFilename <- "MicroarrayExpression.csv"
} else if (sourceExpression == "allen_human_fetal_brain") {
  folderPattern <- "lmd_matrix_.*"
  sampleFilename <- "columns_metadata.csv"
  probeFilename <- "rows_metadata.csv"
  expressionFilename <- "expression_matrix.csv"
}

strip_left_right <- function(structure_name) {
  tokens <- trimws(unlist(strsplit(structure_name, ",")))
  tokens <- tokens[tokens != "left"]
  tokens <- tokens[tokens != "right"]
  cleaned_name <- paste(tokens, collapse = ", ")
  cleaned_name
}

#Replace with Gabe's re-registered files
sampleFilename <- "SampleAnnot_nlin.csv"

#fetal brain has more/better gene symbol mappings
probeInfo <- read_csv(paste0("./data/raw/allen_human_fetal_brain/lmd_matrix_12566/rows_metadata.csv"))
probeInfo %<>% rename(probe_id = probeset_id, probe_name = probeset_name)

allsampleAnnot = NULL
allExpression = NULL

for (donorFolder in list.files(paste0("./data/raw/",sourceExpression,"/"), pattern = folderPattern)) {
  sampleAnnot <- read_csv(paste0("./data/raw/",sourceExpression,"/",donorFolder,"/",sampleFilename))
  
  expressionMatrix <- read_csv(paste0("./data/raw/",sourceExpression,"/",donorFolder, "/",expressionFilename), col_names=F) 
    
  expressionMatrix %<>% rename(probe_id = X1)
    
  sampleAnnot %<>% rowwise() %>% mutate(structure_name_left_right_stripped = strip_left_right(structure_name))
  sampleAnnot %<>% mutate(donorID = donorFolder)
  
  if (sourceExpression == "allen_HBA") {
    sampleAnnot %<>% mutate(uniqueID = paste("ID", structure_id, slab_num, well_id, polygon_id, donorID, sep=".")) %>% select(uniqueID, everything())
    } else if (sourceExpression == "allen_human_fetal_brain") {
      sampleAnnot %<>% mutate(uniqueID = paste("ID", structure_id, well_id, donorID, sep=".")) %>% select(uniqueID, everything())
    }

    colnames(expressionMatrix) <- c("probe_id", sampleAnnot$uniqueID)
    
    expressionMatrix <- inner_join(probeInfo %>% select(probe_id, probe_name), expressionMatrix) %>% select(-probe_id)
    
    #bind cols of expression matrix
    allExpression <- bind_cols(allExpression, expressionMatrix)
    
    #bind rows of sample annot
    allsampleAnnot <- bind_rows(allsampleAnnot, sampleAnnot)
  
}

colnames(allExpression)[1]<- "probe_name"

##################################################################

sampleAnnot <- allsampleAnnot

#Only cortical & left hemisphere
sampleAnnot %<>% filter(grepl("gyri|pole|gyrus|cuneus|planum|operculum|lobule", structure_name)) %>% filter(!grepl("dentate gyrus|right", structure_name))

###Mni points, CIVET resampled on 09c sym template
mni <- read.csv('data/mni_sym09c_mid_surface_rsl_left_81920_points.csv', sep=" ")

distance_vertex <- matrix(nrow=nrow(sampleAnnot), ncol=2)

library(proxy)

for (i in 1:nrow(sampleAnnot)) {
	
	distances <- dist(sampleAnnot[i,15:17], mni[,2:4]) ###distance calculations
	
	distance_vertex[i, 1] <- min(distances)
	distance_vertex[i, 2] <- which.min(distances)
}
colnames(distance_vertex) <- c("Distance", "Vertex")
distance_vertex %<>% as_tibble()

sampleAnnot <- bind_cols(distance_vertex, sampleAnnot)

###Quality control: Exclude with distance greater than 10mm--wrong matches, on right side.
sampleAnnot %<>% filter(Distance < 10)

###Deal with duplicate vertices

length(sampleAnnot$Vertex) #1317
length(unique(sampleAnnot$Vertex)) #1236

sampleAnnot %<>% group_by(Vertex) %>% slice(which.min(Distance))

###Ontology mapping
ontology <- read_csv('data/raw/allen_HBA/normalized_microarray_donor9861/Ontology.csv')
ontology %<>% rowwise() %>% mutate(Level5=as.double(strsplit(structure_id_path, split="/")[[1]][6]))

level5 <- ontology %>% select(Level5) %>% left_join(ontology, by=c("Level5"="id")) %>% select(Level5_name="name", Level5_colour="color_hex_triplet")

ontology <- bind_cols(ontology, level5)

write.table(ontology, file="data/Ontology.csv", sep="\t", row.names=F)

sampleAnnot <- left_join(sampleAnnot, ontology, by=c("structure_id"="id"))

sampleAnnot_MNI <- sampleAnnot

###CIVET surface check for visual quality control

CIVET_surface <- as.data.frame(1:40962)
colnames(CIVET_surface)<- "Vertex"

donors_vertices <- as.data.frame(sampleAnnot$Vertex)
colnames(donors_vertices)<- "Vertex"
donors_vertices$Sample <- 1

CIVET_surface <- merge(CIVET_surface, donors_vertices, select="Vertex", all.x=TRUE)
CIVET_surface[is.na(CIVET_surface)] <- 0

write.table(CIVET_surface, file="results/imaging_align/CIVET_surface_AHBA_samples.txt", row.names=F, col.names=F, sep=" ")

#Set in same order
expressionMatrix <- allExpression[, c("probe_name", sampleAnnot_MNI$uniqueID)]

print(paste("Number of unique regions:", length(unique(sampleAnnot$structure_name_left_right_stripped))))
print(paste("Number of samples:", nrow(sampleAnnot)))
print(paste("Number of donors:", length(unique(sampleAnnot$donorID))))

expressionMatrix <- as.data.frame(expressionMatrix)
rownames(expressionMatrix) <- expressionMatrix$probe_name
expressionMatrix$probe_name <- NULL

#create probe mapping file here
#qc_table <- read_xlsx("./data/gene_symbol_annotations/Miller et al. doi.org_10.1186_1471-2164-15-154 12864_2013_7016_MOESM8_ESM.xlsx",skip=1)
qc_table <- read_xlsx("./data/gene_symbol_annotations/Miller_formatted.xlsx")
qc_table %<>% select(probe_name = "Probe", gene_symbol="Gene", is_qc_pass = `Pass?`, R=`R`, pval=`P-value`, qval=`Q-value`, AvgInt=`Average Intensity`)

#use gene symbol to get NCBI ID, then get updated symbol
symbolToID <- probeInfo %>% select(probe_id, probe_name, gene_symbol, entrez_id) %>% distinct()
symbolToID %<>% filter(!is.na(entrez_id)) #Remove probes without Entrez ID
symbolToID %<>% mutate(new_symbol = getSYMBOL(as.character(entrez_id), data='org.Hs.eg')) 

qc_table <- left_join(qc_table, symbolToID)
qc_table %<>% mutate(legacySymbol = gene_symbol) #save old symbol
qc_table %<>% mutate(gene_symbol = if_else(is.na(new_symbol), gene_symbol, new_symbol)) %>% dplyr::select(gene_symbol, entrez_id,  everything(), -new_symbol)

#remove those not mapping
qc_table <- qc_table %>% filter(!grepl("A_", gene_symbol)) %>% filter(!grepl("CUST_", gene_symbol)) 
qc_table %<>% filter(gene_symbol != "na" & entrez_id!="na")

#Select QC passing probes. For gene with multiple probes, select one with highest correlation=15120 genes.
qc_pass <- qc_table %>% filter(is_qc_pass==1) %>% group_by(gene_symbol) %>% slice(which.max(R))

#Selecting one probe per gene
expressionMatrix_MNI <- expressionMatrix[qc_pass$probe_name,]
identical(rownames(expressionMatrix_MNI), qc_pass$probe_name)

print(paste("Probes after filtering for _ and numeric geneSymbols", nrow(qc_pass)))
print(paste("Gene count",length(unique(qc_pass$gene_symbol))))

#Write out expression and sample files, handle rownames
#dir.create(paste0("./data/processed/", inclusionPattern, ".neocortex.", neoCortexOnly), showWarnings = F)
write_tsv(expressionMatrix_MNI %>% mutate(probe_name = rownames(expressionMatrix_MNI)) %>% select(probe_name, everything()), 
          paste0("./data/processed/Cortical_expressionMatrix_qc_15120genes.tsv"))
write_tsv(sampleAnnot_MNI, paste0("./data/processed/Cortical_expressionMatrix_annotations.tsv"))

#sampleAnnot %>% filter(grepl(inclusionPattern, structure_name_left_right_stripped)) %>% group_by(structure_name_left_right_stripped) %>% summarize(n=dplyr::n())
#as.data.frame(sampleAnnot %>% filter(grepl(inclusionPattern, structure_name_left_right_stripped)) %>% select(structure_name))

###Correlating imaging to genetic
identical(colnames(expressionMatrix_MNI), sampleAnnot_MNI$uniqueID)
identical(rownames(expressionMatrix_MNI), qc_pass$probe_name)

rownames(expressionMatrix_MNI) <- qc_pass$gene_symbol

cortical <- read_tsv(file = "data/Left_cortical_cors.csv")
sampleAnnot_MNI <- left_join(sampleAnnot_MNI, cortical, by="Vertex")

library(lmerTest)

ahba_mat <- matrix(data=NA, nrow=nrow(expressionMatrix_MNI), ncol=3)

for (i in 1:nrow(expressionMatrix_MNI)){
	model <- summary(lmer(sampleAnnot_MNI$NBM_cor ~ t(expressionMatrix_MNI[i,]) + (1|sampleAnnot_MNI$donorID), verbose=0))

	ahba_mat[i,1] <- rownames(expressionMatrix_MNI)[i]
	ahba_mat[i,2:3] <- model$coefficients[2,4:5] ###T-statistic and p-value
}

ahba_mat2 <- matrix(data=NA, nrow=nrow(expressionMatrix_MNI), ncol=3)

for (i in 1:nrow(expressionMatrix_MNI)){
	model <- summary(lmer(sampleAnnot_MNI$Ch123_cor ~ t(expressionMatrix_MNI[i,]) + (1|sampleAnnot_MNI$donorID), verbose=0))
	ahba_mat2[i,1] <- rownames(expressionMatrix_MNI)[i]
	ahba_mat2[i,2:3] <- model$coefficients[2,4:5] ###T-statistic and p-value
}

colnames(ahba_mat) <- c("Gene", "tstat", "pval")
colnames(ahba_mat2) <- c("Gene", "tstat", "pval")

write.table(ahba_mat, file="results/AHBA_NBM-cortical_backup_20210409.csv", sep=",", row.names=F)
write.table(ahba_mat2, file="results/AHBA_DB-cortical_backup_20210409.csv", sep=",", row.names=F)

#########################################################################

ahba_mat %<>% as_tibble()
ahba_mat2 %<>% as_tibble()

ahba_mat %<>% mutate_at(c("tstat", "pval"), as.numeric)
ahba_mat$qval <- p.adjust(ahba_mat$pval, method="fdr")
ahba_mat %<>% mutate(logpvalsign= -log10(pval) * sign(tstat))
ahba_mat %<>% arrange(desc(logpvalsign))

ahba_mat2 %<>% mutate_at(c("tstat", "pval"), as.numeric)
ahba_mat2$qval <- p.adjust(ahba_mat2$pval, method="fdr")
ahba_mat2 %<>% mutate(logpvalsign= -log10(pval) * sign(tstat))
ahba_mat2 %<>% arrange(desc(logpvalsign))

#Write out ranked list, FDR 10, 5, 1% corrected genelists

write.table(ahba_mat %>% filter(tstat > 0 & qval < 0.10) %>% select(Gene), file="AHBA_NBM-cortical_FDR10.txt", sep="\t", row.names=F, col.names=F, quote=FALSE)
write.table(ahba_mat %>% filter(tstat > 0 & qval < 0.05) %>% select(Gene), file="AHBA_NBM-cortical_FDR5.txt", sep="\t", row.names=F, col.names=F, quote=FALSE)
write.table(ahba_mat %>% filter(tstat > 0 & qval < 0.01) %>% select(Gene), file="AHBA_NBM-cortical_FDR1.txt", sep="\t", row.names=F, col.names=F, quote=FALSE)

write.table(ahba_mat2 %>% filter(tstat > 0 & qval < 0.10) %>% select(Gene), file="AHBA_DB-cortical_FDR10.txt", sep="\t", row.names=F, col.names=F, quote=FALSE)
write.table(ahba_mat2 %>% filter(tstat > 0 & qval < 0.05) %>% select(Gene), file="AHBA_DB-cortical_FDR5.txt", sep="\t", row.names=F, col.names=F, quote=FALSE)
write.table(ahba_mat2 %>% filter(tstat > 0 & qval < 0.01) %>% select(Gene), file="AHBA_DB-cortical_FDR1.txt", sep="\t", row.names=F, col.names=F, quote=FALSE)

write.table(ahba_mat %>% select(Gene, logpvalsign), file="AHBA_NBM-cortical_scores.txt", sep="\t", row.names=F, col.names=F, quote=FALSE)
write.table(ahba_mat2 %>% select(Gene, logpvalsign), file="AHBA_DB-cortical_scores.txt", sep="\t", row.names=F, col.names=F, quote=FALSE)

##################################################################################
###Looking at cholinergic genes

cholgenes <- read.csv('Genes_cholinergic.csv', sep="\t")

expressionMatrix_cholinergic <- expressionMatrix_MNI[cholgenes$gene_symbol,]

ahba_mat <- left_join(ahba_mat, cholgenes, by=c("Gene"="gene_symbol"))
ahba_mat2 <- left_join(ahba_mat2, cholgenes, by=c("Gene"="gene_symbol"))

ahba_mat$Group <- replace_na(ahba_mat$Group, "All")
ahba_mat2$Group <- replace_na(ahba_mat2$Group, "All")

ahba_mat %<>% mutate(size=ifelse(Group=="Cholinergic", 3, 1))
ahba_mat2 %<>% mutate(size=ifelse(Group=="Cholinergic", 3, 1))

plot_nbm <-ggplot(ahba_mat %>% arrange(Group), aes(x=tstat, y=-log10(pval), color=Group, size=size)) + geom_point() + theme_Publication() + scale_colour_manual(values=c("#4682B4", "red3")) + scale_size_continuous(range=c(1,3)) + geom_hline(yintercept=-log10(ahba_mat %>% filter(qval < 0.01) %>% select(pval) %>% max()), alpha=0.5, linetype="dashed")+ geom_label_repel(data=subset(ahba_mat, Group=="Cholinergic" & qval < 0.01), aes(label=Gene), size=5) + theme(legend.position="none") + xlab("t-statistic") + ylab("-log10(p value)")

plot_db<-ggplot(ahba_mat2 %>% arrange(Group), aes(x=tstat, y=-log10(pval), color=Group, size=size)) + geom_point() + theme_Publication() + scale_colour_manual(values=c("#4682B4", "red3")) + scale_size_continuous(range=c(1,3)) + geom_hline(yintercept=-log10(ahba_mat2 %>% filter(qval < 0.01) %>% select(pval) %>% max()), alpha=0.5, linetype="dashed")+ geom_label_repel(data=subset(ahba_mat2, Group=="Cholinergic" & qval < 0.01), aes(label=Gene), size=5) + theme(legend.position="none") + xlab("t-statistic") + ylab("-log10(p value)") 
 
for_plot<- bind_cols(NBM_cor=sampleAnnot_MNI$NBM_cor, Ch123_cor=sampleAnnot_MNI$Ch123_cor, t(expressionMatrix_MNI["CHRNA3",]), t(expressionMatrix_MNI["CHRNA2",]))

plot1 <- ggplot(for_plot, aes(x=Ch123_cor, y=CHRNA2)) + theme_Publication()+  geom_point() + geom_smooth(method=lm) + xlab("Ch1-3-cortical correlations") + ylab("CHRNA2 expression")
plot2<- ggplot(for_plot, aes(x=NBM_cor, y=CHRNA3)) + theme_Publication()+  geom_point() + geom_smooth(method=lm) + xlab("NBM-cortical correlations") + ylab("CHRNA3 expression")

plot_cors<- ggarrange(plot_db, plot_nbm, plot1, plot2, ncol=2, nrow=2)

ggsave(filename="Figure_cors.png", plot=plot_cors, scale=1, width=12, height=10, units="in", dpi=300)


###
for_reg<- bind_cols(NBM_cor=sampleAnnot_MNI$NBM_cor, Ch123_cor=sampleAnnot_MNI$Ch123_cor, t(expressionMatrix_MNI[cholgenes$gene_symbol,]))

expressionMatrix <- expressionMatrix[qc_pass$probe_name,]

#Modified to run differential expression across cortical samples now
regionByGene <- NULL

for (targetRegion in sort(unique(sampleAnnot$uniqueID))) {
  print(targetRegion)
  sampleAnnot %<>% mutate(isTarget = (targetRegion == uniqueID))
  sampleAnnot %>% filter(isTarget) %>% select(donorID, everything())
  
  designMatrix <- model.matrix(~ sampleAnnot$isTarget + sampleAnnot$donorID)
  
  fit <- lmFit(expressionMatrix,designMatrix)
  fit <- eBayes(fit)
  
  limmaResults <- inner_join(as_tibble(fit$p.value[,"sampleAnnot$isTargetTRUE", drop=F], rownames="probe_name"), as_tibble(fit$t[,"sampleAnnot$isTargetTRUE", drop=F], rownames="probe_name"), by="probe_name")
  limmaResults %<>% rename( p.value=`sampleAnnot$isTargetTRUE.x`, t = `sampleAnnot$isTargetTRUE.y`)
  
  #Only probes passing QC
  limmaResults <- inner_join(limmaResults, qc_pass, by= "probe_name")
  #may have slight bias for longer genes with more probes
  gene_summary <- limmaResults %>% group_by(gene_symbol) %>% arrange(p.value) %>% 
    summarize(p.value = first(p.value), direction=sign(first(t)), tval=first(t))
  #convert to ranks
  gene_summary %<>% mutate(pValueWithDirection = direction * (nrow(gene_summary) - rank(p.value)))

  #write out for target region
  #if(grepl(inclusionPattern, targetRegion)) {
  #  print("Writing out")
  #  dir.create(paste0("./results/",inclusionPattern, ".neocortex.", neoCortexOnly), showWarnings = F)
  #  dir.create(paste0("./results/",inclusionPattern, ".neocortex.", neoCortexOnly, "/limma/"), showWarnings = F)
  #  write_csv(limmaResults, paste0("./results/",inclusionPattern, ".neocortex.", neoCortexOnly, "/limma/", targetRegion,".",sourceExpression,".csv"))
  #  write_csv(gene_summary, paste0("./results/",inclusionPattern, ".neocortex.", neoCortexOnly, "/limma/", targetRegion,".",sourceExpression,".geneSummary.csv"))
  #}
  
  gene_summary %<>% select(gene_symbol, !! targetRegion := pValueWithDirection)
  
  #join to a region by gene matrix
  if (is.null(regionByGene)) { 
    regionByGene <- gene_summary 
  } else { 
    regionByGene <- inner_join(regionByGene, gene_summary, by= "gene_symbol") 
  }
}

write_tsv(regionByGene, paste0("./data/processed/Cortical_1236samples_vs_15120genes.tsv"))
