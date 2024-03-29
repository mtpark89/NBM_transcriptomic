library(AnnotationDbi)
library(annotate)
library(org.Hs.eg.db)
library(readxl)
detach("package:dplyr", unload=TRUE)
library(dplyr)
library(magrittr)
library(limma)
library(readr)

inclusionPattern <- "meynert"

#inclusionPattern <- "band"

neoCortexOnly <- FALSE

#for fetal data
stripLayers <- TRUE

#uncomment line below to set adult or fetal expression source
sourceExpression <- "allen_HBA"
#sourceExpression <- "allen_human_fetal_brain"

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

#this just checks to see how many samples match the inclusion pattern
for (donorFolder in list.files(paste0("./data/raw/",sourceExpression,"/"), pattern = folderPattern)) {
  sampleAnnot <- read_csv(paste0("./data/raw/",sourceExpression,"/",donorFolder,"/",sampleFilename))
  
  #check if donor has target region
  if (nrow(filter(sampleAnnot, grepl(inclusionPattern, structure_name))) > 0) {
    print(donorFolder)
    print(nrow(filter(sampleAnnot, grepl(inclusionPattern, structure_name))))
    print(sort(filter(sampleAnnot, grepl(inclusionPattern, structure_name))$structure_name))
    print(nrow(sampleAnnot))
  }
}

#fetal brain has more/better gene symbol mappings
probeInfo <- read_csv(paste0("./data/raw/allen_human_fetal_brain/lmd_matrix_12566/rows_metadata.csv"))
probeInfo %<>% rename(probe_id = probeset_id, probe_name = probeset_name)

#Replace with Gabe's re-registered files
sampleFilename <- "SampleAnnot_nlin.csv"

allsampleAnnot = NULL
allExpression = NULL

for (donorFolder in list.files(paste0("./data/raw/",sourceExpression,"/"), pattern = folderPattern)) {
  sampleAnnot <- read_csv(paste0("./data/raw/",sourceExpression,"/",donorFolder,"/",sampleFilename))
  
  #check if donor has target region
  if (nrow(filter(sampleAnnot, grepl(inclusionPattern, structure_name))) > 0) {
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
}

#Numbers don't line up, because of duplicate probe_name...1, etc columns
colnames(allExpression)[1]<- "probe_name"

#DB specific grouping
if (inclusionPattern == "band") {
  allsampleAnnot %<>% mutate(structure_name_left_right_stripped = if_else(structure_name_left_right_stripped == "nucleus of the diagonal band, vertical division", "diagonal band", structure_name_left_right_stripped))
  allsampleAnnot %<>% mutate(structure_name_left_right_stripped = if_else(structure_name_left_right_stripped == "nucleus of the diagonal band, horizontal division", "diagonal band", structure_name_left_right_stripped))
  allsampleAnnot %>% filter(grepl(inclusionPattern, structure_name_left_right_stripped)) %>% group_by(structure_name_left_right_stripped) %>% summarize(n=dplyr::n())
  allsampleAnnot %>% filter(grepl(inclusionPattern, structure_name_left_right_stripped))
}

################################

sampleAnnot <- allsampleAnnot

#Remove diagonal band sample that has wrong placement

filter(sampleAnnot, grepl("band", structure_name))

as.data.frame(filter(sampleAnnot, grepl("band", structure_name))) %>% filter(mni_nlin_z > 36)

#unique ID  ID.4310.6.156906860.126532035.normalized_microarray_donor15496

sampleAnnot %<>% subset(uniqueID!="ID.4310.6.156906860.126532035.normalized_microarray_donor15496")

################################

expressionMatrix <- allExpression[, c("probe_name", sampleAnnot$uniqueID)]

print(paste("Number of unique regions:", length(unique(sampleAnnot$structure_name_left_right_stripped))))
print(paste("Number of samples:", nrow(sampleAnnot)))
print(paste("Number of donors:", length(unique(sampleAnnot$donorID))))
print("Region breakdown:")
sampleAnnot %>% mutate(Region = if_else(grepl(inclusionPattern, structure_name_left_right_stripped), structure_name_left_right_stripped, "remaining structures")) %>% group_by(Region) %>% summarize(n=dplyr::n())

expressionMatrix <- as.data.frame(expressionMatrix)
rownames(expressionMatrix) <- expressionMatrix$probe_name
expressionMatrix$probe_name <- NULL

#create probe mapping file here
#qc_table <- read_xlsx("./data/gene_symbol_annotations/Miller et al. doi.org_10.1186_1471-2164-15-154 12864_2013_7016_MOESM8_ESM.xlsx",skip=1)
#qc_table %<>% select(probe_name = "...1", gene_symbol = "...2", is_qc_pass = `Pass?`)
qc_table <- read_xlsx("./data/gene_symbol_annotations/Miller_formatted.xlsx")

qc_table %<>% select(probe_name = "Probe", gene_symbol="Gene", is_qc_pass = `Pass?`, R=`R`, pval=`P-value`, qval=`Q-value`, AvgInt=`Average Intensity`)

#fix numeric gene symbols from excel date conversions
#qc_table <- inner_join( qc_table, probeInfo %>% select(probe_name, for_excel_problem = gene_symbol))
#qc_table %<>% mutate(gene_symbol = if_else(!is.na(as.numeric(gene_symbol)), for_excel_problem, gene_symbol)) %>% select(-for_excel_problem)

#use gene symbol to get NCBI ID, then get updated symbol
symbolToID <- probeInfo %>% select(probe_id, probe_name, gene_symbol, entrez_id) %>% distinct()
symbolToID %<>% filter(!is.na(entrez_id)) 
symbolToID %<>% mutate(new_symbol = getSYMBOL(as.character(entrez_id), data='org.Hs.eg')) 

qc_table <- left_join(qc_table, symbolToID)
qc_table %<>% mutate(legacySymbol = gene_symbol) #save old symbol
qc_table %<>% mutate(gene_symbol = if_else(is.na(new_symbol), gene_symbol, new_symbol)) %>% dplyr::select(gene_symbol, entrez_id,  everything(), -new_symbol)

#remove those not mapping
qc_table <- qc_table %>% filter(!grepl("A_", gene_symbol)) %>% filter(!grepl("CUST_", gene_symbol)) 
qc_table %<>% filter(gene_symbol != "na" & entrez_id!="na")

#Select QC passing probes. For gene with multiple probes, select one with highest correlation=15120 genes.
qc_pass <- qc_table %>% filter(is_qc_pass==1) %>% group_by(gene_symbol) %>% slice(which.max(R))

expressionMatrix <- expressionMatrix[qc_pass$probe_name,]

print(paste("Probes after filtering for _ and numeric geneSymbols", nrow(qc_pass)))
print(paste("Gene count",length(unique(qc_pass$gene_symbol))))

#handle rownames
dir.create(paste0("./data/processed/", inclusionPattern, ".neocortex.", neoCortexOnly), showWarnings = F)
write_tsv(expressionMatrix %>% mutate(probe_name = rownames(expressionMatrix)) %>% select(probe_name, everything()), 
          paste0("./data/processed/",inclusionPattern, ".neocortex.", neoCortexOnly, "/", sourceExpression, "_sampleMatrix_qcNames.tsv"))
write_tsv(sampleAnnot, paste0("./data/processed/",inclusionPattern, ".neocortex.", neoCortexOnly, "/", sourceExpression, "_sampleMatrix_colAnnotations_qcNames.tsv"))

sampleAnnot %>% filter(grepl(inclusionPattern, structure_name_left_right_stripped)) %>% group_by(structure_name_left_right_stripped) %>% summarize(n=dplyr::n())
as.data.frame(sampleAnnot %>% filter(grepl(inclusionPattern, structure_name_left_right_stripped)) %>% select(structure_name))

#####

regionByGene <- NULL
for (targetRegion in sort(unique(sampleAnnot$structure_name_left_right_stripped))) {
  print(targetRegion)
  sampleAnnot %<>% mutate(isTarget = (targetRegion == structure_name_left_right_stripped))
  sampleAnnot %>% filter(isTarget) %>% select(donorID, everything())
  
  designMatrix <- model.matrix(~ sampleAnnot$isTarget + sampleAnnot$donorID)
  
  fit <- lmFit(expressionMatrix,designMatrix)
  fit <- eBayes(fit)
  
  limmaResults <- inner_join(as_tibble(fit$p.value[,"sampleAnnot$isTargetTRUE", drop=F], rownames="probe_name"), as_tibble(fit$t[,"sampleAnnot$isTargetTRUE", drop=F], rownames="probe_name"), by="probe_name")
  limmaResults %<>% rename( p.value=`sampleAnnot$isTargetTRUE.x`, t = `sampleAnnot$isTargetTRUE.y`)
  
  limmaResults <- inner_join(limmaResults, qc_table, by= "probe_name")
  #may have slight bias for longer genes with more probes
  gene_summary <- limmaResults %>% group_by(gene_symbol) %>% arrange(p.value) %>% 
    summarize(p.value = first(p.value), direction=sign(first(t)), tval=first(t))
  #convert to ranks
  gene_summary %<>% mutate(pValueWithDirection = direction * (nrow(gene_summary) - rank(p.value)))

  #write out for target region
  if(grepl(inclusionPattern, targetRegion)) {
    print("Writing out")
    dir.create(paste0("./results/",inclusionPattern, ".neocortex.", neoCortexOnly), showWarnings = F)
    dir.create(paste0("./results/",inclusionPattern, ".neocortex.", neoCortexOnly, "/limma/"), showWarnings = F)
    write_csv(limmaResults, paste0("./results/",inclusionPattern, ".neocortex.", neoCortexOnly, "/limma/", targetRegion,".",sourceExpression,".csv"))
    write_csv(gene_summary, paste0("./results/",inclusionPattern, ".neocortex.", neoCortexOnly, "/limma/", targetRegion,".",sourceExpression,".geneSummary.csv"))
  }
  
  gene_summary %<>% select(gene_symbol, !! targetRegion := pValueWithDirection)
  
  
  #join to a region by gene matrix
  if (is.null(regionByGene)) { 
    regionByGene <- gene_summary 
  } else { 
    regionByGene <- inner_join(regionByGene, gene_summary, by= "gene_symbol") 
  }
}

write_tsv(regionByGene, paste0("./data/processed/",inclusionPattern, ".neocortex.", neoCortexOnly, "/", sourceExpression, "_brainarea_vs_genes_exp_qcNames.tsv"))

################################################################################

results_nbm <- read_csv("results/meynert.neocortex.FALSE/limma/basal nucleus of meynert.allen_HBA.geneSummary.csv")
results_band <- read_csv("results/band.neocortex.FALSE/limma/diagonal band.allen_HBA.geneSummary.csv")

results_nbm %<>% arrange(desc(tval)) %>% mutate(qval = p.adjust(p.value, method="fdr"))
results_band %<>% arrange(desc(tval)) %>% mutate(qval = p.adjust(p.value, method="fdr"))

write.table(results_nbm %>% filter(qval < 0.05 & tval > 0) %>% select(gene_symbol), file="results/NBM_FDR5.csv", sep="\t", row.names=F, col.names=F, quote=FALSE)
write.table(results_band %>% filter(qval < 0.05 & tval > 0) %>% select(gene_symbol), file="results/DB_FDR5.csv", sep="\t", row.names=F, col.names=F, quote=FALSE)

write.table(results_nbm %>% filter(qval < 0.10 & tval > 0) %>% select(gene_symbol), file="results/NBM_FDR10.csv", sep="\t", row.names=F, col.names=F, quote=FALSE)
write.table(results_band %>% filter(qval < 0.10 & tval > 0) %>% select(gene_symbol), file="results/DB_FDR10.csv", sep="\t", row.names=F, col.names=F, quote=FALSE)

write.table(results_nbm %>% filter(qval < 0.01 & tval > 0) %>% select(gene_symbol), file="results/NBM_FDR1.csv", sep="\t", row.names=F, col.names=F, quote=FALSE)
write.table(results_band %>% filter(qval < 0.01 & tval > 0) %>% select(gene_symbol), file="results/DB_FDR1.csv", sep="\t", row.names=F, col.names=F, quote=FALSE)

write.table(results_nbm[1:1512,] %>% select(gene_symbol), file="results/NBM_top10.csv", sep="\t", row.names=F, col.names=F, quote=FALSE)
write.table(results_band[1:1512,] %>% select(gene_symbol), file="results/DB_top10.csv", sep="\t", row.names=F, col.names=F, quote=FALSE)

write.table(results_nbm[1:756,] %>% select(gene_symbol), file="results/NBM_top5.csv", sep="\t", row.names=F, col.names=F, quote=FALSE)
write.table(results_band[1:756,] %>% select(gene_symbol), file="results/DB_top5.csv", sep="\t", row.names=F, col.names=F, quote=FALSE)

write.table(results_nbm %>% select(gene_symbol, tval), file="results/NBM_ranked.csv", sep="\t", row.names=F, col.names=F, quote=FALSE)
write.table(results_band %>% select(gene_symbol, tval), file="results/DB_ranked.csv", sep="\t", row.names=F, col.names=F, quote=FALSE)

