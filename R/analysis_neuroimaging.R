library(readr)
library(dplyr)
library(magrittr)

summary(lm(gf$Ch123 ~ Age + Sex + CAST_total + Smoker + DX, data=gf))
summary(lm(gf$L_Ch4 ~ Age + Sex + CAST_total + Smoker + DX, data=gf))
summary(lm(gf$R_Ch4 ~ Age + Sex + CAST_total + Smoker + DX, data=gf))

summary(lm(gf$AD_Ch123 ~ Age + Sex + CAST_total + Smoker + DX, data=gf))
summary(lm(gf$AD_L_Ch4 ~ Age + Sex + CAST_total + Smoker + DX, data=gf))
summary(lm(gf$AD_R_Ch4 ~ Age + Sex + CAST_total + Smoker + DX, data=gf))

#######
summary(lm(gf$Ch123 ~ Age + Sex + CAST_total + Smoker + ICV + DX, data=gf))
summary(lm(gf$L_Ch4 ~ Age + Sex + CAST_total + Smoker + ICV + DX, data=gf))
summary(lm(gf$R_Ch4 ~ Age + Sex + CAST_total + Smoker + ICV + DX, data=gf))

summary(lm(gf$AD_Ch123 ~ Age + Sex + CAST_total + Smoker + ICV + DX, data=gf))
summary(lm(gf$AD_L_Ch4 ~ Age + Sex + CAST_total + Smoker + ICV + DX, data=gf))
summary(lm(gf$AD_R_Ch4 ~ Age + Sex + CAST_total + Smoker + ICV + DX, data=gf))

#######Cognition post-hoc

#########################FINAL COGNITIVE SCORES TEST
gf <- gf %>% mutate(DSST=(SymbolSubWritten + SymbolSubOral)/2)

summary(lm(DSST ~ DX, data=gf))
summary(lm(TrailMakingTime ~ DX, data=gf))

tapply(gf$DSST, gf$DX, mean, na.rm=TRUE)
tapply(gf$TrailMakingTime, gf$DX, mean, na.rm=TRUE)

tapply(gf$DSST, gf$DX, sd, na.rm=TRUE)
tapply(gf$TrailMakingTime, gf$DX, sd, na.rm=TRUE)

###FEP
cor.test(gf_FEP$L_Ch4, gf_FEP$DSST)
cor.test(gf_FEP$R_Ch4, gf_FEP$DSST)
cor.test(gf_FEP$Ch123, gf_FEP$DSST)

cor.test(gf_FEP$L_Ch4, gf_FEP$TrailMakingTime)
cor.test(gf_FEP$R_Ch4, gf_FEP$TrailMakingTime)
cor.test(gf_FEP$Ch123, gf_FEP$TrailMakingTime)

cor.test(gf_FEP$AD_L_Ch4, gf_FEP$DSST)
cor.test(gf_FEP$AD_R_Ch4, gf_FEP$DSST)
cor.test(gf_FEP$AD_Ch123, gf_FEP$DSST)

cor.test(gf_FEP$AD_L_Ch4, gf_FEP$TrailMakingTime)
cor.test(gf_FEP$AD_R_Ch4, gf_FEP$TrailMakingTime)
cor.test(gf_FEP$AD_Ch123, gf_FEP$TrailMakingTime)

#HC
cor.test(gf_HC$L_Ch4, gf_HC$DSST)
cor.test(gf_HC$R_Ch4, gf_HC$DSST)
cor.test(gf_HC$Ch123, gf_HC$DSST)

cor.test(gf_HC$L_Ch4, gf_HC$TrailMakingTime)
cor.test(gf_HC$R_Ch4, gf_HC$TrailMakingTime)
cor.test(gf_HC$Ch123, gf_HC$TrailMakingTime)

cor.test(gf_HC$AD_L_Ch4, gf_HC$DSST)
cor.test(gf_HC$AD_R_Ch4, gf_HC$DSST)
cor.test(gf_HC$AD_Ch123, gf_HC$DSST)

cor.test(gf_HC$AD_L_Ch4, gf_HC$TrailMakingTime)
cor.test(gf_HC$AD_R_Ch4, gf_HC$TrailMakingTime)
cor.test(gf_HC$AD_Ch123, gf_HC$TrailMakingTime)
########################################

#Variance testing

library(cars)

leveneTest(Ch123 ~ DX, gf)
leveneTest(L_Ch4 ~ DX, gf)
leveneTest(R_Ch4 ~ DX, gf)

leveneTest(AD_Ch123 ~ DX, gf)
leveneTest(AD_L_Ch4 ~ DX, gf)
leveneTest(AD_R_Ch4 ~ DX, gf)

#Plots

gf$Ch123_scaled <- scale(gf$Ch123)
gf$MD_Ch123_scaled <- scale(gf$MD_Ch123)
gf$AD_Ch123_scaled <- scale(gf$AD_Ch123)

melted <- melt(gf, id.vars=c("DX"), measure.vars=c("Ch123_scaled", "MD_Ch123_scaled", "AD_Ch123_scaled"))

melted$variable<- as.factor(gsub("MD_Ch123_scaled", "MD", melted$variable))
melted$variable<- as.factor(gsub("AD_Ch123_scaled", "AD", melted$variable))
melted$variable<- as.factor(gsub("Ch123_scaled", "qT1", melted$variable))
melted$variable2 = factor(melted$variable, levels=c('qT1', 'MD', 'AD'))

plot1<- ggplot(melted, aes(x=DX, y=value, group=DX, colour=DX)) + theme_Publication() + facet_wrap(~variable2) + geom_point(size=2) + geom_boxplot() + geom_jitter(position=position_jitter(0.2)) + theme(axis.text.x=element_blank(), axis.ticks=element_blank(), axis.title.x=element_blank()) + ylab("Z-score")+ scale_colour_brewer(palette="Set1", direction=-1)

ggsave(filename="Figure_boxplots.png", plot=plot1, scale=1, width=8, height=6, units="in", dpi=300)

#MRS analysis (choline)

cor.test(gf_HC$L_Ch4, gf_HC$Choline)
cor.test(gf_HC$R_Ch4, gf_HC$Choline)
cor.test(gf_HC$Ch123, gf_HC$Choline)

cor.test(gf_FEP$L_Ch4, gf_FEP$Choline)
cor.test(gf_FEP$R_Ch4, gf_FEP$Choline)
cor.test(gf_FEP$Ch123, gf_FEP$Choline)

###Correlations in healthy controls taking into account covariates
summary(lm(Choline ~ Age + Sex + CAST_total + Ch123, data=gf_HC))
summary(lm(Choline ~ Age + Sex + CAST_total + L_Ch4, data=gf_HC))
summary(lm(Choline ~ Age + Sex + CAST_total + R_Ch4, data=gf_HC))

###Interaction models
summary(lm(Choline ~ Age+ Sex + Smoker + Ch123*DX, data=gf))
summary(lm(Choline ~ Age+ Sex + Smoker + L_Ch4*DX, data=gf))
summary(lm(Choline ~ Age+ Sex + Smoker + R_Ch4*DX, data=gf))

summary(lm(Choline ~ Age+ Sex + CAST_total + Smoker + Ch123*DX, data=gf))
summary(lm(Choline ~ Age+ Sex + CAST_total + Smoker + L_Ch4*DX, data=gf))
summary(lm(Choline ~ Age+ Sex + CAST_total + Smoker + R_Ch4*DX, data=gf))

summary(lm(Choline ~ Age+ Sex + Ch123*DX, data=gf))
summary(lm(Choline ~ Age+ Sex + L_Ch4*DX, data=gf))
summary(lm(Choline ~ Age+ Sex + R_Ch4*DX, data=gf))

gf$L_Ch4_scaled <- scale(gf$L_Ch4)

melted <- melt(gf, id.vars=c("DX", "Choline"), measure.vars=c("Ch123_scaled", "L_Ch4_scaled"))
melted$variable<- as.factor(gsub("Ch123_scaled", "Ch123", melted$variable))
melted$variable<- as.factor(gsub("L_Ch4_scaled", "Left NBM", melted$variable))

plot2 <- ggplot(melted, aes(x=value, y=Choline, group=DX, colour=DX)) + geom_point(size=2) + geom_smooth(method=lm, size=2) + facet_wrap(~variable) + theme_Publication() + xlab("qT1 (Z-scores)") + scale_colour_brewer(palette="Set1", direction=-1)

ggsave(filename="Figure_MRS.png", plot=plot2, scale=1, width=8, height=6, units="in", dpi=300)

###Relation to clinical scores

cor.test(gf_FEP$L_Ch4, gf_FEP$totalN)
cor.test(gf_FEP$Ch123, gf_FEP$totalN)

cor.test(gf_FEP$L_Ch4, gf_FEP$YMRS_total)
cor.test(gf_FEP$Ch123, gf_FEP$YMRS_total)

cor.test(gf_FEP$AD_L_Ch4, gf_FEP$totalN)
cor.test(gf_FEP$AD_Ch123, gf_FEP$totalN)

cor.test(gf_FEP$AD_L_Ch4, gf_FEP$YMRS_total)
cor.test(gf_FEP$AD_Ch123, gf_FEP$YMRS_total)

gf_FEP_YMRS <- subset(gf_FEP, YMRS_total!="NA" & CAST_total!="NA")

model <- lm(YMRS_total ~ CAST_total + Smoker, data=gf_FEP_YMRS)
gf_FEP_YMRS$YMRS_resid<- model$resid+ mean(gf_FEP$YMRS_total, na.rm=TRUE)

gf_FEP$L_Ch4_scaled <- scale(gf_FEP$L_Ch4)
gf_FEP$AD_L_Ch4_scaled <- scale(gf_FEP$AD_L_Ch4)

plot1<- ggplot(gf_FEP_YMRS, aes(x=L_Ch4_scaled, y=YMRS_resid)) + geom_point(size=2) + geom_smooth(method=lm, size=2) + theme_Publication() + ylab("Total YMRS (residual)") + xlab("Left NBM qT1 (Z-score)")

plot2<- ggplot(gf_FEP, aes(x=AD_L_Ch4_scaled, y=totalN)) + geom_point(size=2) + geom_smooth(method=lm, size=2) + theme_Publication() + ylab("Negative symptoms") + xlab("Left NBM AD (Z-score)")

plot_scores<- ggarrange(plot1, plot2, ncol=2)

ggsave(filename="Figure_clinical.png", plot=plot_scores, scale=1, width=8, height=6, units="in", dpi=300)

###NBM correlations

get_cors <- function(data1, data2) {
	test <- cor.test(data1, data2)
	return(cbind(test$estimate, test$p.value))
}

get_vertexcors <- function(column, vdata) {
	output <- NULL
	#Assuming rows are subjects and columns are vertices
	for (i in 1:ncol(vdata)){
		output <- rbind(output, get_cors(column, vdata[,i]))
	}
	output <- cbind(output, p.adjust(output[,2], method="fdr"))	
	return(output)
}

left_hipp_labels <- read.csv('/projects/Resources_new/Hippocampus_model_labels/new_lefthc_labels.txt', header=FALSE)
right_hipp_labels <- read.csv('/projects/Resources_new/Hippocampus_model_labels/new_righthc_labels.txt', header=FALSE)

write.table(get_vertexcors(gf$L_Ch4, left_striatum_T1), file="cor_NBM/left_striatum_cor.txt", col.names=F, row.names=F)
write.table(get_vertexcors(gf$L_Ch4, left_thalamus_T1), file="cor_NBM/left_thal_cor.txt", col.names=F, row.names=F)
write.table(get_vertexcors(gf$L_Ch4, left_gp_T1), file="cor_NBM/left_gp_cor.txt", col.names=F, row.names=F)
write.table(cbind(get_vertexcors(gf$L_Ch4, left_hipp_T1), left_hipp_labels), file="cor_NBM/left_hipp_cor.txt", col.names=F, row.names=F)
write.table(get_vertexcors(gf$L_Ch4, left_amyg_T1), file="cor_NBM/left_amyg_cor.txt", col.names=F, row.names=F)

write.table(get_vertexcors(gf$R_Ch4, right_striatum_T1), file="cor_NBM/right_striatum_cor.txt", col.names=F, row.names=F)
write.table(get_vertexcors(gf$R_Ch4, right_thalamus_T1), file="cor_NBM/right_thal_cor.txt", col.names=F, row.names=F)
write.table(get_vertexcors(gf$R_Ch4, right_gp_T1), file="cor_NBM/right_gp_cor.txt", col.names=F, row.names=F)
write.table(cbind(get_vertexcors(gf$R_Ch4, right_hipp_T1), right_hipp_labels), file="cor_NBM/right_hipp_cor.txt", col.names=F, row.names=F)
write.table(get_vertexcors(gf$R_Ch4, right_amyg_T1), file="cor_NBM/right_amyg_cor.txt", col.names=F, row.names=F)

###Ch123 cor
write.table(get_vertexcors(gf$Ch123, left_striatum_T1), file="cor_Ch123/left_striatum_cor_Ch123.txt", col.names=F, row.names=F)
write.table(get_vertexcors(gf$Ch123, left_thalamus_T1), file="cor_Ch123/left_thal_cor_Ch123.txt", col.names=F, row.names=F)
write.table(get_vertexcors(gf$Ch123, left_gp_T1), file="cor_Ch123/left_gp_cor_Ch123.txt", col.names=F, row.names=F)
write.table(cbind(get_vertexcors(gf$Ch123, left_hipp_T1), left_hipp_labels), file="cor_Ch123/left_hipp_cor_Ch123.txt", col.names=F, row.names=F)
write.table(get_vertexcors(gf$Ch123, left_amyg_T1), file="cor_Ch123/left_amyg_cor_Ch123.txt", col.names=F, row.names=F)

write.table(get_vertexcors(gf$Ch123, right_striatum_T1), file="cor_Ch123/right_striatum_cor_Ch123.txt", col.names=F, row.names=F)
write.table(get_vertexcors(gf$Ch123, right_thalamus_T1), file="cor_Ch123/right_thal_cor_Ch123.txt", col.names=F, row.names=F)
write.table(get_vertexcors(gf$Ch123, right_gp_T1), file="cor_Ch123/right_gp_cor_Ch123.txt", col.names=F, row.names=F)
write.table(cbind(get_vertexcors(gf$Ch123, right_hipp_T1), right_hipp_labels), file="cor_Ch123/right_hipp_cor_Ch123.txt", col.names=F, row.names=F)
write.table(get_vertexcors(gf$Ch123, right_amyg_T1), file="cor_Ch123/right_amyg_cor_Ch123.txt", col.names=F, row.names=F)

gf_subcortical <- gf

median(cors_subcortical)

gf <- subset(gf_backup, QC_CIVET=="Pass")

gf$left_T1_depth3 <- paste("/projects/TOPSY_final_20200504/5_CIVET_t1map/",gf$Subject,"/depth_50_left_native_T1map_fwhm20_rsl.txt", sep="")
gf$right_T1_depth3 <- paste("/projects/TOPSY_final_20200504/5_CIVET_t1map/",gf$Subject,"/depth_50_right_native_T1map_fwhm20_rsl.txt", sep="")

left_T1_depth3<- t(vertexTable(gf$left_T1_depth3))
right_T1_depth3 <- t(vertexTable(gf$right_T1_depth3))

write.table(get_vertexcors(gf$L_Ch4, left_T1_depth3), file="left_cortical_cor.txt", col.names=F, row.names=F)
write.table(get_vertexcors(gf$R_Ch4, right_T1_depth3), file="right_cortical_cor.txt", col.names=F, row.names=F)

write.table(get_vertexcors(gf$Ch123, left_T1_depth3), file="left_cortical_cor_Ch123.txt", col.names=F, row.names=F)
write.table(get_vertexcors(gf$Ch123, right_T1_depth3), file="right_cortical_cor_Ch123.txt", col.names=F, row.names=F)

cors_cortical<- rbind(
t(cor(gf$L_Ch4, left_T1_depth3)),
t(cor(gf$R_Ch4, right_T1_depth3)),
t(cor(gf$Ch123, left_T1_depth3)),
t(cor(gf$Ch123, right_T1_depth3))
)

median(cors_cortical) #-0.1187

gf_cortical <- gf

###Gene expression to cortical

###AHBA

ahba <- read.csv('../../Allen_to_CIVET/Analysis_mapping/Donors_cortex_20200823.csv', sep="\t")

###Left cortical correlations

cor_BF <- cbind(as.data.frame(t(cor(gf$L_Ch4, left_T1_depth3))), as.data.frame(t(cor(gf$Ch123, left_T1_depth3))))
colnames(cor_BF) <- c("NBM_cor", "Ch123_cor")

cor_BF$Vertex <- 1:nrow(cor_BF) #Need vertex column to merge

write.table(cor_BF, file="Left_cortical_cors.csv", sep="\t", row.names=F)

ahba_1 <- merge(left_cors, ahba, select="Vertex", all.y=TRUE)

nrow(ahba_1) #1236 vertices to samples

###NBM-cortical qT1 correlations to correlate with gene expression
library(lmerTest)

ahba_mat <- matrix(data=NA, nrow=15126, ncol=3)

for (i in 25:ncol(ahba_1)){
	model <- summary(lmer(NBM_cor ~ ahba_1[,i] + (1|Donor), data=ahba_1, verbose=0))

	j=i-24
	ahba_mat[j,1] <- colnames(ahba_1)[i]
	ahba_mat[j,2:3] <- model$coefficients[2,4:5] ###T-statistic
	
}

ahba_mat <- as.data.frame(ahba_mat)
colnames(ahba_mat) <- c("Gene", "tstat", "pval")
ahba_mat$tstat <- as.numeric(ahba_mat$tstat)
ahba_mat$pval <- as.numeric(ahba_mat$pval)
ahba_mat$qval <- p.adjust(ahba_mat$pval, method="fdr")

ahba_mat <- ahba_mat[order(ahba_mat$tstat, decreasing=TRUE),]

write.table(ahba_mat, file="AHBA_mat_backup_20210211.csv", sep=",", row.names=F)

write.table(ahba_mat, file="AHBA_leftNBM_lmerTstats.txt", sep="\t", row.names=F, col.names=F)
write.table(ahba_mat$Gene, file="AHBA_leftNBM_lmerTstats_genes.txt", sep="\t", row.names=F, col.names=F)
write.table(ahba_mat$Gene[1:1512], file="AHBA_leftNBM_lmerTstats_genestop10.txt", sep="\t", row.names=F, col.names=F)
write.table(ahba_mat$Gene[13614:15126], file="AHBA_leftNBM_lmerTstats_genesbottom10.txt", sep="\t", row.names=F, col.names=F)

cholgenes <- read.csv('Genes_cholinergic.csv', sep="\t")
ahba_mat <- merge(ahba_mat, cholgenes, select="Gene", all.x=TRUE)

ahba_mat$Group <- replace_na(ahba_mat$Group, "Rest")
ahba_mat$Group <- ifelse(ahba_mat$Group=="Rest", 'All', 'Cholinergic')

ahba_mat$size <- ifelse(ahba_mat$Group=="Cholinergic", '3', '1')
ahba_mat$size <- as.numeric(ahba_mat$size)
ahba_mat$size<- replace_na(ahba_mat$size, 1)

ggplot(ahba_mat %>% arrange(Group), aes(x=tstat, y=-log10(pval), color=Group)) + geom_point() + theme_Publication() + scale_colour_manual(values=c("#4682B4", "red3"))

plot_vol<- ggplot(ahba_mat %>% arrange(Group), aes(x=tstat, y=-log10(pval), color=Group, size=size)) + geom_point() + theme_Publication() + scale_colour_manual(values=c("#4682B4", "red3")) + scale_size_continuous(range=c(1,3)) + geom_label_repel(data=subset(ahba_mat, Group=="Cholinergic" & qval < 0.05), aes(label=Gene), size=5) + theme(legend.position="none") + xlab("t-statistic") + ylab("-log10(p value)")

ggsave(filename="Figure_volcano.png", plot=plot_vol, scale=1, width=8, height=6, units="in", dpi=300)

###Layer for arranging points, not used
layer1 <- ahba_mat[ahba_mat$Group=="Rest",]
layer2 <- ahba_mat[ahba_mat$Group=="Cholinergic",]

ggplot() + theme_Publication() +  
    geom_point(
        data=layer1,
        aes(x=tstat, y=-log10(pval)), 
        colour="#4682B4", 
        size=layer1$size) +
    geom_point(
        data=layer2, 
        aes(x=tstat, y=-log10(pval)), 
        colour="red3", 
        size=layer2$size) + geom_text()


###

ahba_1_cors <- cor(ahba_1$NBM_cor, ahba_1[,25:ncol(ahba_1)])
ahba_1_cors<- as.data.frame(t(ahba_1_cors))

ahba_1_cors$Gene <- row.names(ahba_1_cors)
ahba_1_cors <- ahba_1_cors[order(ahba_1_cors$V1, decreasing=TRUE),]

write.table(ahba_1_cors$Gene, file="AHBA_leftNBM_cor_genes.txt", row.names=F, col.names=F)
write.table(ahba_1_cors$Gene[1:1512], file="AHBA_leftNBM_cor_genes_top10.txt", row.names=F, col.names=F)
write.table(ahba_1_cors$Gene[1:756], file="AHBA_leftNBM_cor_genes_top5.txt", row.names=F, col.names=F)

write.table(cbind(ahba_1_cors$Gene, ahba_1_cors$V1), file="AHBA_leftNBM_allcors.txt", row.names=F, col.names=F, quote=FALSE, sep="\t")


###histograms

melted <- melt(gf_master, id.vars=c("Subject"), measure.vars=c("L_Ch4", "R_Ch4", "Ch123"))

melted$variable<- gsub("L_", "Left ", melted$variable)
melted$variable<- gsub("R_", "Right ", melted$variable)

fig_histo<- ggplot(melted, aes(x=value, group=variable, colour=variable, fill=variable)) + theme_Publication() + geom_density(alpha=0.5) + ylab("Density") + xlab("qT1")

ggsave(filename="histograms.png", plot=fig_histo, scale=1, width=8, height=6, units="in", dpi=300)

###Demographics table

table(gf$DX, gf$Sex)
chisq.test(table(gf$DX, gf$Sex))

tapply(gf$Age, gf$DX, summary)
tapply(gf$Age, gf$DX, sd)
summary(lm(Age ~ DX, data=gf))

tapply(gf$PANSSTOTAL, gf$DX, mean, na.rm=TRUE)
tapply(gf$totalP, gf$DX, mean, na.rm=TRUE)
tapply(gf$totalN, gf$DX, mean, na.rm=TRUE)

tapply(gf$PANSSTOTAL, gf$DX, sd, na.rm=TRUE)
tapply(gf$totalP, gf$DX, sd, na.rm=TRUE)
tapply(gf$totalN, gf$DX, sd, na.rm=TRUE)

tapply(gf$CDS_total, gf$DX, mean, na.rm=TRUE)
tapply(gf$CDS_total, gf$DX, sd, na.rm=TRUE)

tapply(gf$SOFAS, gf$DX, mean, na.rm=TRUE)
tapply(gf$SOFAS, gf$DX, sd, na.rm=TRUE)

demo <- read_tsv("Demographics_DDD_2021.csv")

gf_demo <- left_join(as_tibble(gf), demo, by="ID")

tapply(gf_demo$DDD_LifeTime, gf_demo$DX, mean, na.rm=TRUE)
tapply(gf_demo$DDD_LifeTime, gf_demo$DX, sd, na.rm=TRUE)

tapply(gf_demo$DUP_Weeks, gf_demo$DX, mean, na.rm=TRUE)
tapply(gf_demo$DUP_Weeks, gf_demo$DX, sd, na.rm=TRUE)

table(gf_demo$DX, gf_demo$ConsensusDx)

summary(gf_cortical$DX)
summary(gf_subcortical$DX)

new_master <- as_tibble(read.spss("TOPSY_Database_MASTER_20210511.sav"))

test <- subset(gf_master, Ch123!="NA")

summary(as.factor(test$QC_sag))
summary(as.factor(test$QC_axial))
summary(as.factor(test$QC_NBM))
summary(as.factor(test$QC_CIVET))

###Addressing Translational reviews###

choline <- read_tsv('Choline_with_CRLB.csv')

DDD <- read_tsv("Demographics_DDD_2021.csv")
gf_DDD <- left_join(gf, DDD, by="ID")
gf_DDD_fep <- subset(gf_DDD, DX=="First Episode Patient")

tapply(gf_demo$SES.y, gf_demo$DX, mean, na.rm=TRUE)
tapply(gf_demo$SES.y, gf_demo$DX, sd, na.rm=TRUE)
summary(lm(gf_demo$SES.y ~ gf_demo$DX))

###Radial diffusivity tests

summary(lm(gf$RD_Ch123 ~ Age + Sex + CAST_total + Smoker + DX, data=gf))
summary(lm(gf$RD_L_Ch4 ~ Age + Sex + CAST_total + Smoker + DX, data=gf))
summary(lm(gf$RD_R_Ch4 ~ Age + Sex + CAST_total + Smoker + DX, data=gf))

cor.test(gf_HC$RD_L_Ch4, gf_HC$Choline)
cor.test(gf_HC$RD_R_Ch4, gf_HC$Choline)
cor.test(gf_HC$RD_Ch123, gf_HC$Choline)

cor.test(gf_FEP$RD_L_Ch4, gf_FEP$Choline)
cor.test(gf_FEP$RD_R_Ch4, gf_FEP$Choline)
cor.test(gf_FEP$RD_Ch123, gf_FEP$Choline)

###
summary(lm(Choline ~ Age + Sex + CAST_total + RD_Ch123, data=gf_HC))
summary(lm(Choline ~ Age + Sex + CAST_total + RD_L_Ch4, data=gf_HC))
summary(lm(Choline ~ Age + Sex + CAST_total + RD_R_Ch4, data=gf_HC))

summary(lm(Choline ~ Age + Sex + CAST_total + RD_Ch123, data=gf_FEP))
summary(lm(Choline ~ Age + Sex + CAST_total + RD_L_Ch4, data=gf_FEP))
summary(lm(Choline ~ Age + Sex + CAST_total + RD_R_Ch4, data=gf_FEP))

summary(lm(Choline ~ Age + Sex + CAST_total + Smoker + RD_Ch123, data=gf_FEP))
summary(lm(Choline ~ Age + Sex + CAST_total + Smoker + RD_L_Ch4, data=gf_FEP))
summary(lm(Choline ~ Age + Sex + CAST_total + Smoker + RD_R_Ch4, data=gf_FEP))

summary(lm(Choline ~ Age+ Sex + RD_Ch123*DX, data=gf))
summary(lm(Choline ~ Age+ Sex + RD_L_Ch4*DX, data=gf))
summary(lm(Choline ~ Age+ Sex + RD_R_Ch4*DX, data=gf))

###Relation to clinical scores

cor.test(gf_FEP$RD_L_Ch4, gf_FEP$totalN)
cor.test(gf_FEP$RD_R_Ch4, gf_FEP$totalN)
cor.test(gf_FEP$RD_Ch123, gf_FEP$totalN)

cor.test(gf_FEP$RD_L_Ch4, gf_FEP$YMRS_total)
cor.test(gf_FEP$RD_R_Ch4, gf_FEP$YMRS_total)
cor.test(gf_FEP$RD_Ch123, gf_FEP$YMRS_total)
