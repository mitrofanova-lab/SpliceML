# Install packages
library(survival)
library(survminer)
library(pheatmap)
library(RColorBrewer)
library(data.table)
library(stringr)
# install.packages('ggfortify')
library(ggfortify)
# install.packages("timeROC")
library(timeROC)


# Read data
target_psi_df = readRDS('psi_data_AC_ICAM.rds')
final_clinical_df = readRDS('clinical_data_AC_ICAM.rds')

# Extract identified AS events from TCGA analysis
targetSpliceEvents <- target_psi_df[target_psi_df$V1 %in% 
                                      c('ENSG00000107779;A3:10:86924641-86929297:86924641-86929300:+', 
                                        'ENSG00000167792;A5:11:67611978-67612120:67611569-67612120:+', 
                                        'ENSG00000204271;A3:X:56978498-56978580:56978421-56978580:-', 
                                        'ENSG00000241553;A5:3:9793181-9797659:9793174-9797659:+', 
                                        'ENSG00000196778;A3:11:4483176-4488573:4483176-4488953:+'), ]

targetSpliceEvents$ENSG <- NULL
targetSpliceEvents$event <- NULL

# Preprocess targetSpliceEvents PSI data
t.targetSpliceEvents <- as.data.frame(t(targetSpliceEvents[, -1]))
colnames(t.targetSpliceEvents) <- targetSpliceEvents$V1
colnames(t.targetSpliceEvents)[1] <- 'BMPR1A'
colnames(t.targetSpliceEvents)[2] <- 'NDUFV1'
colnames(t.targetSpliceEvents)[3] <- 'SPIN3'
colnames(t.targetSpliceEvents)[4] <- 'ARPC4'
colnames(t.targetSpliceEvents)[5] <- 'OR52K1'

t.targetSpliceEvents <- tibble::rownames_to_column(t.targetSpliceEvents, var='Run1')


# Merge PSI data and clincal meta data
merged_df <- merge(final_clinical_df, t.targetSpliceEvents, by.x = 'Run', by.y = 'Run1')

# SCALING
merged_df$BMPR1A <- scale(merged_df$BMPR1A)[,1]
merged_df$SPIN3 <- scale(merged_df$SPIN3)[,1]
merged_df$NDUFV1 <- scale(merged_df$NDUFV1)[,1]
merged_df$ARPC4 <- scale(merged_df$ARPC4)[,1]
merged_df$OR52K1 <- scale(merged_df$OR52K1)[,1]


cox_model <- coxph(Surv(time, status) ~ BMPR1A + NDUFV1 + SPIN3 + ARPC4 + OR52K1, data = merged_df)
summ <- summary(cox_model)
# non_scaled_coeffs = summ$coefficients[, 1]
coeffs <- summary(cox_model)$coefficients[, 1]
risk_score <- as.numeric(as.matrix(merged_df[, names(coeffs)]) %*% coeffs)
# risk_score = as.numeric(as.matrix(merged_df[, names(coeffs)]) %*% non_scaled_coeffs)

merged_df$risk_score <- risk_score

riskscore_model <- coxph(Surv(PFS_MONTHS, status) ~ risk_score, data = merged_df)
summary(riskscore_model)

# Risk score Kaplan Meier curve
merged_df$scaled_riskscore <- (merged_df$risk_score - min(merged_df$risk_score)) / (max(merged_df$risk_score) - min(merged_df$risk_score))

merged_df$risk_group <- cut(merged_df$scaled_riskscore,
                            breaks = c(-Inf, mean(merged_df$scaled_riskscore), Inf),
                            labels = c("Low", "High"))
survdiff(Surv(time, status) ~ risk_group, data=merged_df)
km_fit <- survfit(Surv(time, status) ~ risk_group, data = merged_df)
km_fit

km_plot <- ggsurvplot(km_fit, data = merged_df,
                      pval = TRUE,
                      risk.table = TRUE,
                      ggtheme = theme_minimal(),
                      risk.table.height = 0.25,
                      legend.title = "",
                      risk.table.fontsize = 3,
                      pval.method = T,
                      palette = c('red', 'green')
) + labs(x = 'Days', y = 'SurvProb (PFS)')

km_plot

# Density plot
# Plot risk score distribution and division into risk groups
dens <- density(merged_df$scaled_riskscore)
dens <- data.frame(x = dens$x, y = dens$y)
# dens <- dens[dens$x >= 0, ]
# dens <- dens[dens$x <= 1, ]
dens$band <- ifelse(dens$x < mean(merged_df$scaled_riskscore), 0, 1)

dens_plot = ggplot(dens, aes(x = x, y = y, fill = factor(band))) +
  geom_area(alpha = 0.9) +
  geom_line(color = "black", size = 0.5) +
  geom_vline(xintercept = c(mean(merged_df$scaled_riskscore)), linetype = "dashed", color = "black", size = 0.5) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  scale_fill_manual(values = c("red", "green"),
                    name = "Risk groups",
                    labels = c("Low", "High")) +
  labs(title = "Density Plot of Risk Score", x = "Risk score", y = "Density") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

dens_plot



# Heatmaps
heatmap_data <- as.data.frame(merged_df[, c('OR52K1', 'SPIN3', 'NDUFV1', 'BMPR1A', 'ARPC4', 'risk_group')])
heatmap_data_ordered <- heatmap_data[order(heatmap_data$risk_group), ]
heatmap_data_ordered$risk_group <- NULL
heatmap_data_ordered <- t(heatmap_data_ordered)

patient_info <- data.frame(risk_group = merged_df$risk_group)
row.names(patient_info) <- rownames(merged_df)

# set up color scaling to be symmetrical around zero value
min_val <- min(heatmap_data_ordered)
max_val <- max(heatmap_data_ordered)
abs_max <- max(abs(min_val), abs(max_val))
abs_max = 2.5
min_val <- -abs_max
max_val <- abs_max
num_colors <- 51 
my_color_palette <- colorRampPalette(c("blue", "white", "red"))
colors <- my_color_palette(num_colors)
breaks <- seq(min_val, max_val, length.out = num_colors + 1)

annotation_colors <- list(risk_group = c("Low" = "red", "High" = "green"))

heat_plot = pheatmap(heatmap_data_ordered, annotation_col = patient_info,
                     cluster_cols = FALSE, cluster_rows = FALSE, 
                     show_colnames = FALSE, border_color = 'black', 
                     color = colors,
                     breaks = breaks,
                     main = 'Risk groups', 
                     legend_breaks = c(-2, 0, 2), 
                     annotation_colors = annotation_colors
)

heat_plot

low_group_data <- merged_df[merged_df$risk_group == "Low", ]
heatmap_data_low <- as.data.frame(low_group_data[, c('OR52K1', 'SPIN3', 'NDUFV1', 'BMPR1A', 'ARPC4')])
heatmap_data_low <- t(heatmap_data_low)

low_heatplot = pheatmap(heatmap_data_low, 
                        cluster_cols = T, cluster_rows = F, 
                        show_colnames = F, border_color = NA,   
                        color = colors, 
                        breaks = breaks,
                        legend_breaks = c(-2, 0, 2), 
                        main = 'Low risk group (ICAM data)'
)

low_heatplot


high_group_data <- merged_df[merged_df$risk_group == "High", ]
heatmap_data_high <- as.data.frame(high_group_data[, c('OR52K1', 'SPIN3', 'NDUFV1', 'BMPR1A', 'ARPC4')])
heatmap_data_high <- t(heatmap_data_high)

high_heatplot = pheatmap(heatmap_data_high, 
                         cluster_cols = T, cluster_rows = F, 
                         show_colnames = F, border_color = NA,   
                         color = colors, 
                         breaks = breaks,
                         legend_breaks = c(-2, 0, 2), 
                         main = 'High risk group (ICAM data)'
)

high_heatplot

# Comparison to known transcriptomic markers
msi_info = read.table('AC-ICAM/MSI_AC_ICAM.csv', header=T, sep = ';')

msi_info = merge(msi_info, final_clinical[, c('patientId', 'PFS_MONTHS', 'PFS_STATUS')], by.x = 'Patient_ID', by.y = 'patientId')
msi_info$MANTIS.score = sub(',', '.', msi_info$MANTIS.score)
msi_info$MANTIS.score <- as.numeric(msi_info$MANTIS.score)
msi_info$status = ifelse(msi_info$PFS_STATUS == '0:DiseaseFree', 0, 1)

# MANTIS/MSI score without risk score
summary(coxph(Surv(PFS_MONTHS, status) ~ MANTIS.score, data = msi_info))

# MANTIS/MSI score with risk score
# merged_df = readRDS(file = 'AC_ICAM_merged_df_4other_markers_analysis_and_4plots.rds')
merged_df = merged_df[, c('Run', 'status', 'time', 'risk_score', 'scaled_riskscore')]
merged_df = merge(merged_df, final_clinical[, c('Run', 'patientId')], by = 'Run')

msi_info_ext = merge(merged_df, msi_info, by.x = 'patientId', by.y = 'Patient_ID')
summary(coxph(Surv(PFS_MONTHS, status.x) ~ MANTIS.score + risk_score, data = msi_info_ext))


transcriptomic_markers = c('RUNX1', 'MMP2', 'TP53', 'VEGFA', 'CXCL12', 'CXCR4', 'MMP7', 'MMP9', 'CDX2', 'EGFR', 
                           'SOX9', 'ERBB2', 'MYC', 'CDH1', 'BCL2', 'CD44', 'SMAD4', 'EPHA2', 'FOXM1')

# Take vsd normalized counts of AC-ICAM dataset
vst_mat = readRDS('mrna_expr_AC_ICAM.rds')
vst_df = as.data.frame(t(vst_mat))
# Convert ENSG_ids to gene_names
library(biomaRt)
my_list = colnames(vst_df)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
bm = getBM(attributes = c('ensembl_gene_id','hgnc_symbol'), filters = 'ensembl_gene_id', values = my_list, mart = ensembl)

mapping = setNames(bm$hgnc_symbol, bm$ensembl_gene_id)
new_colnames = mapping[colnames(vst_df)]
colnames(vst_df) = new_colnames

# Selected our target gene_names from vsd counts
vst_df = vst_df[, transcriptomic_markers]
# and Merge with survival data
markers_df = merge(merged_df[, c('Run', 'AGE_AT_DX', 'stage_groups', 'gender', 'time', 'status', 'risk_score', 'scaled_riskscore')], 
                   vst_df, 
                   by.x = 'Run', 
                   by.y = 0)

# Run univariavle cox analysis on transcriptomic markers and transc.markers + risk_score
markerExpr_coxoutput = data.frame()
for (marker in transcriptomic_markers){
  tmp_df = markers_df[, c('time', 'status', marker)]
  cph1 = coxph(Surv(time, status) ~ ., data = tmp_df)
  summm <- summary(cph1)
  waldtest <- as_tibble(t(summm$waldtest)[3]) %>% rename(wald_pvalue=value)
  summm1 <- bind_cols(waldtest, as_tibble(summm$conf.int)[, c(3, 4)], as_tibble(summm$coefficients)) %>% mutate(gene=marker)
  markerExpr_coxoutput <- rbind(markerExpr_coxoutput, summm1[1, ])
}

# Run univariable cox + RiskScore
markerExpr_riskscore_coxoutput = data.frame()
for (marker in transcriptomic_markers){
  tmp_df = markers_df[, c('time', 'status', 'scaled_riskscore', marker)]
  cph1 = coxph(Surv(time, status) ~ ., data = tmp_df)
  summm = summary(cph1)
  waldtest = as_tibble(t(summm$waldtest)[3]) %>% rename(wald_pvalue=value)
  tstat = as_tibble(t(summm$waldtest)[1]) %>% rename(tstat=value)
  summm1 = bind_cols(tstat, waldtest, as_tibble(summm$conf.int)[, c(3, 4)], as_tibble(summm$coefficients)) %>% mutate(gene=marker)
  markerExpr_riskscore_coxoutput = rbind(markerExpr_riskscore_coxoutput, summm1[1,])
}


plotdata_othermarkers_ICAM <- data.frame(
  Method = c("RUNX1", "MMP2", "TP53", "VEGF", "CXCL12", "CXCR4", "MMP7", "MMP9", "CDX2", "EGFR", "SOX9", "ERBB2", "MYC", 
             "CDH1","BCL2", "CD44", "SMAD4", "EPHA2", "FOXM1", 'MSI'),
  without_riskscore = c(markerExpr_coxoutput$wald_pvalue, 0.4404541),
  with_riskscore = c(markerExpr_riskscore_coxoutput$wald_pvalue, 0.00364844))

desired_order <- rev(c("RUNX1", "MMP2", "TP53", "VEGF", "CXCL12", "CXCR4", "MMP7", "MMP9", "CDX2", "EGFR", "SOX9", "ERBB2", "MYC", 
                       "CDH1", "BCL2", "CD44", "SMAD4", "EPHA2", "FOXM1", 'MSI'))
data_long <- melt(plotdata_othermarkers_ICAM, id.vars = "Method", variable.name = "Model", value.name = "p_value")
data_long$Method <- factor(data_long$Method, levels = desired_order)
data_long$Model <- factor(data_long$Model, levels = c("with_riskscore", "without_riskscore"))

ggplot(data_long, aes(x=Method, y=-log2(p_value), fill = Model)) + 
  geom_col(position = 'dodge') +
  coord_flip() +
  labs(title = "P-values for markers with and without Risk Score", 
       x = "Marker", y = "-log2(p-value)") +
  theme_minimal() +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, face = 'bold'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
  ) +
  geom_hline(yintercept = -log2(0.05), linetype = 'dashed', color = 'black') +
  scale_fill_manual(values = c("without_riskscore" = "blue", "with_riskscore" = "red")) +
  scale_y_continuous(breaks = seq(0, 10, by=5))
# scale_y_continuous(breaks = c(0, 5, 10, 15))


# Comparison to genomic markers
icam_cna <- read.table('coad_silu_2022/data_cna.txt', header = T, sep='\t')
colnames(icam_cna) = sub('.PT.01', '', colnames(icam_cna))
colnames(icam_cna) = gsub('\\.', '-', colnames(icam_cna))

genomic_markers = c('TP53', 'TGFBR2', 'SOX9', 'SMAD4', 'PTEN', 'PIK3CA', 'NRAS', 'KRAS', 'IDH2', 'IDH1', 'FBXW7', 'ERBB2', 
                    'CTNNB1', 'CDK12', 'BRCA2', 'BRCA1', 'APC', 'AKT1')

icam_cna = icam_cna[icam_cna$Hugo_Symbol %in% genomic_markers, ]
rownames(icam_cna) = NULL
icam_cna = tibble::column_to_rownames(icam_cna, var='Hugo_Symbol')
icam_cna = as.data.frame(t(icam_cna))

icam_survival_cna = merge(msi_info_ext, icam_cna, by.x = 'patientId', by.y = 0)


cna_diploid_wo_riskscore <- data.frame(Gene = character(), Diploid_pvalue = numeric(), 
                                       Diploid_HR = numeric(),
                                       Diploid_Lower95 = numeric(), 
                                       Diploid_Upper95 = numeric(),
                                       statistic = numeric(),
                                       stringsAsFactors = FALSE)

for (marker in genomic_markers) {
  tmp = icam_survival_cna[, c('time', 'status.x', marker)]
  tmp$diploid <- ifelse(as.numeric(tmp[[marker]]) == 0, 'diploid', 'no_diploid')
  cox_diploid <- summary(coxph(Surv(time, status.x) ~ diploid, tmp))
  cna_diploid_wo_riskscore <- rbind(cna_diploid_wo_riskscore, data.frame(Gene = marker,
                                                                         Diploid_pvalue = round(cox_diploid$waldtest[3], 5),
                                                                         Diploid_HR = cox_diploid$coefficients[2],
                                                                         Diploid_Lower95 = cox_diploid$conf.int[3], 
                                                                         Diploid_Upper95 = cox_diploid$conf.int[4],
                                                                         statistic = cox_diploid$waldtest[1]))
}

cna_diploid_w_riskscore <- data.frame(Gene = character(), Diploid_pvalue = numeric(), 
                                      Diploid_HR = numeric(),
                                      Diploid_Lower95 = numeric(), 
                                      Diploid_Upper95 = numeric(), 
                                      statistic = numeric(), 
                                      stringsAsFactors = FALSE)

for (marker in genomic_markers) {
  tmp = icam_survival_cna[, c('time', 'status.x', 'scaled_riskscore', marker)]
  tmp$diploid <- ifelse(as.numeric(tmp[[marker]]) == 0, 'diploid', 'no_diploid')
  tmp[[marker]] = NULL
  cox_diploid <- summary(coxph(Surv(time, status.x) ~ ., tmp))
  cna_diploid_w_riskscore <- rbind(cna_diploid_w_riskscore, data.frame(Gene = marker,
                                                                       Diploid_pvalue = round(cox_diploid$waldtest[3], 5),
                                                                       Diploid_HR = cox_diploid$coefficients[2, 2],
                                                                       Diploid_Lower95 = cox_diploid$conf.int[2, 3], 
                                                                       Diploid_Upper95 = cox_diploid$conf.int[2, 4], 
                                                                       statistic = cox_diploid$waldtest[1]))
}

### GAIN w/o risk score
cna_gain_wo_riskscore <- data.frame(Gene = character(), gain_pvalue = numeric(), 
                                    gain_HR = numeric(),
                                    gain_Lower95 = numeric(), 
                                    gain_Upper95 = numeric(), 
                                    statistic = numeric(),
                                    stringsAsFactors = FALSE)
for (marker in genomic_markers) {
  tmp = icam_survival_cna[, c('time', 'status.x', marker)]
  tmp$gain <- ifelse(as.numeric(tmp[[marker]]) %in% c(1, 2), 'gain', 'no_gain')
  cox_gain <- summary(coxph(Surv(time, status.x) ~ gain, tmp))
  cna_gain_wo_riskscore <- rbind(cna_gain_wo_riskscore, data.frame(Gene = marker,
                                                                   gain_pvalue = round(cox_gain$waldtest[3], 5),
                                                                   gain_HR = cox_gain$coefficients[2],
                                                                   gain_Lower95 = cox_gain$conf.int[3], 
                                                                   gain_Upper95 = cox_gain$conf.int[4], 
                                                                   statistic = cox_gain$waldtest[1]))
}
### GAIN w risk score
cna_gain_w_riskscore <- data.frame(Gene = character(), gain_pvalue = numeric(), 
                                   gain_HR = numeric(),
                                   gain_Lower95 = numeric(), 
                                   gain_Upper95 = numeric(), 
                                   statistic = numeric(),
                                   stringsAsFactors = FALSE)
for (marker in genomic_markers) {
  tmp = icam_survival_cna[, c('time', 'status.x', 'scaled_riskscore', marker)]
  tmp$gain <- ifelse(as.numeric(tmp[[marker]]) %in% c(1, 2), 'gain', 'no_gain')
  tmp[[marker]] = NULL
  cox_gain <- summary(coxph(Surv(time, status.x) ~ ., tmp))
  cna_gain_w_riskscore <- rbind(cna_gain_w_riskscore, data.frame(Gene = marker,
                                                                 gain_pvalue = round(cox_gain$waldtest[3], 5),
                                                                 gain_HR = cox_gain$coefficients[2, 2],
                                                                 gain_Lower95 = cox_gain$conf.int[2, 3], 
                                                                 gain_Upper95 = cox_gain$conf.int[2, 4], 
                                                                 statistic = cox_gain$waldtest[1]))
}

### DELETION w/o risk score
cna_deletion_wo_riskscore <- data.frame(Gene = character(), del_pvalue = numeric(), 
                                        del_HR = numeric(),
                                        del_Lower95 = numeric(), 
                                        del_Upper95 = numeric(), 
                                        statistic = numeric(),
                                        stringsAsFactors = FALSE)
for (marker in genomic_markers) {
  tmp = icam_survival_cna[, c('time', 'status.x', marker)]
  tmp$deletion <- ifelse(as.numeric(tmp[[marker]]) %in% c(-1, -2), 'deletion', 'no_deletion')
  cox_deletion <- summary(coxph(Surv(time, status.x) ~ deletion, tmp))
  cna_deletion_wo_riskscore <- rbind(cna_deletion_wo_riskscore, data.frame(Gene = marker,
                                                                           del_pvalue = round(cox_deletion$waldtest[3], 5),
                                                                           del_HR = cox_deletion$coefficients[2],
                                                                           del_Lower95 = cox_deletion$conf.int[3], 
                                                                           del_Upper95 = cox_deletion$conf.int[4], 
                                                                           statistic = cox_deletion$waldtest[1]
  ))
}

### DELETION w risk score
cna_deletion_w_riskscore <- data.frame(Gene = character(), del_pvalue = numeric(), 
                                       del_HR = numeric(),
                                       del_Lower95 = numeric(), 
                                       del_Upper95 = numeric(), 
                                       statistic = numeric(),
                                       stringsAsFactors = FALSE)
for (marker in genomic_markers) {
  tmp = icam_survival_cna[, c('time', 'status.x', 'scaled_riskscore', marker)]
  tmp$deletion <- ifelse(as.numeric(tmp[[marker]]) %in% c(-1, -2), 'deletion', 'no_deletion')
  tmp[[marker]] = NULL
  cox_deletion <- summary(coxph(Surv(time, status.x) ~ ., tmp))
  cna_deletion_w_riskscore <- rbind(cna_deletion_w_riskscore, data.frame(Gene = marker,
                                                                         del_pvalue = round(cox_deletion$waldtest[3], 5),
                                                                         del_HR = cox_deletion$coefficients[2, 2],
                                                                         del_Lower95 = cox_deletion$conf.int[2, 3], 
                                                                         del_Upper95 = cox_deletion$conf.int[2, 4], 
                                                                         statistic = cox_deletion$waldtest[1]))
}


### DIPLOID
plotdata_diploid_ICAM <- data.frame(
  Method = c('TP53', 'TGFBR2', 'SOX9', 'SMAD4', 'PTEN', 'PIK3CA', 'NRAS', 'KRAS', 'IDH2', 'IDH1', 'FBXW7', 'ERBB2', 
             'CTNNB1', 'CDK12', 'BRCA2', 'BRCA1', 'APC', 'AKT1'),
  without_riskscore = c(cna_diploid_wo_riskscore$Diploid_pvalue),
  with_riskscore = c(cna_diploid_w_riskscore$Diploid_pvalue))

data_long <- melt(plotdata_diploid_ICAM, id.vars = "Method", variable.name = "Model", value.name = "p_value")
data_long$Model <- factor(data_long$Model, levels = c("with_riskscore", "without_riskscore"))


ggplot(data_long, aes(x=Method, y=-log2(p_value), fill = Model)) + 
  geom_col(position = 'dodge') +
  coord_flip() +
  labs(title = "P-values for genomic markers (diploid) with and without Risk Score", 
       x = "Marker", y = "-log2(p-value)") +
  theme_minimal() +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, face = 'bold'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
  ) +
  geom_hline(yintercept = -log2(0.05), linetype = 'dashed', color = 'black') +
  scale_fill_manual(values = c("without_riskscore" = "blue", "with_riskscore" = "red")) 
# scale_y_continuous(breaks = seq(0, 10, by=5))
# scale_y_continuous(breaks = c(0, 5, 10, 15))



### GAIN
plotdata_gain_ICAM <- data.frame(
  Method = c('TP53', 'TGFBR2', 'SOX9', 'SMAD4', 'PTEN', 'PIK3CA', 'NRAS', 'KRAS', 'IDH2', 'IDH1', 'FBXW7', 'ERBB2', 
             'CTNNB1', 'CDK12', 'BRCA2', 'BRCA1', 'APC', 'AKT1'),
  without_riskscore = c(cna_gain_wo_riskscore$gain_pvalue),
  with_riskscore = c(cna_gain_w_riskscore$gain_pvalue))

data_long <- melt(plotdata_gain_ICAM, id.vars = "Method", variable.name = "Model", value.name = "p_value")
data_long$Model <- factor(data_long$Model, levels = c("with_riskscore", "without_riskscore"))


ggplot(data_long, aes(x=Method, y=-log2(p_value), fill = Model)) + 
  geom_col(position = 'dodge') +
  coord_flip() +
  labs(title = "P-values for genomic markers (gain) with and without Risk Score", 
       x = "Marker", y = "-log2(p-value)") +
  theme_minimal() +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, face = 'bold'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
  ) +
  geom_hline(yintercept = -log2(0.05), linetype = 'dashed', color = 'black') +
  scale_fill_manual(values = c("without_riskscore" = "blue", "with_riskscore" = "red")) 


### DELETION
plotdata_deletion_ICAM <- data.frame(
  Method = c('TP53', 'TGFBR2', 'SOX9', 'SMAD4', 'PTEN', 'PIK3CA', 'NRAS', 'KRAS', 'IDH2', 'IDH1', 'FBXW7', 'ERBB2', 
             'CTNNB1', 'CDK12', 'BRCA2', 'BRCA1', 'APC', 'AKT1'),
  without_riskscore = c(cna_deletion_wo_riskscore$del_pvalue),
  with_riskscore = c(cna_deletion_w_riskscore$del_pvalue))

data_long <- melt(plotdata_deletion_ICAM, id.vars = "Method", variable.name = "Model", value.name = "p_value")
data_long$Model <- factor(data_long$Model, levels = c("with_riskscore", "without_riskscore"))

ggplot(data_long, aes(x=Method, y=-log2(p_value), fill = Model)) + 
  geom_col(position = 'dodge') +
  coord_flip() +
  labs(title = "P-values for genomic markers (deletion) with and without Risk Score", 
       x = "Marker", y = "-log2(p-value)") +
  theme_minimal() +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, face = 'bold'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
  ) +
  geom_hline(yintercept = -log2(0.05), linetype = 'dashed', color = 'black') +
  scale_fill_manual(values = c("without_riskscore" = "blue", "with_riskscore" = "red")) 
