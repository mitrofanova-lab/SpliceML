# install.packages() - command for any library to be installed on your computer

library(readr)
library(dplyr)
library(survival)
library(survminer)
library(glmnet)
library(tidyverse) 
library(tidyr)
library(car)
library(usdm)
library(ggplot2)
library(reshape2)
library(knitr)
library(DT) 
library(caret) 
library(rsample) 
library(gplots)
library(UpSetR)
library(fastDummies)
library(rbsurv)
library(survivalROC)
options(warn = -1)

risk_score_calc_function <- function(coef_df, multi_cox_df){
  coef_df_columns <- gsub("`", "", grep('ENSG', coef_df$Event, value=T))
  multi_cox_df_columns <- grep('ENSG', colnames(multi_cox_df), value=T)
  
  same_columns <- setequal(coef_df_columns, multi_cox_df_columns)
  # Check whether we use same columns in coef_df and multi_cox_df
  if (same_columns) { 
    print("GOOD: The columns are the same.")
    print(length(coef_df_columns))
  } else { 
    print("WARNING:The columns are not the same") 
  }
  
  # Risk Score calculation
  sum_risk_score <- rep(0, nrow(multi_cox_df))
  
  for (event in grep('ENSG', coef_df$Event, value=T)){
    coef <- coef_df[coef_df$Event == event, 'coef']
    print(coef)
    psi_zvalue <- multi_cox_df[, gsub("`", "", event)]
    risk_score <- coef * psi_zvalue
    sum_risk_score <- sum_risk_score + risk_score
  }
  
  multi_cox_df$sum_risk_score <- as.vector(sum_risk_score)
  return(multi_cox_df)
}

# Read dataset of 
data <- readRDS('dataset.Rds')

# Remove AS events whose SD < 0.1
tmp_data <- data[, grep("^ENSG", names(data), value=TRUE)]
sd_values = sapply(tmp_data, function(x) if (is.numeric(x)) sd(x) else -Inf)

data_after_sd = tmp_data[, sd_values >= 0.1]
print(paste("Number of AS events more than cutoff threshold is", length(data_after_sd)))

data_after_sd_none = tmp_data[, sd_values < 0.1]
print(paste("Number of AS events less than cutoff threshold is", length(data_after_sd_none)))

# Plot number of AS events (Figure 2)
#pdf("AS_events_after_lowSD_removing.pdf", width=10)
plot_data <- data.frame(genes = colnames(data_after_sd))
plot_data <- tidyr::separate(plot_data, "genes", sep='_', into=c('ENSG', 'AS_type', 'chromosome'), remove=FALSE)
plot_data <- as.data.frame(table(plot_data$AS_type))
plot_data$Dataset <- "Data after SD"
ggplot(plot_data, aes(x = Freq, y = Var1)) + 
  geom_bar(stat='identity', position = "dodge", width=0.6, fill = 'skyblue4') +
  theme_bw() +
  labs(title = "",
       x = "Number of AS events",
       y = "AS event type") +
  geom_text(aes(label = Freq), position = position_dodge(width=0.7), vjust=0.5, hjust=-0.1, size=4) + 
  xlim(0, max(plot_data$Freq) + 100)
#dev.off()

# Scale PSI dataset and merge with clinical info
data_after_sd_scaled <- scale(data_after_sd)
data_to_split <- cbind(data[, c('PFI.time', 'PFI', 
                                'gender', 
                                'age_at_initial_pathologic_diagnosis', 
                                'stage',
                                'stage_groups_discrepance')],
                       data_after_sd_scaled)

# Split into training/testing cohort
tmp <- cbind(data[c('bcr_patient_barcode', 'PFI.time', 'PFI', 
                    # 'stage',
                    'stage_groups_discrepance', 'gender', 
                    'age_at_initial_pathologic_diagnosis')], 
             data_after_sd_scaled)
set.seed(957) 
tmp$strata <- with(tmp, paste(stage_groups_discrepance, sep='_'))
indexes <- createDataPartition(tmp$strata, p=0.5, list=F, times = 1)
train <- tmp[indexes, ]
test <- tmp[-indexes,  ]
dim(train)
dim(test)

print('Train PFI events amount:')
table(train$PFI)
print('Test PFI events amount:')
table(test$PFI)

print('Train stage groups amount:')
table(train$stage_groups_discrepance)
print('Test stage groups amount:')
table(test$stage_groups_discrepance)

print('Train gender amount:')
table(train$gender)
print('Test gender amount:')
table(test$gender)

train$strata <- NULL
test$strata <- NULL

paste("Mean age of train patients", round(mean(train$age_at_initial_pathologic_diagnosis), 2))
paste("Mean age of test patients", round(mean(test$age_at_initial_pathologic_diagnosis), 2))

paste("Standard Deviation of age in train patients", round(sd(train$age_at_initial_pathologic_diagnosis), 2))
paste("Standard Deviation of age in test patients", round(sd(test$age_at_initial_pathologic_diagnosis), 2))

# Cox PH analysis on training cohort
print("Cox PH analysis on training cohort")
train_4cox <- train
dim(train_4cox)
train_cox_results <- data.frame()
events <- grep('^ENSG', names(train_4cox), value = TRUE)

for (event in events){
  tmp = train_4cox[c('PFI.time', 'PFI',
                     'stage_groups_discrepance',
                     'gender',
                     'age_at_initial_pathologic_diagnosis',
                     event)]
  
  cph1 <- coxph(Surv(PFI.time, PFI) ~ ., data = tmp)
  summm <- summary(cph1)
  waldtest <- as_tibble(t(summm$waldtest)[3]) %>% rename(wald_pvalue=value)
  summm1 <- bind_cols(waldtest, as_tibble(summm$conf.int)[, c(3, 4)], as_tibble(summm$coefficients)) %>% mutate(gene=event)
  train_cox_results <- rbind(train_cox_results, summm1[4, ])
}

sign_train_events_cox <- subset(train_cox_results, `Pr(>|z|)` < 0.01) 
print(as_tibble(sign_train_events_cox), n=Inf) 
paste(nrow(sign_train_events_cox), '- Number of significant AS events after Cox PH (p<0.01)')

# Factorize binary clinical variables
train$gender_bin <- ifelse(train$gender == 'MALE', 1, 0) 
train$stage_groups_discrepance_bin <- ifelse(train$stage_groups_discrepance == 'Stages 1&2', 1, 0)
train$PFI.time <- ifelse(train$PFI.time == 0, 1, train$PFI.time)
test$PFI.time <- ifelse(test$PFI.time == 0, 1, test$PFI.time)
test$stage_groups_discrepance_bin <- ifelse(test$stage_groups_discrepance == 'Stages 1&2', 1, 0)
test$gender_bin <- ifelse(test$gender == 'MALE', 1, 0) 

# RBSURV modeling
print('RBSURV modeling')
our_events <- sign_train_events_cox$gene
# length(our_events)
n_genes <- length(our_events)
lasso_train <- cbind(train[, c('PFI.time', 'PFI')], train[, our_events])
z <- cbind(train$stage_groups_discrepance_bin, train$gender_bin, train$age_at_initial_pathologic_diagnosis)
colnames(z) <- c('stage', 'gender', 'age')

fit1 <- rbsurv(time = lasso_train$PFI.time,
               status = lasso_train$PFI,
               x = t(as.matrix(train[, our_events])),
               max.n.genes = n_genes,
               z=z,
               gene.ID = our_events,
               n.iter = 10,
               n.fold = 3,
               n.seq = 1)

s <- fit1$model
all_rbsurv_events <- s[((s$Selected == "*       ") & (s$Seq==1)), 'Gene']
print("AS events set with lowest AIC score in RBSURV modeling")
print(all_rbsurv_events)

# Plot RBSURV AIC plot
pdf("rbsurv_aic_plot.pdf", width=10)
# par(fig=c(0.46,0.98,0.48,1), pty="s", new=TRUE)
tmp <- which(fit1$model$Seq==1)
d <- fit1$model$AIC[tmp][1]
plot(fit1$model$Order[tmp], fit1$model$AIC[tmp] / d,
     type="l", pch=1,
     xlab="Genes", ylab="AIC")
title(paste("RBSURV"))
jj <- which.min(fit1$model$AIC[tmp])
text(fit1$model$Order[tmp][jj],fit1$model$AIC[tmp][jj] / d, "*", cex=2 )
dev.off()

# VIF analysis on defined AS events after rbsuvr modeling
usdm::vif(train[, all_rbsurv_events])


# Multivariable Cox PH analysis on training cohort using defined AS events after rbsurv modeling
train_multi_cox = train[, c('PFI.time', 'PFI', 'stage_groups_discrepance', 'gender_bin', 'age_at_initial_pathologic_diagnosis', all_rbsurv_events)]
cox_result = coxph(Surv(PFI.time, PFI) ~ ., data=train_multi_cox)
multicox_train_summary = data.frame(summary(cox_result)$coefficients)
multicox_train_summary = rownames_to_column(multicox_train_summary, var='Event')

sign_final_events <- gsub("`", "", subset(multicox_train_summary, `Pr...z..` < 0.05)$Event)
sign_final_events <- sign_final_events[grepl("^ENSG", sign_final_events)]
print('Final AS events after Multivariable Cox PH analysis')
print(sign_final_events)

# Multivariable Cox PH analysis on the testing cohort using only clinical variables
test_multi_cox_clin = test[, c('PFI.time', 'PFI', 'stage_groups_discrepance', 'gender','age_at_initial_pathologic_diagnosis')]
cox_model <- coxph(Surv(PFI.time, PFI) ~ ., data=test_multi_cox_clin)
# summary(cox_model)
risk_scores <- predict(cox_model, type = 'risk')
survcomp::concordance.index(x=risk_scores, surv.time = test_multi_cox_clin$PFI.time,
                            surv.event = test_multi_cox_clin$PFI, outx=F)$c.index
survival::concordance(cox_model, ranks = T)


# Multivariable Cox PH analysis on the testing cohort using only significant final defined list of AS events
test_multi_cox_no_adj = test[, c('PFI.time', 'PFI', sign_final_events)]
cox_model_no_adj = coxph(Surv(PFI.time, PFI) ~ ., data=test_multi_cox_no_adj)
# summary(cox_model_no_adj)
risk_scores_no_adj <- predict(cox_model_no_adj, type = 'risk')
# survcomp::concordance.index(x=risk_scores_no_adj, surv.time = test_multi_cox_no_adj$PFI.time, surv.event = test_multi_cox_no_adj$PFI, outx=F)$c.index

# Multivariable Cox PH analysis on the testing cohort using clinical variables + final defined list of AS events
test_multi_cox_adj = test[, c('PFI.time', 'PFI', 'stage_groups_discrepance', 'gender', 'age_at_initial_pathologic_diagnosis', sign_final_events)]
cox_model_adj = coxph(Surv(PFI.time, PFI) ~ ., data=test_multi_cox_adj)
summary(cox_model_adj)
risk_scores_adj <- predict(cox_model_adj, type = 'risk')
# survcomp::concordance.index(x=risk_scores_adj, surv.time = test_multi_cox_adj$PFI.time, surv.event = test_multi_cox_adj$PFI, outx=F)$c.index

# Multivariable Cox PH. Training cohort. Calculation of risk score
train_multi_cox_adj = train[, c('PFI.time', 'PFI', sign_final_events)]
cox_model <- coxph(Surv(PFI.time, PFI) ~., data=train_multi_cox_adj)
# summary(cox_model)

multicox_train_summary = data.frame(summary(cox_model)$coefficients)
multicox_train_summary = rownames_to_column(multicox_train_summary, var='Event')

train_multi_cox_adj <- risk_score_calc_function(multicox_train_summary, train_multi_cox_adj)
train_multi_cox_adj$risk_score_scaled <- (train_multi_cox_adj$sum_risk_score - min(train_multi_cox_adj$sum_risk_score)) / (max(train_multi_cox_adj$sum_risk_score) - min(train_multi_cox_adj$sum_risk_score))

# Plot risk score distribution
dens <- density(train_multi_cox_adj$risk_score_scaled)
dens <- data.frame(x = dens$x, y = dens$y)
dens <- dens[dens$x >= 0, ]
dens <- dens[dens$x <= 1, ]
dens$band <- ifelse(dens$x < 0.4, 0, ifelse(dens$x > 0.7, 2, 1))

pdf('density_plot_risk_score_train.pdf', width=8, height=5)
ggplot(dens, aes(x = x, y = y, fill = factor(band))) +
  geom_area(alpha = 0.9) +
  geom_line(color = "black", size = 0.5) +
  geom_vline(xintercept = c(0.4, 0.7), linetype = "dashed", color = "black", size = 0.5) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  scale_fill_manual(values = c("red", "green", "blue"),
                    name = "Risk groups",
                    labels = c("Low", "Intermediate", "High")) +
  labs(title = "Density Plot of Risk Score (training set)", x = "Risk score", y = "Density") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()

# Find p-values cutoffs
cutoffs <- numeric()
pvalues <- numeric()

train_min_riskscore <- min(train_multi_cox_adj$risk_score_scaled)
train_max_riskscore <- max(train_multi_cox_adj$risk_score_scaled) - 0.002

for (cutoff in seq(train_min_riskscore, train_max_riskscore, by=0.01)){
  # print(cutoff)
  train_multi_cox_adj$risk_group <- ifelse(train_multi_cox_adj$risk_score_scaled > cutoff, 1, 0)

  pvalue <- survdiff(Surv(PFI.time, PFI) ~ risk_group, data = train_multi_cox_adj)$pvalue
  # print(pvalue)

  cutoffs <- c(cutoffs, cutoff)
  pvalues <- c(pvalues, pvalue)
}

# Print plot of cutoffs and p-values from KM
results_df <- data.frame(Cutoff = cutoffs, P_Value = pvalues)
results_df$log2p <- -log2(results_df$P_Value)

pdf('risk_score_plot_of_pvalue_vs_cutoff_for_train.pdf', width=8, height = 5)
ggplot(results_df, aes(x = Cutoff, y = log2p)) +
  geom_line() + # Plot lines
  # geom_point() + # Add points
  labs(title = "P-Value by Risk Score Cutoff for train set, step=0.01", x = "Risk Score Cutoff", y = "-log2(P-Value)") +
  theme_minimal() +
  geom_hline(yintercept = -log2(0.05))
dev.off()

# Remove 16 sample which are high risk and run split on rest of data
train_multi_cox_adj_removed <- subset(train_multi_cox_adj, risk_group!=1)

cutoffs <- numeric()
pvalues <- numeric()

train_min_riskscore <- min(train_multi_cox_adj_removed$risk_score_scaled)
train_max_riskscore <- max(train_multi_cox_adj_removed$risk_score_scaled) - 0.002

for (cutoff in seq(train_min_riskscore, train_max_riskscore, by=0.01)){
  # print(cutoff)
  train_multi_cox_adj_removed$risk_group <- ifelse(train_multi_cox_adj_removed$risk_score_scaled > cutoff, 1, 0)

  pvalue <- survdiff(Surv(PFI.time, PFI) ~ risk_group, data = train_multi_cox_adj_removed)$pvalue
  # print(pvalue)

  cutoffs <- c(cutoffs, cutoff)
  pvalues <- c(pvalues, pvalue)
}

# Print plot of cutoffs and p-values from KM
results_df <- data.frame(Cutoff = cutoffs, P_Value = pvalues)
results_df$log2p <- -log2(results_df$P_Value)

# pdf('risk_score_plot_of_pvalue_vs_cutoff_for_train_removing_16highrisks.pdf', width=8, height = 5)
# ggplot(results_df, aes(x = Cutoff, y = log2p)) +
#   geom_line() + # Plot lines
#   labs(title = "P-Value by Risk Score Cutoff for train set, step=0.01. \nremove 16 high risk groups. and run again analysis", x = "Risk Score Cutoff", y = "-log2(P-Value)") +
#   theme_minimal() +
#   geom_hline(yintercept = -log2(0.05))
# dev.off()

test_multi_cox_adj = test[, c('PFI.time', 'PFI', 
                              # 'stage_groups_discrepance',
                              # 'gender',
                              # 'age_at_initial_pathologic_diagnosis',
                              sign_final_events)]

cox_model_adj = coxph(Surv(PFI.time, PFI) ~ ., data=test_multi_cox_adj)
# summary(cox_model_adj)
multicox_test_summary = data.frame(summary(cox_model_adj)$coefficients)
multicox_test_summary = rownames_to_column(multicox_test_summary, var='Event')

test_multi_cox_adj <- risk_score_calc_function(multicox_test_summary, test_multi_cox_adj)

test_multi_cox_adj$risk_score_scaled <- 
  (test_multi_cox_adj$sum_risk_score - min(test_multi_cox_adj$sum_risk_score)) /
  (max(test_multi_cox_adj$sum_risk_score) - min(test_multi_cox_adj$sum_risk_score))


# Plot risk score distribution and division into risk groups
dens <- density(test_multi_cox_adj$risk_score_scaled)
dens <- data.frame(x = dens$x, y = dens$y)
dens <- dens[dens$x >= 0, ]
dens <- dens[dens$x <= 1, ]
dens$band <- ifelse(dens$x < 0.4, 0, ifelse(dens$x > 0.7, 2, 1))

pdf('density_plot_risk_score_test.pdf', width=8, height=5)
ggplot(dens, aes(x = x, y = y, fill = factor(band))) +
  geom_area(alpha = 0.9) +
  geom_line(color = "black", size = 0.5) +
  geom_vline(xintercept = c(0.4, 0.7), linetype = "dashed", color = "black", size = 0.5) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  scale_fill_manual(values = c("red", "green", "blue"),
                    name = "Risk groups",
                    labels = c("Low", "Intermediate", "High")) +
  labs(title = "Density Plot of Risk Score (testing set)", x = "Risk score", y = "Density") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
dev.off()

