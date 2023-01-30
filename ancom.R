if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ANCOMBC")

library(ANCOMBC)
library(tidyverse)
library(caret)
library(DT)
out = ancombc2(data =TM1_TM4_Fdata_Gut, tax_level = "Genus", fix_formula = "Treatment..", p_adj_method = "holm") 

res = out$res


data("atlas1006")
tse = atlas1006[, atlas1006$time == 0]

tse$bmi = recode(tse$bmi_group,
                 obese = "obese",
                 severeobese = "obese",
                 morbidobese = "obese")

tse = tse[, tse$bmi %in% c("lean", "overweight", "obese")]

tse$bmi = factor(tse$bmi, levels = c("obese", "overweight", "lean"))

tse$region = recode(as.character(tse$nationality),
                    Scandinavia = "NE", UKIE = "NE", SouthEurope = "SE", 
                    CentralEurope = "CE", EasternEurope = "EE",
                    .missing = "unknown")

tse = tse[, ! tse$region %in% c("EE", "unknown")]

set.seed(123)
test_output = ancombc2(data = tse, assay_name = "counts", tax_level = "Family",
                  fix_formula = "age + region + bmi", rand_formula = NULL,
                  p_adj_method = "holm", pseudo = 0, pseudo_sens = TRUE,
                  prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                  group = "bmi", struc_zero = TRUE, neg_lb = TRUE,
                  alpha = 0.05, n_cl = 2, verbose = TRUE,
                  global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = TRUE,
                  iter_control = list(tol = 1e-2, max_iter = 20, 
                                      verbose = TRUE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  lme_control = lme4::lmerControl(),
                  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                  trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                              nrow = 2, 
                                                              byrow = TRUE),
                                                       matrix(c(-1, 0, 1, -1),
                                                              nrow = 2, 
                                                              byrow = TRUE)),
                                       node = list(2, 2),
                                       solver = "ECOS",
                                       B = 100))

test_tab_zero = test_output$zero_ind
test_tab_zero %>%
  datatable(caption = "The detection of structural zeros")


test_tab_sens = test_output$pseudo_sens_tab
test_tab_sens %>%
  datatable(caption = "Sensitivity Scores") %>%
  formatRound(colnames(test_tab_sens)[-1], digits = 2)

test_res_prim = test_output$res

#Select columns that only have age in it along with the taxon column
test_df_age = test_res_prim %>%
  dplyr::select(taxon, ends_with("age")) 

#Select rows that have TRUE (1) in the diff_age column i.e. differentially expressed
test_df_fig_age = test_df_age %>%
  filter(diff_age == 1) %>% 
  #Arrange the table by highest to lowest log fold change
  arrange(desc(lfc_age)) %>%
  #Create a new column called direct that says Positive LFC for and Negative LFC for appropriate lfc_age values
  mutate(direct = ifelse(lfc_age > 0, "Positive LFC", "Negative LFC"))

#Change taxon to a factor and keep the levels the same
test_df_fig_age$taxon = factor(test_df_fig_age$taxon, levels = test_df_fig_age$taxon)

#Change direct to a factor and change the levels
test_df_fig_age$direct = factor(test_df_fig_age$direct, 
                           levels = c("Positive LFC", "Negative LFC"))

#Create a bar graph with taxon by lfc_age filled by direction
test_fig_age = test_df_fig_age %>%
  ggplot(aes(x = taxon, y = lfc_age, fill = direct)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  #Add error bars based upon standard error
  geom_errorbar(aes(ymin = lfc_age - se_age, ymax = lfc_age + se_age), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "Log fold changes as one unit increase of age") + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1))
test_fig_age

test_res_dunn = test_output$res_dunn

#Keep items that were differentially expressed (TRUE)
test_df_fig_dunn = test_res_dunn %>%
  dplyr::filter(diff_bmilean == 1 | diff_bmioverweight == 1) %>%
  #Add columns with corresponding lfc value if differentially expressed, and 0 if it is not
  mutate(lfc_lean = ifelse(diff_bmilean == 1, 
                           lfc_bmilean, 0),
         lfc_overweight = ifelse(diff_bmioverweight == 1, 
                                 lfc_bmioverweight, 0)) %>%
  #keep the taxon column, convert two previous columns to new columns based upon them with two sig figs
  transmute(taxon, 
            `Lean vs. Obese` = round(lfc_lean, 2), 
            `Overweight vs. Obese` = round(lfc_overweight, 2)) %>%
  #Convert from wide to long converting the two columns made above to a group column and their values to a value column
  pivot_longer(cols = `Lean vs. Obese`:`Overweight vs. Obese`, 
               names_to = "group", values_to = "value") %>%
  arrange(taxon)

#Find the miniumum value and keep the nearest whole number below it
lo = floor(min(test_df_fig_dunn$value))

#Find the maximum value and keep the nearest whole number above it
up = ceiling(max(test_df_fig_dunn$value))

#Find the value that is halfway between the two above values
mid = (lo + up)/2


test_fig_dunn = test_df_fig_dunn %>%
  ggplot(aes(x = group, y = taxon, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(group, taxon, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL, title = "Log fold changes as compared to obese subjects") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
test_fig_dunn






set.seed(123)
output = ancombc2(data = Fdata_Gut, tax_level = "Genus",
                  fix_formula = "Treatment..", rand_formula = NULL,
                  p_adj_method = "holm", pseudo = 0, pseudo_sens = TRUE,
                  prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                  group = "Treatment..", struc_zero = TRUE, neg_lb = TRUE,
                  alpha = 0.05, n_cl = 2, verbose = TRUE,
                  global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = FALSE,
                  iter_control = list(tol = 1e-2, max_iter = 20, 
                                      verbose = TRUE),
                  em_control = list(tol = 1e-5, max_iter = 100),
                  lme_control = lme4::lmerControl(),
                  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                  trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                              nrow = 2, 
                                                              byrow = TRUE),
                                                       matrix(c(-1, 0, 1, -1),
                                                              nrow = 2, 
                                                              byrow = TRUE)),
                                       node = list(2, 2),
                                       solver = "ECOS",
                                       B = 100))


res_dunn = output$res_dunn

#Keep items that were differentially expressed (TRUE)
df_fig_dunn = res_dunn %>%
  dplyr::filter(diff_Treatment..TM2 == 1 | diff_Treatment..TM3 == 1 | diff_Treatment..TM4 == 1) %>%
  #Add columns with corresponding lfc value if differentially expressed, and 0 if it is not
  mutate(lfc_TM2 = ifelse(diff_Treatment..TM2 == 1, 
                                     lfc_Treatment..TM2, 0),
         lfc_TM3 = ifelse(diff_Treatment..TM3 == 1, 
                                     lfc_Treatment..TM3, 0),
         lfc_TM4 = ifelse(diff_Treatment..TM4 == 1, 
                                     lfc_Treatment..TM4, 0)) %>%
  #keep the taxon column, convert two previous columns to new columns based upon them with two sig figs
  transmute(taxon, 
            `50% DMM vs. 100% FM` = round(lfc_TM2, 2), 
            `100% DMM vs. 100% FM` = round(lfc_TM3, 2),
            `50% WMM vs. 100% FM` = round(lfc_TM4, 2)) %>%
  #Convert from wide to long converting the two columns made above to a group column and their values to a value column
  pivot_longer(cols = `50% DMM vs. 100% FM`:`50% WMM vs. 100% FM`, 
               names_to = "group", values_to = "value") %>%
  arrange(taxon)

df_fig_dunn$taxon<-gsub("Genus:","",as.character(df_fig_dunn$taxon))


#Find the miniumum value and keep the nearest whole number below it
lo = floor(min(df_fig_dunn$value))

#Find the maximum value and keep the nearest whole number above it
up = ceiling(max(df_fig_dunn$value))

#Find the value that is halfway between the two above values
mid = (lo + up)/2

#Change levels of groups columnn
df_fig_dunn$group = factor(df_fig_dunn$group, levels = c("50% DMM vs. 100% FM, 100% DMM vs. 100% FM, 50% WMM vs. 100% FM"))


fig_dunn = df_fig_dunn %>%
  ggplot(aes(x = group, y = taxon, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(group, taxon, label = value), color = "black", size = 4, fontface = "bold") +
  scale_x_discrete(limits=c('50% DMM vs. 100% FM', '100% DMM vs. 100% FM', '50% WMM vs. 100% FM')) +
  labs(x = NULL, y = NULL, title = "Log fold changes as compared to 100% FM") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x  = element_text(angle=45, hjust=1))
  
  
fig_dunn





#Keep items that were differentially expressed (TRUE)
df_fig_dunn_bars = res_dunn %>%
  dplyr::filter(diff_Treatment..TM2 == 1 | diff_Treatment..TM3 == 1 | diff_Treatment..TM4 == 1) %>%
  #Add columns with corresponding lfc value if differentially expressed, and 0 if it is not
  mutate(lfc_TM2 = ifelse(diff_Treatment..TM2 == 1, 
                          lfc_Treatment..TM2, 0),
         lfc_TM3 = ifelse(diff_Treatment..TM3 == 1, 
                          lfc_Treatment..TM3, 0),
         lfc_TM4 = ifelse(diff_Treatment..TM4 == 1, 
                          lfc_Treatment..TM4, 0)) %>%
  #keep the taxon column, convert two previous columns to new columns based upon them with two sig figs
  transmute(taxon, se_Treatment..TM2, se_Treatment..TM3, se_Treatment..TM4,
            `50% DMM vs. 100% FM` = round(lfc_TM2, 2), 
            `100% DMM vs. 100% FM` = round(lfc_TM3, 2),
            `50% WMM vs. 100% FM` = round(lfc_TM4, 2)) %>%
  #Convert from wide to long converting the two columns made above to a group column and their values to a value column
  pivot_longer(cols = `50% DMM vs. 100% FM`:`50% WMM vs. 100% FM`, 
               names_to = "group", values_to = "value") %>%
  filter(value !=0) %>%
  arrange(taxon)

df_fig_dunn_bars$se <- c(1.11, 1.13, 0.969)

df_fig_dunn_bars <- df_fig_dunn_bars %>% select(taxon, se, group, value )

df_fig_dunn_bars$taxon<-gsub("Genus:","",as.character(df_fig_dunn_bars$taxon))


fig_dunn_bars = df_fig_dunn_bars %>%
  ggplot(aes(x = taxon, y = value, fill = group)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  #Add error bars based upon standard error
  geom_errorbar(aes(ymin = value - se, ymax = value + se), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "Log fold changes as compared to 100% FM") + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1))
fig_dunn_bars
