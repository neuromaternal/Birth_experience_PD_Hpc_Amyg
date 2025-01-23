# Clear workspace and set working directory
rm(list = ls())
setwd("~/BEATA-MADRID/hippocampus_amygdala_peripartum")
# Load required libraries
library(dplyr)
library(tidyr)
library(feather)
library(stringr)
library(reshape2)
library(fslmer)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(ggcorrplot)
library(ppcor)
library(gtsummary)
library(writexl)




## Load data
data <- feather::read_feather("data.feather")

# Define the structure and substructure names
hippo_subunit_names <- c("CA1_L", "CA1_R", "CA3_L", "CA3_R", "CA4_L", "CA4_R", "SubComplex_L", "SubComplex_R", "DentateGyrus_L", 
                         "DentateGyrus_R", "HATA_L", "HATA_R", "Fimbria_L", "Fimbria_R", "Fissure_L", "Fissure_R", "Tail_L", "Tail_R")
amyg_subunit_names <- c("Superficial_L", "Superficial_R", "CentroMedial_L", "CentroMedial_R", "LateroBasal_L", "LateroBasal_R")
subunit_names <- c(hippo_subunit_names, amyg_subunit_names)
global_names <- c("hippocampus_L", "hippocampus_R", "amygdala_L", "amygdala_R")


# ########### ########### ###########
# ########### ########### ###########
# ###########   RESULTS   ###########
# ########### ########### ###########
# ########### ########### ###########


# #############################################
# Amygdalar and hippocampal volume differences 
# #############################################
# MOTHERS VS. CONTROLS
# ---------------------
subunits_wide <- data %>%
  select(ID, Group, ICV, Age.ses.3, PostPartum.days, GestationalBirthWeeks,labor, BirthType,
         paste0(subunit_names,".ses.3"), paste0(subunit_names,".ses.4"),
         paste0(global_names,".ses.3"), paste0(global_names,".ses.4")) 
subunits_long <- melt(subunits_wide, id.vars = c("ID", "Group","ICV","Age.ses.3", "PostPartum.days", "GestationalBirthWeeks", "labor", "BirthType"),
                      variable.name = "subunits") %>%
  separate(subunits, c("subunits", "ses"), c(".ses.")) %>%
  mutate_at(vars(ses), list(~ plyr::revalue(as.factor(.), c("3" = "ses.3", "4" = "ses.4")))) %>%
  arrange(ID)

nps_long <- merge(melt(data %>% select(ID, starts_with("Total.PSQI.ses.")), id.vars = c("ID"), value.name = "PSQI") %>% 
        separate(variable, c("variable", "ses"), c("Total.PSQI.")) %>% select(-variable), 
      melt(data %>% select(ID, starts_with( "Total.PSS.14.ses.")), id.vars = c("ID"), value.name = "PSS") %>% 
        separate(variable, c("variable", "ses"), c("Total.PSS.14.")) %>% select(-variable),
      by = c("ID", "ses"))
subunits_long <- merge(subunits_long, nps_long, by = c("ID", "ses"))

source("lme_subunits_table.R")

# 1. Null Model: 
  # Volume ~ Group * Session

# Define inputs
C <- matrix(c(0, 1, 0, 0,
              0, 1, 0, 1,
              0, 0, 0, 1), nrow = 3, byrow = TRUE)
C_names <- c("Prg", "Post", "Group*Session")
X <- model.matrix(~Group*ses, subunits_long %>% filter(subunits == subunit_names[1]))

# Compute
table_1_and_S1 <- lme_subunits_table(C, C_names, X, subunits_long, c(subunit_names, global_names), descriptives = TRUE, "Group")

# Save results
table_1_and_S1 <- table_1_and_S1 %>% select(everything(),  -Prg, -Post, -`Group*Session`, Prg, Post, `Group*Session`)

write.table(table_1_and_S1 %>% 
              filter(subunits %in% c("Left Whole Hippocampus", "Right Whole Hippocampus", "Left Whole Amygdala", "Right Whole Amygdala")) %>%
              mutate(subunits = ifelse(statistic == "Mean", subunits, NA)),
             "Table1.txt", quote = F, row.names = F,  na = "", sep=";")
write.table(table_1_and_S1 %>% 
              filter(!subunits %in% c("Left Whole Hippocampus", "Right Whole Hippocampus", "Left Whole Amygdala", "Right Whole Amygdala")) %>%
              mutate(subunits = ifelse(statistic == "Mean", subunits, NA)),
                     "SMTable1.txt", quote = F, row.names = F,  na = "", sep=";")
rm(C, X)

# 2. Complete Model: 
# Volume ~ Group * Session + Age + ICV

# Define inputs
C = matrix(c(0,1,0,0,0,0,
             0,1,0,0,0,1,
             0,0,0,0,0,1),nrow=3,byrow=T)
X <- model.matrix(~Group*ses+Age.ses.3+ICV, subunits_long %>% filter(subunits == subunit_names[1]))

# Compute
table_S2 <- lme_subunits_table(C, C_names, X, subunits_long, c(subunit_names, global_names), descriptives = FALSE)

# Save results
table_S2 <- table_S2 %>% select(everything(),  -Prg, -Post, -`Group*Session`, Prg, Post, `Group*Session`)
write.table(table_S2 %>%
              mutate(subunits = ifelse(statistic == "F-stat.", subunits, NA)),
            "SMTable2.txt", quote = F, row.names = F,  na = "", sep=";")
rm(C, X)


# 3. Complete Model: 
# Volume ~ Group * Session + Age + ICV +Total.PSQI + Total.PSS

# Define inputs
C = matrix(c(0,1,0,0,0,0, 0, 0,
             0,1,0,0,0,0,0,1,
             0,0,0,0,0,0,0,1),nrow=3,byrow=T)
X <- model.matrix(~Group*ses+Age.ses.3+ICV+PSQI+PSS, subunits_long %>% filter(subunits == subunit_names[1]))

# Compute
table_S3 <- lme_subunits_table(C, C_names, X, subunits_long, c(subunit_names, global_names), descriptives = FALSE)
write.table(table_S3 %>%
              mutate(subunits = ifelse(statistic == "F-stat.", subunits, NA)),
            "SMTable_2.3.txt", quote = F, row.names = F,  na = "", sep=";")
rm(C, X)

# Save results
table_S2 <- table_S2 %>% select(everything(),  -Prg, -Post, -`Group*Session`, Prg, Post, `Group*Session`)
write.table(table_S2 %>%
              mutate(subunits = ifelse(statistic == "F-stat.", subunits, NA)),
            "SMTable2.txt", quote = F, row.names = F,  na = "", sep=";")
rm(C, X)

# PRG VS. POST in MOTHERS 
# ------------------------
subunits_long <- subunits_long %>% 
  filter( Group == "mother")
# 3. Complete Model: 
# Volume ~ Session + Age + ICV + pp-time + gestational.birth

# Define inputs
C <- matrix(c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 1, byrow = TRUE)
# # Define the design matrix
X <- model.matrix(~ ses*scale(Age.ses.3, center=T, scale=F)+scale(ICV, center=T, scale=F)*ses +
                    scale(PostPartum.days, center=T, scale=F)*ses + scale(GestationalBirthWeeks, center=T, scale=F)*ses, 
                  subunits_long %>% filter(subunits == subunit_names[1]))

sum_table <- as.data.frame(matrix(, nrow = length(c(subunit_names, global_names)), ncol = 6))
names(sum_table) <- c("subunits", "Fstatistic", "uncP", "dfWithin", "dfBetween", "sgn")
cont <- 0 
# LME statistics
for (sub in c(subunit_names, global_names)) {
  cont <- cont + 1  
  ni <- matrix(unname(table(subunits_long$ID[subunits_long$subunits == sub])), ncol = 1)
  # Extract response variable 'Y' for the given 'sub'
  Y <- subunits_long %>% filter(subunits == sub) %>% select(value) %>% as.matrix()
  # Fit the linear mixed-effects model and retrieve statistics
  stats <- lme_fit_FS(X, Zcols = c(1), Y, ni)
  # Compute F-statistics based on the contrast matrix 'C'
  F_C <- lme_F(stats, C)
  
  # Update the corresponding row in 'sum_table' with summary statistics
  sum_table[cont, ] <- c(sub, F_C$F[1], F_C$pval[1], F_C$df[1:2], F_C$sgn)
}
# make it prettier
sum_table <- sum_table %>%
  mutate_at(c("Fstatistic", "uncP", "dfWithin", "dfBetween", "sgn"), as.numeric) %>%
  mutate(
    partial_etaSq = Fstatistic * dfWithin / (Fstatistic * dfWithin + dfBetween),
    fdr = p.adjust(uncP, method = "fdr"),
    fdr_signif = case_when(
      uncP < 0.05 & fdr > 0.05 ~ "*", 
      uncP < 0.05 & fdr < 0.05 ~ "**", 
      TRUE ~ ""),
    Fstatistic = round(Fstatistic, 2),
    sgn_effect = gsub("0*$", "", format(signif(sgn * partial_etaSq, 2), scientific = FALSE)),
    uncP = paste0(gsub("0*$", "", format(signif(uncP, 2), scientific = FALSE)), fdr_signif)
    )

table_S3 <- sum_table %>%
  arrange(match(subunits, c("hippocampus_L", "hippocampus_R", hippo_subunit_names,
                            "amygdala_L", "amygdala_R", amyg_subunit_names))) %>%
  separate(subunits, c("subunits", "hemi"), c("_")) %>%
  mutate(hemi = ifelse(hemi == "L", "Left", "Right")) %>%
  mutate(subunits = case_when(subunits == "SubComplex" ~ "Subicular Complex",
                              subunits == "DentateGyrus" ~ "Dentate Gyrus",
                              subunits == "Tail" ~ "Hippocampal Tail",
                              subunits == "CentroMedial" ~ "Centromedial",
                              subunits == "LateroBasal" ~ "Laterobasal",
                              subunits == "hippocampus" ~ "Whole Hippocampus",
                              subunits == "amygdala" ~ "Whole Amygdala",
                              TRUE ~ subunits)
  ) %>%
  mutate(subunits = paste(hemi, subunits)) %>%
  select("subunits", "Fstatistic", "uncP", "dfWithin", "dfBetween", "sgn_effect" )

write.table(table_S3, "SMTable3.txt",quote = F,row.names = F, na = "", sep=";")



# 4. Complete Model: 
# Volume ~ BirthType*Session + Age*Session + ICV*Session + pp-time*Session + gestational.birth*Session

# Define inputs

C <- matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
              0, 0, 0, 0, 0, 0, 0, 0, 1, -1, 0, 0, 0, 0), nrow = 3, byrow = TRUE)
X <- model.matrix(~BirthType * ses + ses*scale(Age.ses.3, center=T, scale=F)+scale(ICV, center=T, scale=F)*ses +
                    scale(PostPartum.days, center=T, scale=F)*ses + scale(GestationalBirthWeeks, center=T, scale=F)*ses, subunits_long %>% filter(subunits == subunit_names[1]))
C_names <- c("V_E", "V_S", "E_S")

# Compute
table_S4 <- lme_subunits_table(C, C_names, X, subunits_long, c(subunit_names, global_names), descriptives = TRUE, group_var = "BirthType")

# Save results
write.table(table_S4 %>%
              mutate(subunits = ifelse(statistic == "Mean", subunits, NA)), "SMTable4.txt",quote = F,row.names = F, na = "", sep=";")


# 5. Complete Model: 
# Volume ~ Labor*Session + Age + ICV + pp-time + gestational.birth

# Define inputs
C_names <- c("Labor*Session")
C <- matrix(c(
  0, 0, 0, 0, 0, 0 , 0, 1, 0, 0, 0, 0), nrow = 1, byrow = TRUE)
X <- model.matrix(~labor*ses + ses*scale(Age.ses.3, center=T, scale=F)+scale(ICV, center=T, scale=F)*ses +
                    scale(PostPartum.days, center=T, scale=F)*ses + scale(GestationalBirthWeeks, center=T, scale=F)*ses, subunits_long %>% filter(subunits == subunit_names[1]))

# Compute
table_S5 <- lme_subunits_table(C, C_names, X, subunits_long, c(subunit_names, global_names), descriptives = TRUE, "labor")

# Save results
write.table(table_S5 %>% 
              mutate(subunits = ifelse(statistic == "Mean", subunits, NA)),
            "SMTable5.txt", quote = F, row.names = F,  na = "", sep=";")
rm(C, X)

# #############################################################
# Perinatal depression symptoms and postpartum birth experience
# #############################################################

# 1. Descriptive statistics
data %>%
  filter(Group == "mother") %>% 
  summarize(
    across(c(contains("Total.EDS."), Total.BEQ), 
           list(mean = ~mean(.), sd = ~sd(.), min = ~min(.), max = ~max(.)), 
           .names = "{.col}_{.fn}"))

data %>%
  filter(Group == "mother") %>% 
  filter(Total.EDS.ses.3 > 10) %>%
  count()

data %>%
  filter(Group == "mother") %>% 
  filter(Total.EDS.ses.4 > 10) %>%
  count()


# 3. Longitudinal changes in EDS
EDS_long <- melt(data %>%
                   filter(Group == "mother") %>% 
                   select(ID, contains("Total.EDS.")) ,
                   id.vars = c("ID"), variable.name = "ses", value.name = "Total.EDS_value") %>%
  mutate(ses = str_replace(ses, "Total.EDS.","")) %>%
  drop_na() %>%
  arrange(ID) 

# Create an empty data frame for storing summary statistics
sum_table <- as.data.frame(matrix(, nrow = 1, ncol = 5))
names(sum_table) <- c("Fstatistic", "uncP", "dfWithin", "dfBetween", "sgn")

# LME statistics
X <- model.matrix(~ ses, EDS_long)
C <- matrix(c(0, 1), nrow = 1, byrow = TRUE)
ni <- matrix(unname(table(EDS_long$ID)), ncol = 1)
Y <- EDS_long %>% select(Total.EDS_value) %>% as.matrix()
stats <- lme_fit_FS(X, Zcols = c(1), Y, ni)
F_C <- lme_F(stats, C)

sum_table[1, ] <- c(F_C$F[1], F_C$pval[1], F_C$df[1:2], F_C$sgn)
sum_table <- sum_table %>%
  mutate_at(c("Fstatistic", "uncP", "dfWithin", "dfBetween", "sgn"), as.numeric) %>%
  mutate(
    partial_etaSq = Fstatistic * dfWithin / (Fstatistic * dfWithin + dfBetween),
    Fstatistic = format(round(Fstatistic, 2), nsmall = 2),
    sgn_effect = format(round(sgn * partial_etaSq, 3), nsmall = 3),
    uncP = format(as.numeric(uncP), digits = 4, scientific = F)
  )

Fig1.A <- ggplot(EDS_long)+
  aes(x = ses, y = Total.EDS_value)+
  geom_boxplot(color = "black", fill = "lightgray", alpha = 0.7)+
  geom_point(position = position_jitter(width = 0.35), ) + 
  scale_x_discrete(labels = c("Prg", "Post")) +
  labs(y = "Depression", x = "") +
  theme_bw()

# 3. Pearson Correlation
Fig1.B <- ggscatter(data %>% filter(Group == "mother"), x = "Total.BEQ", y = "dif_Total.EDS",
          add = "reg.line", conf.int = TRUE,    
          add.params = list(fill = "lightgray"),
          ggtheme = theme_bw(),
          size=2,alpha=0.8)+
  labs(x = "Negative Birth Experience", y = "Depression (Post - Prg)")+
  stat_cor(method = "pearson", 
           label.x = 4.75, label.y = 16) 
  

Fig1 <- plot_grid(Fig1.A, Fig1.B, labels = c("A", "B"))
ggsave(Fig1, filename = "Fig1.pdf", width = 12, height = 6, units = "in", bg = "white", dpi = 600)


# #########################################################################
# hippocampal-amygdalar volumes, depression symptoms, and birth experience
# #########################################################################

# 1. WHOLE STRUCTURE
# -------------------
covariates <- c("PostPartum.days", "GestationalBirthWeeks")
quest_names <- c("Total.BEQ", "dif_Total.EDS") #dif_Total.MAS
corr_pchglobal_wide <- data %>%
  filter(Group == "mother") %>%
  select(ID, all_of(covariates), all_of(global_names), all_of(quest_names)) 

corr_pchlong_long <- melt(
  corr_pchglobal_wide, id.vars = c("ID", global_names, covariates), variable.name = "scale") %>%
  rename(scale_value = value) %>%
  melt(., id.vars = c("ID", "scale", "scale_value", covariates), variable.name = "region") %>%
  rename(region_value = value) %>%
  separate(region, c("region","hemi"),"_")

# Calculate correlation coefficients and p-values
r_WB <- corr_pchlong_long %>%
  group_by(region, scale, hemi) %>%
  summarize(
    cor = cor.test(region_value, scale_value, method = "pearson")$estimate,
    p = cor.test(region_value, scale_value, method = "pearson")$p.value
  ) %>%
  mutate(fdr.sig = ifelse(p.adjust(p, "fdr")<0.05, "*",""),
         label = paste0("R= ", round(cor, 3), "\nP= ", round(p, 3), fdr.sig))


col_rh_h <- "#F09E6A"; col_lh_h <- "#E9C46A";col_lh_a <- "#49BF9B"; col_rh_a <- "#5CA9C8"
source ("plot_region_corr.R")

legend_values <- r_WB$label[r_WB$scale == "dif_Total.EDS" & r_WB$region == "amygdala"]
p_a_eds = plot_region_corr(data = corr_pchlong_long,
                      scale_name = "dif_Total.EDS",
                      region_name = "amygdala",
                      legend_values = legend_values,
                      scale_name_text = "Depression (Post - Prg)",
                      col_lh = col_lh_a,
                      col_rh = col_rh_a)

legend_values <- r_WB$label[r_WB$scale == "Total.BEQ" & r_WB$region == "hippocampus"]
p_h_beq = plot_region_corr(data = corr_pchlong_long,
                      scale_name = "Total.BEQ",
                      region_name = "hippocampus",
                      legend_values = legend_values,
                      scale_name_text = "Negative Birth Experience",
                      col_lh = col_lh_h,
                      col_rh = col_rh_h)

legend_values <- r_WB$label[r_WB$scale == "Total.BEQ" & r_WB$region == "amygdala"]
p_a_beq = plot_region_corr(data = corr_pchlong_long,
                      scale_name = "Total.BEQ",
                      region_name = "amygdala",
                      legend_values = legend_values,
                      scale_name_text = "Negative Birth Experience",
                      col_lh = col_lh_a,
                      col_rh = col_rh_a)

legend_values <- r_WB$label[r_WB$scale == "dif_Total.EDS" & r_WB$region == "hippocampus"]
p_h_eds = plot_region_corr(data = corr_pchlong_long,
                      scale_name = "dif_Total.EDS",
                      region_name = "hippocampus",
                      legend_values = legend_values,
                      scale_name_text = "Depression (Post - Prg)",
                      col_lh = col_lh_h,
                      col_rh = col_rh_h)


# Partial correlations
source("compute_partial_cor.R")

table_S6 <- data.frame()
for (region in c("hippocampus_L", "hippocampus_R", "amygdala_L", "amygdala_R")) {
  for (quest in quest_names) {
    result <- compute_partial_cor(corr_pchglobal_wide%>%select(region),
                                  corr_pchglobal_wide%>%select(quest),
                                  corr_pchglobal_wide%>%select(covariates),
                                  corr_method = "spearman")
    table_S6 <- rbind(table_S6, 
                     data.frame(Structure = region, Questionnaire = quest, partial_correlation = result$cor, p_value = result$p.value) 
                     )
  }
}
table_S6 <- table_S6 %>%
  separate(Structure, c("Structure", "hemi"), c("_")) %>%
  mutate(hemi = case_when(hemi == "L" ~ "lh", hemi == "R" ~ "rh")) %>%
  mutate(fdr = p.adjust(p_value, method = "fdr"),
         fdr_signif = case_when(
           p_value < 0.05 & fdr > 0.05 ~ "*", 
           p_value < 0.05 & fdr < 0.05 ~ "**", 
           TRUE ~ ""),
         partial_correlation = signif(partial_correlation, 3),
         p_value = paste0(signif(p_value, 3), fdr_signif),
         Questionnaire = case_when(Questionnaire == "Total.BEQ" ~ "Birth Experience", Questionnaire == "dif_Total.EDS" ~ "Depression (Post - Prg)")
  ) %>% 
  select(-starts_with("fdr")) %>%
  melt(., id.vars = c("Structure", "Questionnaire", "hemi"),
               measure.vars = c("partial_correlation", "p_value")) %>%
  dcast(., Structure + Questionnaire ~ hemi + variable, 
                    value.var = "value") %>%
  select(Structure, Questionnaire, starts_with("lh"), starts_with("rh")) %>%
  rename(`Psychological Scale` = Questionnaire)

write.table(table_S6, "SMTable6.txt", quote = F, row.names = F,  na = "", sep=";")

# 2. POST-HOC SUBSTRUCTURES
# ----------------

# corr matrix
corr_pch_long <- melt( data %>% select(ID, all_of(subunit_names), Total.BEQ, dif_Total.EDS, PostPartum.days, GestationalBirthWeeks), 
                       id.vars = c("ID", subunit_names, "PostPartum.days", "GestationalBirthWeeks"), variable.name = "scale", value.name = "scale_value") %>%
  melt(., id.vars = c("ID", "scale", "scale_value","PostPartum.days"), variable.name = "region", value.name = "region_value") %>%
  separate(region, c("region","hemi"), c("_"))

correlation_subunits <- corr_pch_long %>%
  group_by(region, scale, hemi) %>%
  summarize(cor = cor.test(region_value, scale_value, method = "pearson")$estimate,
            p = cor.test(region_value, scale_value, method = "pearson")$p.value) %>%
  mutate(region = paste0(region, "_", hemi)) %>%
  arrange(match(region, subunit_names))

r_M = dcast(correlation_subunits %>% select(scale, region, cor), scale ~ region, value.var = "cor") %>%
  tibble::column_to_rownames(var = "scale") %>%
  select(all_of(hippo_subunit_names), all_of(amyg_subunit_names))
p_M = dcast(correlation_subunits %>% select(scale, region, p), scale ~ region, value.var = "p") %>%
  tibble::column_to_rownames(var = "scale") %>%
  select(all_of(hippo_subunit_names), all_of(amyg_subunit_names))  

row.names(r_M) <- c(" ", "")#c("Depression EDS (Post-Prg)","Birth Experience BEQ (Post)")
row.names(p_M) <- row.names(r_M)
colnames(r_M) <- str_replace(colnames(r_M), pattern = "_L", replacement = "")
colnames(r_M) <- str_replace(colnames(r_M), pattern = "_R", replacement = "")
colnames(p_M) <- colnames(r_M)

r_M_hpc <- rbind(r_M[2,seq(1,length(hippo_subunit_names),by=2)],
                 r_M[2,seq(2,length(hippo_subunit_names),by=2)])
p_M_hpc <- rbind(p_M[2,seq(1,length(hippo_subunit_names),by=2)],
                 p_M[2,seq(2,length(hippo_subunit_names),by=2)])
rownames(r_M_hpc) <- c("lh","rh")
rownames(p_M_hpc) <- rownames(r_M_hpc) 

corrplot_hpc <- ggcorrplot(t(r_M_hpc), p.mat = t(p_M_hpc),  insig = "blank") +
  ggtitle("Negative Birth Experience")+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.text.x = element_text(color = "#3A3A3A",size=10),
        axis.text.y = element_text(color = "#3A3A3A",size=10),
        legend.key.size = unit(0.75, "lines"),
        legend.title = element_text(angle = 0, size=8),
        legend.text = element_text(size=6.5)) +
  guides(fill=guide_legend(title='R', title.position = "top"))

r_M_amg <- r_M[1, seq(length(hippo_subunit_names)+2,ncol(r_M), by=2)]
p_M_amg <- p_M[1, seq(length(hippo_subunit_names)+2,ncol(r_M), by=2)]
rownames(r_M_amg) <- c("rh")
rownames(p_M_amg) <-rownames(r_M_amg) 

corrplot_amg <- ggcorrplot(t(r_M_amg), p.mat = t(p_M_amg),  insig = "blank") + 
  ggtitle("Depression (Post - Prg)")+
  theme(plot.title = element_text(size = 10,hjust = 0.5),
        axis.text.x = element_text(color = "#3A3A3A",size=10),
        axis.text.y = element_text(color = "#3A3A3A",size=10),
        legend.key.size = unit(0.75, "lines"),
        legend.title = element_text(angle = 0, size=8),
        legend.text = element_text(size=6.5)) +
  guides(fill=guide_legend(title='R', title.position = "top"))

# SAVE FIGURES
Fig2 <- plot_grid(
  plot_grid(p_h_beq, corrplot_hpc+ theme(plot.margin = unit(c(1,0,2,1), "lines")), 
            ncol = 1, labels = c("A", "C"), 
            rel_heights = c(1, 0.5)),
  plot_grid(p_a_eds, corrplot_amg+ theme(plot.margin = unit(c(1,0,3.5,1), "lines")),
            ncol=1, labels = c("B", "D"),
            rel_heights = c(1, 0.5)),
  ncol=2, align = "v")
ggsave("Fig2.pdf", Fig2, width = 10, height = 7.5, units = "in", bg = "white", dpi = 600)

Fig_S1 <- plot_grid(p_h_eds, p_a_beq, ncol=2, align = "v", labels = c("A", "B"))
ggsave("FigSM1.tiff", Fig_S1, width = 10, height = 5, units = "in", bg = "white", dpi = 600)


r_subunits <- corr_pch_long %>%
  group_by(region, scale, hemi) %>%
  summarize(
    cor = cor.test(region_value, scale_value, method = "pearson")$estimate,
    p = cor.test(region_value, scale_value, method = "pearson")$p.value
  ) %>%
  mutate(label = paste0("R= ", round(cor, 4), 
                        "\np= ", format(p, digits = 2, scientific = F)) ) 

p_hippo <- list()
for (sub in unique(sapply(str_split(hippo_subunit_names, "_"), function(x) x[[1]]))){
  legend_values <- r_subunits$label[r_subunits$scale == "Total.BEQ" & r_subunits$region == sub]
  p_hippo[[length(p_hippo) + 1]] <- plot_region_corr(data = corr_pch_long,
                                                     scale_name = "Total.BEQ",
                                                     region_name = sub,
                                                     legend_values = legend_values,
                                                     scale_name_text = "Negative Birth Experience",
                                                     col_lh = col_lh_h,
                                                     col_rh = col_rh_h)
}
FigSM2 <- plot_grid(plotlist = p_hippo, ncol=3)
ggsave("FigSM2.tiff", FigSM2, width = 14, height = 13, units = "in")

p_amyg <- list()
for (sub in unique(sapply(str_split(amyg_subunit_names, "_"), function(x) x[[1]]))){
  legend_values <- r_subunits$label[r_subunits$scale == "dif_Total.EDS" & r_subunits$region == sub]
  p_amyg[[length(p_amyg) + 1]] <- plot_region_corr(data = corr_pch_long,
                                                   scale_name = "dif_Total.EDS",
                                                   region_name = sub,
                                                   legend_values = legend_values,
                                                   scale_name_text = "Depression (Post - Prg)",
                                                   col_lh = col_lh_a,
                                                   col_rh = col_rh_a)
}
FigSM3 <- plot_grid(plotlist = p_amyg, ncol=3, labels = unique(sapply(str_split(amyg_subunit_names, "_"), function(x) x[[1]])))
ggsave("FigSM3.tiff", FigSM3, width = 14, height = 5, units = "in", bg = "white", dpi = 600)


# ########### ########### ###########
# ########### ########### ###########
# ########### DEMOGRAPHIC ###########
# ########### ########### ###########
# ########### ########### ###########
demos <- data %>%
  select(Group,Age.ses.3,InterSessions.days, Education)%>%
  tbl_summary(by=Group,
              type = list(c(Age.ses.3,InterSessions.days) ~ "continuous2"),
              statistic = list(all_continuous() ~ c("{mean} ({sd})","{min}-{max}")),
              digits = all_continuous() ~ 2,
              label = list(Age.ses.3~"Age at pregnancy session [years]",
                           InterSessions.days~"Time between sessions [days]"),
              missing = "no") 

t.test(Age.ses.3~Group, data_wide.all,var.equal=T)
t.test(InterSessions.days~Group, data_wide.all,var.equal=T)
t.test(Total.Digits.WAIS.ses.3~Group, data_wide.all,var.equal=T)
chisq.test(table(data_wide.all$Education,data_wide.all$Group))


demos <- data %>%
  filter(Group == "mother") %>%
  select(Age.ses.3,PrgWeeks,GestationalBirthWeeks,BirthType, PostPartum.days,
         Education, BirthPROM, BirthInduction,InstrumentalBirth,Episiotomy,Kristeller,SkinToSkin,FeedingType, DeliveryEscort) %>%
  tbl_summary(
    type = list(c(Age.ses.3,PrgWeeks,GestationalBirthWeeks,PostPartum.days,InterSessions.days) ~ "continuous2",
                c( Education, BirthPROM, BirthInduction,InstrumentalBirth,Episiotomy,Kristeller,SkinToSkin,FeedingType, DeliveryEscort)~"categorical"),
    statistic = list(all_continuous() ~ c("{mean} ({sd})","{min}-{max}")),
    digits = all_continuous() ~ 2,
    missing = "no"  )  
  
data%>%
  select(Conception, AssistedConception, BirthType,Baby.n,sexbaby, NeonatalCare,BloodLoss)%>%
  tbl_summary()


