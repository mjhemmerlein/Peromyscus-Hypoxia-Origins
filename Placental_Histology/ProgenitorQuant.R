library(tidyverse)
library(readxl)
library(lmerTest)
library(emmeans)
library(broom.mixed)
library(performance)

# Read and clean data

EP <- read_xlsx("Placental_Histology/EP BW_ME Quantification.xlsx") %>%
  select(-Qualitative_Notes) %>%
  mutate(across(c(Mom, Pop, O2), as.factor),
         LZ_Area = as.numeric(LZ_Area),
         Prog_Sum = Progenitor_Area_1 + Progenitor_Area_2,
         Prog_LZ = Prog_Sum / LZ_Area,
         Pop = fct_recode(Pop, "Lowland" = "BW", "Highland" = "ME"),
         Group = paste(Pop, O2),
         Group = factor(Group, levels = c("Lowland 1N", "Lowland 2H", "Highland 1N", "Highland 2H")),
         Imp_Number = substr(`Slide ID`, 1, 3),
         Imp_Uniq = paste0(Mom,"_",Imp_Number),
         Imp_Uniq = as.factor(Imp_Uniq))


# write.csv(EP, "Placental_Histology/Histology_Data.csv", row.names = FALSE)

# LZ area -------

LZ = lmer(LZ_Area ~ Pop*O2 + (1|Imp_Uniq), data = EP)
anova(LZ)
summary(LZ)

contrast(emmeans(LZ, ~ Pop*O2), method = "pairwise")

plot(resid(LZ))
qqnorm(resid(LZ));qqline(resid(LZ))

LZ_anova = anova(LZ)
LZ_anova = as.data.frame(LZ_anova)
LZ_anova = rownames_to_column(LZ_anova, var = "Term")

# Prog area ---------

Prog = lmer(Prog_Sum ~ Pop*O2 + (1|Imp_Uniq), data = EP)
anova(Prog)
summary(Prog)

contrast(emmeans(Prog, ~ Pop*O2), method = "pairwise")

contrast(emmeans(Prog, ~ Pop*O2), method = "pairwise", adjust = "none")

plot(resid(Prog))
qqnorm(resid(Prog));qqline(resid(Prog))

Prog_anova = anova(Prog)
Prog_anova = as.data.frame(Prog_anova)
Prog_anova = rownames_to_column(Prog_anova, var = "Term")

# Prog area log ---------

Prog_log = lmer(log(Prog_Sum) ~ Pop*O2 + (1|Imp_Uniq), data = EP)
anova(Prog_log)

contrast(emmeans(Prog_log, ~ Pop*O2), method = "pairwise")

contrast(emmeans(Prog_log, ~ Pop*O2), method = "pairwise", adjust = "none")

plot(resid(Prog_log))
qqnorm(resid(Prog_log));qqline(resid(Prog_log))

Prog_log_anova = anova(Prog_log)
Prog_log_anova = as.data.frame(Prog_log_anova)
Prog_log_anova = rownames_to_column(Prog_log_anova, var = "Term")


# Calculate proportion and model on proportion (Ratio) ----------
Prog_ratio = lmer(log(Prog_LZ) ~ Pop*O2 + (1|Imp_Uniq), data = EP)
anova(Prog_ratio)

contrast(emmeans(Prog_ratio, ~ Pop*O2), method = "pairwise")

contrast(emmeans(Prog_ratio, ~ Pop*O2), method = "pairwise", adjust = "none")

Prog_ratio_anova = anova(Prog_ratio)
Prog_ratio_anova = as.data.frame(Prog_ratio_anova)
Prog_ratio_anova = rownames_to_column(Prog_ratio_anova, var = "Term")


# Progenitor with LZ as a predictor -------

Prog_LZ = lmer(log(Prog_Sum) ~ Pop*O2 + LZ_Area + (1|Imp_Uniq), data = EP)
anova(Prog_LZ)
summary(Prog_LZ)

contrast(emmeans(Prog_LZ, ~ Pop*O2), method = "pairwise")

plot(resid(Prog_LZ))
qqnorm(resid(Prog_LZ));qqline(resid(Prog_LZ))

Prog_LZ_anova = anova(Prog_LZ)
Prog_LZ_anova = as.data.frame(Prog_LZ_anova)
Prog_LZ_anova = rownames_to_column(Prog_LZ_anova, var = "Term")



# Combine tables -------

# Add model/test labels
LZ_anova$Test        <- "LZ_area"
Prog_anova$Test      <- "Prog_area"
Prog_log_anova$Test  <- "Prog_log_area"
Prog_ratio_anova$Test <- "Prog_ratio"
Prog_LZ_anova$Test   <- "Prog_LZ_model"

# Bind all ANOVA outputs
all_anova <- dplyr::bind_rows(
  LZ_anova,
  Prog_anova,
  Prog_log_anova,
  Prog_ratio_anova,
  Prog_LZ_anova
)

all_anova

