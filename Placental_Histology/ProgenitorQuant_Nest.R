library(tidyverse)
library(readxl)
library(lmerTest)
library(emmeans)
library(broom.mixed)

# Read and clean data

EP <- read_xlsx("Placental_Histology/EP BW_ME Quantification.xlsx") %>%
  select(-Qualitative_Notes) %>%
  drop_na() %>%
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

# LZ and Progenitor area ------

# Variables to analyze
vars <- c("LZ_Area", "Prog_Sum")

model_results <- vars %>%
  map_df(~{
    formula <- as.formula(paste(.x, "~ Pop*O2 + (1|Imp_Uniq)"))
    fit <- lmer(formula, data = EP)
    summ <- summary(fit)
    
    tibble(
      Variable = .x,
      ANOVA = list(anova(fit) %>% broom::tidy()),
      Fixed_Effects = list(as.data.frame(coef(summary(fit)))),
      Random_Effects = list(as.data.frame(VarCorr(fit))),
      Emmeans = list(emmeans(fit, ~ Pop*O2) %>% pairs(adjust = "tukey") %>% as.data.frame())
    )
  })


# View ANOVA results as tidy table
ANOVAresults = model_results %>%
  select(Variable, ANOVA) %>%
  unnest(ANOVA)

EMresults = model_results %>%
  select(Variable, Emmeans) %>%
  unnest(Emmeans)

Fixresults = model_results %>%
  select(Variable, Fixed_Effects) %>%
  unnest(Fixed_Effects)


# Check for Type III

LZ = lmer(LZ_Area ~ Pop*O2 + (1|Imp_Uniq), data = EP)
anova(LZ)
summary(LZ)

plot(resid(LZ))
qqnorm(resid(LZ));qqline(resid(LZ))

Prog = lmer(Prog_Sum ~ Pop*O2 + (1|Imp_Uniq), data = EP)
anova(Prog)
summary(Prog)

plot(resid(Prog))
qqnorm(resid(Prog));qqline(resid(Prog))


# Progenitor with LZ as a predictor -------

# Variables to analyze
vars <- c("Prog_Sum")

model_results <- vars %>%
  map_df(~{
    formula <- as.formula(paste(.x, "~ Pop*O2 + LZ_Area + (1|Imp_Uniq)"))
    fit <- lmer(formula, data = EP)
    summ <- summary(fit)
    
    tibble(
      Variable = .x,
      ANOVA = list(anova(fit) %>% broom::tidy()),
      Fixed_Effects = list(as.data.frame(coef(summary(fit)))),
      Random_Effects = list(as.data.frame(VarCorr(fit))),
      Emmeans = list(emmeans(fit, ~ Pop*O2) %>% pairs(adjust = "tukey") %>% as.data.frame())
    )
  })


# View ANOVA results as tidy table
ANOVAresults = model_results %>%
  select(Variable, ANOVA) %>%
  unnest(ANOVA)

EMresults = model_results %>%
  select(Variable, Emmeans) %>%
  unnest(Emmeans)

Fixresults = model_results %>%
  select(Variable, Fixed_Effects) %>%
  unnest(Fixed_Effects)

# Check Type III

Prog_LZ = lmer(Prog_Sum ~ Pop*O2 + LZ_Area + (1|Imp_Uniq), data = EP)
anova(Prog_LZ)
summary(Prog_LZ)

plot(resid(Prog_LZ))
qqnorm(resid(Prog_LZ));qqline(resid(Prog_LZ))


# Progenitor divided by LZ -------

# Variables to analyze
vars <- c("Prog_LZ")

model_results <- vars %>%
  map_df(~{
    formula <- as.formula(paste(.x, "~ Pop*O2 + (1|Imp_Uniq)"))
    fit <- lmer(formula, data = EP)
    summ <- summary(fit)
    
    tibble(
      Variable = .x,
      ANOVA = list(anova(fit) %>% broom::tidy()),
      Fixed_Effects = list(as.data.frame(coef(summary(fit)))),
      Random_Effects = list(as.data.frame(VarCorr(fit))),
      Emmeans = list(emmeans(fit, ~ Pop*O2) %>% pairs(adjust = "tukey") %>% as.data.frame())
    )
  })


# View ANOVA results as tidy table
ANOVAresults = model_results %>%
  select(Variable, ANOVA) %>%
  unnest(ANOVA)

EMresults = model_results %>%
  select(Variable, Emmeans) %>%
  unnest(Emmeans)

Fixresults = model_results %>%
  select(Variable, Fixed_Effects) %>%
  unnest(Fixed_Effects)

# Check Type III

Prog_LZ = lmer(Prog_LZ ~ Pop*O2 + (1|Imp_Uniq), data = EP)
anova(Prog_LZ)
summary(Prog_LZ)

plot(resid(Prog_LZ))
qqnorm(resid(Prog_LZ));qqline(resid(Prog_LZ))


