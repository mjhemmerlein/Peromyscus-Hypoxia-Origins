library(tidyverse)
library(readxl)
library(lmerTest)
library(emmeans)
library(broom.mixed)

# Read and clean data
EP <- read_xlsx("Dream_Raw_Data/EP BW_ME Quantification.xlsx") %>%
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

# Variables to analyze
vars <- c("LZ_Area", "Prog_Sum", "Prog_LZ")

# Nest and model each variable
model_results <- vars %>%
  map_df(~{
    formula <- as.formula(paste(.x, "~ Pop*O2 + (1|Imp_Uniq)"))
    fit <- lmer(formula, data = EP)
    
    tibble(
      Variable = .x,
      Model = list(fit),
      ANOVA = list(anova(fit) %>% broom::tidy()),
      Summary = list(summary(fit)),
      Emmeans = list(emmeans(fit, ~ Pop*O2) %>% pairs() %>% as.data.frame())
    )
  })

# View ANOVA results as tidy table
model_results %>%
  select(Variable, ANOVA) %>%
  unnest(ANOVA)

model_results %>%
  select(Variable, Emmeans) %>%
  unnest(Emmeans)



# Filter and prep data
EP_split <- EP %>%
  filter(Pop %in% c("Lowland", "Highland")) %>%
  select(Pop, Mom, O2, LZ_Area, Prog_Sum, Prog_LZ, Imp_Uniq)

# Pivot longer to group by response variable
EP_long <- EP_split %>%
  pivot_longer(cols = c(LZ_Area, Prog_Sum, Prog_LZ), names_to = "Variable", values_to = "Value")

# Nest by Pop and Variable
EP_nested <- EP_long %>%
  group_by(Pop, Variable) %>%
  nest()

# Fit models and extract ANOVA + emmeans
split_models <- EP_nested %>%
  mutate(
    Model = map(data, ~ lmer(Value ~ O2 + (1|Imp_Uniq), data = .x)),
    ANOVA = map(Model, ~ anova(.x) %>% broom::tidy()),
    Emmeans = map(Model, ~ emmeans(.x, ~ O2) %>% pairs() %>% as.data.frame())
  )

split_models %>%
  select(Pop, Variable, ANOVA) %>%
  unnest(ANOVA)

split_models %>%
  select(Pop, Variable, Emmeans) %>%
  unnest(Emmeans)



