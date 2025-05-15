library(readxl)
library(tidyverse)

# Early pregnancy 

EP = read_xlsx("RNA_Seq_RawData/MetaData_EP.xlsx")

EP %>%
  group_by(Group) %>%
  summarise(MomUnique = n_distinct(Mom))

EP %>%
  group_by(Group) %>%
  summarise(count = n())


# Late pregnancy 

LP = read_xlsx("RNA_Seq_RawData/MetaData_LP.xlsx")
LP = na.omit(LP)

LP %>%
  group_by(Group) %>%
  summarise(MomUnique = n_distinct(Mom))

LP %>%
  group_by(Group) %>%
  summarise(count = n())
