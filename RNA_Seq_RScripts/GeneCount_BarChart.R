# Libraries

library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggbreak)
library(readxl)


# Read and pivot
GeneCounts = read_xlsx("RNA_Seq_RawData/DreamCounts.xlsx")

GeneCounts = GeneCounts %>%
  pivot_longer("Pop":"Highland", names_to = "DE", values_to = "DE_Freq") %>%
  mutate(
    DE = fct_relevel(DE, "Pop", "Trt", "Ixn", "Lowland", "Highland"),
    Category = fct_recode(Category,
                          "EP" = "Early Pregnancy",
                          "LP" = "Late Pregnancy",
                          "Shared" = "Conserved Across Pregnancy"),
    Category = fct_relevel(Category, "EP", "LP", "Shared"),
    )

GeneCounts %>%
  ggplot(aes(x = Category, y = DE_Freq)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ DE, scales = "free_y") +
  scale_fill_identity() +
  scale_y_log10(limits = c(1, 40000)) +  # Adjust limits to your data
  geom_text(aes(label = DE_Freq), vjust = -0.5, size = 8) +
  labs(
    y = expression(bold("Differential Gene Frequency")),
    x = NULL
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title.y = element_text(face = "bold", size = 24, color = "black"),
    axis.text.y = element_text(face = "bold", size = 20, color = "black"),
    axis.text.x = element_text(face = "bold", size = 20, color = "black"),
    plot.title = element_text(face = "bold", size = 24, hjust = 0.5)
  )

ggsave("Plots/GeneCount_BarChart/BarChart.pdf", width = 12, height = 8, units = "in", dpi = 300)


# Read in data

GeneCounts = read_xlsx("RNA_Seq_RawData/DreamCounts.xlsx")

GeneCounts = GeneCounts %>%
  pivot_longer("Pop":"Highland", names_to = "DE", values_to = "DE_Freq")

GeneCounts$DE = fct_relevel(GeneCounts$DE, "Pop", "Trt", "Ixn", "Lowland", "Highland")

GeneCounts$Category = fct_relevel(GeneCounts$Category, "Early Pregnancy", "Late Pregnancy", "Conserved Across Pregnancy")

PopOnly = c("Lowland", "Highland")

GeneCounts = GeneCounts %>%
  mutate(color = case_when(
    DE == "Lowland" ~ "Low", 
    DE == "Highland" ~ "High", 
    TRUE ~ "NA"))


GeneCounts %>%
  ggplot(aes(x = Category, y = log(DE_Freq), fill = color)) +
  facet_wrap(DE) + 
  scale_fill_manual(values = c("skyblue", "goldenrod1", "gray")) +
  geom_bar(stat = "identity") + 
  theme_bw() +
  scale_y_continuous(limits = c(0, 10.5)) +
  geom_text(aes(label = DE_Freq), vjust = -0.5, size = 8) +
  labs(title = "Early Pregnancy - 21,545 Genes", y = expression(bold(atop("Differential Gene Frequency", paste("Log(GeneCounts)")))), x = NULL) +
  theme(
    legend.position = "none",
    axis.title.y = element_text(face = "bold", size = 24, color = "black"),
    axis.text.y = element_text(face = "bold", size = 20, color = "black"),
    axis.text.x = element_text(face = "bold", size = 20, color = "black"),
    plot.title = element_text(face = "bold", size = 24, hjust = 0.5))


# Log Axis

# Plot
GeneCounts %>%
  filter(Category == "Early Pregnancy") %>%
  ggplot(aes(x = DE, y = log(DE_Freq), fill = color)) +
  scale_fill_manual(values = c("skyblue", "goldenrod1", "gray")) +
  geom_bar(stat = "identity") + 
  theme_bw() +
  scale_y_continuous(limits = c(0, 10.5)) +
  geom_text(aes(label = DE_Freq), vjust = -0.5, size = 8) +
  labs(title = "Early Pregnancy - 21,545 Genes", y = expression(bold(atop("Differential Gene Frequency", paste("Log(GeneCounts)")))), x = NULL) +
  theme(
    legend.position = "none",
    axis.title.y = element_text(face = "bold", size = 24, color = "black"),
    axis.text.y = element_text(face = "bold", size = 20, color = "black"),
    axis.text.x = element_text(face = "bold", size = 20, color = "black"),
    plot.title = element_text(face = "bold", size = 24, hjust = 0.5))

ggsave("Dream_Output/Hypoxia_Figs/EP_GeneCount.png", width = 8, height = 6, units = "in", dpi = 300)


GeneCounts %>%
  filter(Category == "Late Pregnancy") %>%
  ggplot(aes(x = DE, y = log(DE_Freq), fill = color)) +
  scale_fill_manual(values = c("skyblue", "goldenrod1", "gray")) +
  geom_bar(stat = "identity") + 
  theme_bw() +
  scale_y_continuous(limits = c(0, 10.5)) +
  geom_text(aes(label = DE_Freq), vjust = -0.5, size = 8) +
  labs(title = "Late Pregnancy - 22,370 Genes", y = expression(bold(atop("Differential Gene Frequency", paste("Log(GeneCounts)")))), x = NULL) +
  theme(
    legend.position = "none",
    axis.title.y = element_text(face = "bold", size = 24, color = "black"),
    axis.text.y = element_text(face = "bold", size = 20, color = "black"),
    axis.text.x = element_text(face = "bold", size = 20, color = "black"),
    plot.title = element_text(face = "bold", size = 24, hjust = 0.5))

ggsave("Dream_Output/Hypoxia_Figs/LP_GeneCount.png", width = 8, height = 6, units = "in", dpi = 300)


GeneCounts %>%
  filter(Category == "Conserved Across Pregnancy") %>%
  ggplot(aes(x = DE, y = log(DE_Freq), fill = color)) +
  scale_fill_manual(values = c("skyblue", "goldenrod1", "gray")) +
  geom_bar(stat = "identity") + 
  theme_bw() +
  scale_y_continuous(limits = c(0, 10.5)) +
  geom_text(aes(label = DE_Freq), vjust = -0.5, size = 8) +
  labs(title = "Conserved Across Pregnancy", y = expression(bold(atop("Differential Gene Frequency", paste("Log(GeneCounts)")))), x = NULL) +
  theme(
    legend.position = "none",
    axis.title.y = element_text(face = "bold", size = 24, color = "black"),
    axis.text.y = element_text(face = "bold", size = 20, color = "black"),
    axis.text.x = element_text(face = "bold", size = 20, color = "black"),
    plot.title = element_text(face = "bold", size = 24, hjust = 0.5))

ggsave("Dream_Output/Hypoxia_Figs/Cons_GeneCount.png", width = 8, height = 6, units = "in", dpi = 300)



# No log axis
# Plot
GeneCounts %>%
  filter(Category == "Early Pregnancy") %>%
  ggplot(aes(x = DE, y = DE_Freq)) +
  geom_bar(stat = "identity") + 
  theme_bw() +
  scale_y_continuous(limits = c(0, 13500),
                     breaks = c(seq(0, 8000, by = 2500), 9000, 11000, 13000),
                     sec.axis = sec_axis(trans = ~., breaks = NULL)) +
  scale_y_break(c(9000, 12000)) +
  geom_text(aes(label = DE_Freq), vjust = -0.5, size = 7) +
  labs(title = "Early Pregnancy", y = "Differential Gene Freqeuency", x = NULL) +
  theme(
    axis.title.y = element_text(face = "bold", size = 20),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 24, hjust = 0.5))

#ggsave("Dream_Output/Hypoxia_Figs/EP_GeneCount.png", width = 8, height = 6, units = "in", dpi = 300)


GeneCounts %>%
  filter(Category == "Late Pregnancy") %>%
  ggplot(aes(x = DE, y = DE_Freq)) +
  geom_bar(stat = "identity") + 
  theme_bw() +
  scale_y_continuous(limits = c(0, 13500),
                     breaks = c(seq(0, 8000, by = 2500), 9000, 11000, 13000),
                     sec.axis = sec_axis(trans = ~., breaks = NULL)) +
  scale_y_break(c(9000, 12000)) +
  geom_text(aes(label = DE_Freq), vjust = -0.5, size = 7) +
  labs(title = "Late Pregnancy", y = "Differential Gene Freqeuency", x = NULL) +
  theme(
    axis.title.y = element_text(face = "bold", size = 20),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 24, hjust = 0.5))

#ggsave("Dream_Output/Hypoxia_Figs/LP_GeneCount.png", width = 8, height = 6, units = "in", dpi = 300)



GeneCounts %>%
  filter(Category == "Conserved Across Pregnancy") %>%
  ggplot(aes(x = DE, y = DE_Freq)) +
  geom_bar(stat = "identity") + 
  theme_bw() +
  scale_y_continuous(limits = c(0, 13500),
                     breaks = c(seq(0, 8000, by = 2500), 9000, 11000, 13000),
                     sec.axis = sec_axis(trans = ~., breaks = NULL)) +
  scale_y_break(c(9000, 12000)) +
  geom_text(aes(label = DE_Freq), vjust = -0.5, size = 7) +
  labs(title = "Conserved Across Pregnancy", y = "Differential Gene Freqeuency", x = NULL) +
  theme(
    axis.title.y = element_text(face = "bold", size = 20),
    axis.text.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(face = "bold", size = 12),
    plot.title = element_text(face = "bold", size = 24, hjust = 0.5))

#ggsave("Dream_Output/Hypoxia_Figs/Cons_GeneCount.png", width = 8, height = 6, units = "in", dpi = 300)


