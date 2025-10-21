
library(tidyverse)
library(ggplot2)
library(readxl)
library(emmeans) 
library(lmerTest)

# Read in data

EP = read_xlsx("Placental_Histology/EP BW_ME Quantification.xlsx")

EP = EP %>%
  select(-Qualitative_Notes) %>%
  na.omit(EP)

EP$Mom = as.factor(EP$Mom)
EP$Pop = as.factor(EP$Pop)
EP$O2 = as.factor(EP$O2)

EP$LZ_Area = as.numeric(EP$LZ_Area)

EP = EP %>%
  mutate(Prog_Sum = Progenitor_Area_1 + Progenitor_Area_2)

EP = EP %>%
  mutate(Prog_LZ = Prog_Sum/LZ_Area) %>%
  mutate(Pop = fct_recode(Pop, 
                          "Lowland" = "BW", 
                          "Highland" = "ME")) %>%
  mutate(Group = paste(Pop,O2)) 

EP$Group = as.factor(EP$Group)

EP = EP %>%
  mutate(Group = factor(Group, levels = c("Lowland 1N", "Lowland 2H", "Highland 1N", "Highland 2H")))

EP_BW = EP %>%
  filter(Pop == "Lowland")

EP_ME = EP %>%
  filter(Pop == "Highland")

# LZ

EP %>%
  group_by(Group) %>%
  ggplot(aes(Group, LZ_Area, fill = Pop)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 3.5, alpha = 0.5, width = 0.1) +
  scale_fill_manual(values = c("Lowland" = "#F4C552", "Highland" = "#247BB1")) +
  theme_classic() +
  labs(y = expression(bold("Labyrinth Zone Area (mm"^2*")")), x = NULL) +
  theme(legend.position = "none") +
  theme(
    axis.title.y = element_text(face = "bold", size = 24, color = "black", margin = margin(r = 15)),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.title.x = element_text(face = "bold", size = 20),
    axis.text.y = element_text(face = "bold", size = 20, color = "black"),
    plot.title = element_text(face = "bold", size = 24, hjust = 0.5), 
    legend.text = element_text(face = "bold", size = 20),
    legend.title = element_blank())

EP %>%
  group_by(Group) %>%
  ggplot(aes(Group, LZ_Area, fill = Pop)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 6.5, width = 0.1, alpha = 0.7) +
  scale_fill_manual(values = c("Lowland" = "#F4C552", "Highland" = "#247BB1")) +
  theme_classic() +
  labs(y = expression(bold("Labyrinth Zone Area (mm"^2*")")), x = NULL) +
  ylim(0.5, 5.4) +
  theme(legend.position = "none") +
  theme(
    axis.title.y = element_text(face = "bold", size = 24, color = "black", margin = margin(r = 15)),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(face = "bold", size = 20, color = "black"),
    plot.title = element_text(face = "bold", size = 24, hjust = 0.5), 
    legend.text = element_blank(),
    legend.title = element_blank())


 ggsave("Plots/Placental_Histology/LZ_Quant.pdf", width = 7, height = 6, units = "in", dpi = 300)


# Progenitor area

EP %>%
  group_by(Group) %>%
  ggplot(aes(Group, Prog_Sum, fill = Pop)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 3.5, alpha = 0.5, width = 0.1) +
  scale_fill_manual(values = c("Lowland" = "#F4C552", "Highland" = "#247BB1")) +
  theme_classic() +
  labs(y = expression(bold("Labyrinth Zone Area (mm"^2*")")), x = NULL) +
  theme(legend.position = "none") +
  theme(
    axis.title.y = element_text(face = "bold", size = 24, color = "black", margin = margin(r = 15)),
    axis.text.x = element_text(face = "bold", size = 10),
    axis.title.x = element_text(face = "bold", size = 20),
    axis.text.y = element_text(face = "bold", size = 20, color = "black"),
    plot.title = element_text(face = "bold", size = 24, hjust = 0.5), 
    legend.text = element_text(face = "bold", size = 20),
    legend.title = element_blank())

EP %>%
  group_by(Group) %>%
  ggplot(aes(Group, Prog_Sum, fill = Pop)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_manual(values = c("Lowland" = "#F4C552", "Highland" = "#247BB1")) +
  theme_classic() +
  geom_jitter(size = 4, width = 0.1, alpha = 0.7) +
  theme(legend.position = "none") +
  labs(y = expression(bold("Fetal Progenitor Area (mm"^2*")")), x = NULL) +
  ylim(0.15, 1.8)+
  theme(
    axis.title.y = element_text(face = "bold", size = 24, color = "black", margin = margin(r = 15)),
    axis.title.x = element_blank(),
    axis.text.y = element_text(face = "bold", size = 20, color = "black"),
    axis.text.x = element_blank(),
    plot.title = element_text(face = "bold", size = 24, hjust = 0.5), 
    legend.text = element_blank(),
    legend.title = element_blank(),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 35))

 ggsave("Plots/Placental_Histology/Prog_Quant.pdf", width = 7, height = 6, units = "in", dpi = 300)

 

# Progenitor / LZ
 EP %>%
   group_by(Group) %>%
   ggplot(aes(Group, Prog_LZ, fill = Pop)) +
   geom_boxplot(outlier.shape = NA) +
   geom_jitter(size = 3.5, alpha = 0.5, width = 0.1) +
   scale_fill_manual(values = c("Lowland" = "#F4C552", "Highland" = "#247BB1")) +
   theme_classic() +
   labs(y = expression(bold("Relative Progenitor Area \n (% of Labyrinth Zone)")), x = NULL) +
   theme(legend.position = "none") +
   theme(
     axis.title.y = element_text(face = "bold", size = 24, color = "black", margin = margin(r = 15)),
     axis.text.x = element_text(face = "bold", size = 10),
     axis.title.x = element_text(face = "bold", size = 20),
     axis.text.y = element_text(face = "bold", size = 20, color = "black"),
     plot.title = element_text(face = "bold", size = 24, hjust = 0.5), 
     legend.text = element_text(face = "bold", size = 20),
     legend.title = element_blank())
 
 EP %>%
   group_by(Group) %>%
   ggplot(aes(Group, Prog_LZ, fill = Pop)) +
   geom_boxplot(outlier.shape = NA) +
   scale_fill_manual(values = c("Lowland" = "#F4C552", "Highland" = "#247BB1")) +
   theme_classic() +
   geom_jitter(size = 6.5, width = 0.1, alpha = 0.7) +
   theme(legend.position = "none") +
   labs(y = expression(bold("Relative Progenitor Area \n (% of Labyrinth Zone)")), x = NULL) +
   ylim(0.01, 0.7) +
   theme(
     axis.title.y = element_text(face = "bold", size = 24, color = "black", margin = margin(r = 15)),
     axis.title.x = element_blank(),
     axis.text.y = element_text(face = "bold", size = 20, color = "black"),
     axis.text.x = element_blank(),
     plot.title = element_text(face = "bold", size = 24, hjust = 0.5), 
     legend.text = element_blank(),
     legend.title = element_blank(),
     plot.margin = margin(t = 10, r = 10, b = 10, l = 35))

 ggsave("Plots/Placental_Histology/ProgLZ_Quant.pdf", width = 7, height = 6, units = "in", dpi = 300)






