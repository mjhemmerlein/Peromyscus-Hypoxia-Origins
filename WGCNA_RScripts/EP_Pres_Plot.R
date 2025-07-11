library(ggplot2)
library(ggrepel)
library(dplyr)

BW = read.csv("EP_WGCNA_Output/EP_BW_Pres_DE.csv")

library(ggrepel)

BW %>%
  ggplot(aes(x = Module_Size, y = Zsummary.pres, size = PopDE, fill = Module_Color)) + 
  geom_point(shape = 21, color = "black", stroke = 0.5, alpha = 0.8) +
  geom_text_repel(aes(label = PopDE), size = 4, nudge_y = 2) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 2, color = "blue", linetype = "dashed") +
  geom_hline(yintercept = 10, color = "red", linetype = "dashed") +
  scale_fill_identity(guide = "legend", name = "Module Color") +
  scale_y_continuous(limits = c(0, 100)) +
  theme_bw() +
  labs(x = "Module Size",
       y = "Zsummary Preservation")

BW %>%
  ggplot(aes(x = Module_Size, y = Zsummary.pres, size = PopDE, fill = Module_Color)) + 
  geom_point(shape = 21, color = "black", stroke = 0.5, alpha = 0.8) +
  geom_text(aes(label = PopDE), size = 4, color = "black") +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 2, color = "blue", linetype = "dashed") +
  geom_hline(yintercept = 10, color = "red", linetype = "dashed") +
  scale_fill_identity(guide = "legend", name = "Module Color") +
  scale_size_continuous(
    name = "Pop DE Genes",
    breaks = c(10, 100, 500, 1000, 2000, 2500),
    range = c(2, 20),
    trans = "sqrt"
  ) +
  scale_y_continuous(limits = c(0, 100)) +
  theme_bw() +
  labs(x = "Module Size", y = "Zsummary Preservation")


### Fishers on a few modules (use whole transcriptome)

### 3 most and least 3 conserved (correct for multiple testing)

### also try on everyone









