# Coverage visualization
## 02.09.2024
### Alen Suljič

# load the libraries
library(tidyverse)

# set working directory
path <- "/path/to/coverage.csv"

setwd(path)

# load the data

c <- read_csv("coverage.csv")

# transform the data
c1 <- c %>% 
  rename(samp = 1,
         pos = 2,
         cov = 3) %>% 
  group_by(pos) %>% 
  mutate(pos_median_cov = median(cov)) %>% 
  ungroup() %>% 
  group_by(samp) %>% 
  mutate(sample_median_cov = median(cov)) %>% 
  ungroup() %>% 
  drop_na() %>% 
  mutate(global_cov = median(cov),
         gene = case_when(pos <= 265 ~ "5'UTR",
                          pos >= 266 & pos <= 21555 ~ "ORF1ab",
                          pos >= 21556 & pos <= 21562 ~ "ITR1",
                          pos >= 21563 & pos <= 25384 ~ "S",
                          pos >= 25385 & pos <= 25392 ~ "ITR2",
                          pos >= 25393 & pos <= 26220 ~ "ORF3a",
                          pos >= 26221 & pos <= 26244 ~ "ITR3",
                          pos >= 26245 & pos <= 26472 ~ "E",
                          pos >= 26473 & pos <= 26522 ~ "ITR4",
                          pos >= 26523 & pos <= 27191 ~ "M",
                          pos >= 27192 & pos <= 27201 ~ "ITR5",
                          pos >= 27202 & pos <= 27387 ~ "ORF6a",
                          pos >= 27388 & pos <= 27393 ~ "ITR6",
                          pos >= 27394 & pos <= 27887 ~ "ORF7ab",
                          pos >= 27888 & pos <= 27893 ~ "ITR7",
                          pos >= 27894 & pos <= 28259 ~ "ORF8",
                          pos >= 28260 & pos <= 28273 ~ "ITR8",
                          pos >= 28274 & pos <= 29533 ~ "N",
                          pos >= 29534 & pos <= 29557 ~ "ITR9",
                          pos >= 29558 & pos <= 29674 ~ "ORF10",
                          pos >= 29675 & pos <= 29903 ~ "3'UTR"))

# plotting
### define the color palette
my_pallete <- c(`5'UTR` = "#000000", ORF1ab = "#3cb44b", ITR1 = "#000000", S = "#4363d8", ITR2 = "#000000",
                ORF3a = "#f58231", ITR3 = "#000000", E = "#ffe119", ITR4 = "#000000", M = "#42d4f4",
                ITR5 = "#000000", ORF6a = "#e6194B", ITR6 = "#000000", ORF7ab = "#808000", ITR7 = "#000000",
                ORF8 ="#469990", ITR8 = "#000000", N = "#f032e6", ITR8 = "#000000", ORF10 = "#ffd8b1",
                `3'UTR` = "#000000")

### define the borders of the genes
mut_pos_lab <- c(266, 21563, 25393, 26245, 26523, 27202, 27394, 27894, 28274, 29558)

#### calculate the mean positions of the genes
x = c(-200, 10910, 23473, 28950)
x1 = c(26050, 30250)
x2 = c(27500)
x3 = c(29615.5)
x4 = c(21558.5, 25388, 29600)
x5 = c(26400)
x6 = c(27750)

#####################################################
# adjust y coordinates according to maximum coverage!
#####################################################

y = c(-810)
y1 = c(-800)
y2 = c(-1250)
y3 = c(-1250)
y4 = c(-2100)
y5 = c(-2400)

label = c("5'UTR", "ORF1ab", "S", "N")
label1 = c("ORF3a, E, M", "3'UTR")
label2 = c("ORF6, 7ab, 8")
label3 = c("ORF10")
label4 = c("ITR1", "ITR2", "ITR9")
label5 = c("ITR 3,4")
label6 = c("ITR 5,6,7,8")

median_cov <- c1$global_cov[1]
  
c1 %>% 
  distinct(pos, .keep_all = TRUE) %>% 
  ggplot(aes(pos, pos_median_cov, color = gene)) +
  geom_line() +
  geom_hline(aes(yintercept = global_cov), color = "red", lty = 2) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank()) +
  scale_color_manual(values = my_pallete) +
  scale_x_continuous(limits = c(-1000, 30500)) +
  scale_y_continuous(limits = c(-2400, 7000), breaks = c(seq(0, 7000, 1000))) +
  annotate("segment", x = 0, xend = 265, y = -300, yend = -300, color = "#000000", size = 4, alpha = 0.7) +
  annotate("segment", x = 266, xend = 21555, y = -300, yend = -300, color = "#3cb44b", size = 4, alpha = 0.7) +
  annotate("segment", x = 21556, xend = 21562, y = -1700, yend = -1700, color = "#000000", size = 4, alpha = 0.7) +
  annotate("segment", x = 21563, xend = 25384, y = -300, yend = -300, color = "#4363d8", size = 4, alpha = 0.7) +
  annotate("segment", x = 25385, xend = 25392, y = -1700, yend = -1700, color = "#000000", size = 4, alpha = 0.7) +
  annotate("segment", x = 25393, xend = 26220, y = -300, yend = -300, color = "#f58231", size = 4, alpha = 0.7) +
  annotate("segment", x = 26221, xend = 26244, y = -1700, yend = -1700, color = "#000000", size = 4, alpha = 0.7) +
  annotate("segment", x = 26245, xend = 26472, y = -300, yend = -300, color = "#ffe119", size = 4, alpha = 0.7) +
  annotate("segment", x = 26473, xend = 26522, y = -1700, yend = -1700, color = "#000000", size = 4, alpha = 0.7) +
  annotate("segment", x = 26523, xend = 27191, y = -300, yend = -300, color = "#42d4f4", size = 4, alpha = 0.7) +
  annotate("segment", x = 27192, xend = 27201, y = -1700, yend = -1700, color = "#000000", size = 4, alpha = 0.7) +
  annotate("segment", x = 27202, xend = 27387, y = -800, yend = -800, color = "#e6194B", size = 4, alpha = 0.7) +
  annotate("segment", x = 27388, xend = 27393, y = -1700, yend = -1700, color = "#000000", size = 4, alpha = 0.7) +
  annotate("segment", x = 27394, xend = 27887, y = -800, yend = -800, color = "#808000", size = 4, alpha = 0.7) +
  annotate("segment", x = 27888, xend = 27893, y = -1700, yend = -1700, color = "#000000", size = 4, alpha = 0.7) +
  annotate("segment", x = 27894, xend = 28259, y = -800, yend = -800, color = "#469990", size = 4, alpha = 0.7) +
  annotate("segment", x = 28260, xend = 28273, y = -1700, yend = -1700, color = "#000000", size = 4, alpha = 0.7) +
  annotate("segment", x = 28274, xend = 29533, y = -300, yend = -300, color = "#f032e6", size = 4, alpha = 0.7) +
  annotate("segment", x = 29534, xend = 29557, y = -1700, yend = -1700, color = "#000000", size = 4, alpha = 0.7) +
  annotate("segment", x = 29558, xend = 29674, y = -800, yend = -800, color = "#ffd8b1", size = 4, alpha = 0.7) +
  annotate("segment", x = 29675, xend = 29903, y = -300, yend = -300, color = "#000000", size = 4, alpha = 0.7) +
  annotate("text", x = x, y = y, label = label, color = "black", size = 3, fontface = 1) +
  annotate("text", x = x1, y = y1, label = label1, color = "black", size = 2, fontface = 2) +
  annotate("text", x = x2, y = y2, label = label2, color = "black", size = 2, fontface = 2) +
  annotate("text", x = x3, y = y3, label = label3, color = "black", size = 2, fontface = 2) +
  annotate("text", x = x4, y = y4, label = label4, color = "black", size = 2, fontface = 2) +
  annotate("text", x = x5, y = y5, label = label5, color = "black", size = 2, fontface = 2) +
  annotate("text", x = x6, y = y4, label = label6, color = "black", size = 2, fontface = 2) +
  #geom_text(aes(y = mean_cov + 200, x = -1000, label = mean_cov), color = "red", size = 4, hjust = 1) +
  labs(x = "SARS-CoV-2 reference genome (NC_045512.2)", y = "Median genome coverage per base")


ggsave("run_coverage.tiff", device = "tiff", path = path,
       width = 24, height = 18, units = "cm", dpi = 300)

### coverage for each sample
samples <- c1 %>% 
  distinct(sample) %>% 
  select(sample) %>% 
  unlist()
  

covplot_sample <- function(sample) {
  c1 %>% 
    filter(samp == sample) %>% 
    ggplot(aes(pos, cov, color = gene)) +
    geom_line() +
    geom_hline(aes(yintercept = sample_median_cov), color = "red", lty = 2) +
    theme_bw() +
    theme(legend.position = "none",
          panel.grid.minor = element_blank(),
          axis.line.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank()) +
    scale_color_manual(values = my_pallete) +
    scale_x_continuous(limits = c(-1000, 30500)) +
    scale_y_continuous(limits = c(-840, 6000), breaks = c(seq(0, 6000, 500))) +
    annotate("segment", x = 0, xend = 265, y = -80, yend = -80, color = "#000000", size = 4, alpha = 0.7) +
    annotate("segment", x = 266, xend = 21555, y = -80, yend = -80, color = "#3cb44b", size = 4, alpha = 0.7) +
    annotate("segment", x = 21556, xend = 21562, y = -600, yend = -600, color = "#000000", size = 4, alpha = 0.7) +
    annotate("segment", x = 21563, xend = 25384, y = -80, yend = -80, color = "#4363d8", size = 4, alpha = 0.7) +
    annotate("segment", x = 25385, xend = 25392, y = -600, yend = -600, color = "#000000", size = 4, alpha = 0.7) +
    annotate("segment", x = 25393, xend = 26220, y = -80, yend = -80, color = "#f58231", size = 4, alpha = 0.7) +
    annotate("segment", x = 26221, xend = 26244, y = -600, yend = -600, color = "#000000", size = 4, alpha = 0.7) +
    annotate("segment", x = 26245, xend = 26472, y = -80, yend = -80, color = "#ffe119", size = 4, alpha = 0.7) +
    annotate("segment", x = 26473, xend = 26522, y = -600, yend = -600, color = "#000000", size = 4, alpha = 0.7) +
    annotate("segment", x = 26523, xend = 27191, y = -80, yend = -80, color = "#42d4f4", size = 4, alpha = 0.7) +
    annotate("segment", x = 27192, xend = 27201, y = -600, yend = -600, color = "#000000", size = 4, alpha = 0.7) +
    annotate("segment", x = 27202, xend = 27387, y = -320, yend = -320, color = "#e6194B", size = 4, alpha = 0.7) +
    annotate("segment", x = 27388, xend = 27393, y = -600, yend = -600, color = "#000000", size = 4, alpha = 0.7) +
    annotate("segment", x = 27394, xend = 27887, y = -320, yend = -320, color = "#808000", size = 4, alpha = 0.7) +
    annotate("segment", x = 27888, xend = 27893, y = -600, yend = -600, color = "#000000", size = 4, alpha = 0.7) +
    annotate("segment", x = 27894, xend = 28259, y = -320, yend = -320, color = "#469990", size = 4, alpha = 0.7) +
    annotate("segment", x = 28260, xend = 28273, y = -600, yend = -600, color = "#000000", size = 4, alpha = 0.7) +
    annotate("segment", x = 28274, xend = 29533, y = -80, yend = -80, color = "#f032e6", size = 4, alpha = 0.7) +
    annotate("segment", x = 29534, xend = 29557, y = -600, yend = -600, color = "#000000", size = 4, alpha = 0.7) +
    annotate("segment", x = 29558, xend = 29674, y = -320, yend = -320, color = "#ffd8b1", size = 4, alpha = 0.7) +
    annotate("segment", x = 29675, xend = 29903, y = -80, yend = -80, color = "#000000", size = 4, alpha = 0.7) +
    annotate("text", x = x, y = y, label = label, color = "black", size = 3, fontface = 1) +
    annotate("text", x = x1, y = y1, label = label1, color = "black", size = 2, fontface = 2) +
    annotate("text", x = x2, y = y2, label = label2, color = "black", size = 2, fontface = 2) +
    annotate("text", x = x3, y = y3, label = label3, color = "black", size = 2, fontface = 2) +
    annotate("text", x = x4, y = y4, label = label4, color = "black", size = 2, fontface = 2) +
    annotate("text", x = x5, y = y5, label = label5, color = "black", size = 2, fontface = 2) +
    annotate("text", x = x6, y = y4, label = label6, color = "black", size = 2, fontface = 2) +
    #geom_text(aes(y = mean_cov + 200, x = -1000, label = mean_cov), color = "red", size = 4, hjust = 1) +
    labs(x = "SARS-CoV-2 reference genome (NC_045512.2)", y = "Median genome coverage per base")
  
  ggsave(filename = paste0("coverage_", sample, ".tiff"), device = "tiff",
         path = path,
         width = 24, height = 18, units = "cm", dpi = 300)
  
}

map(samples, covplot_sample)

## filter out positions with 0 coverage
zeroes <- c1 %>% 
  filter(pos_median_cov < 1) %>% 
  distinct(pos, .keep_all = TRUE)
