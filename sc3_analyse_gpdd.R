
library(tidyverse)
library(foreach)
library(patchwork)

gpdd_data <- read.csv("data_gpdd/data.csv")
gpdd_main <- read.csv("data_gpdd/main.csv")
gpdd_taxa <- read.csv("data_gpdd/taxon.csv")

source("functions_pe.R")

# Filter GPDD data --------------------------------------------------------

# Time series length larger than 10
n <- table(gpdd_data$MainID)
idx_n <- which(n >= 10)

gpdd_data1 <- filter(gpdd_data, MainID %in% names(idx_n))

# No time skips
step_length <- function (x) {
  x[2:length(x)] - x[1:(length(x)-1)]
}

idx_step <-group_by(gpdd_data1, MainID) %>%
  summarise(steps = step_length(SampleYear)) %>%
  summarise(step_skip = any(steps > 1)) %>%
  ungroup() %>%
  filter(step_skip == FALSE) %>%
  select(MainID)

gpdd_data2 <- filter(gpdd_data1, MainID %in% idx_step$MainID)

# No count replicated more than 15%
idx_rep <- group_by(gpdd_data2, MainID) %>%
  summarise(n = table(Population)/sum(table(Population))) %>%
  summarise(rep = any(n > 0.15)) %>%
  filter(rep == FALSE) %>%
  select(MainID)

gpdd_data3 <- filter(gpdd_data2, MainID %in% idx_rep$MainID)

# No data transformation
gpdd_main2 <- filter(gpdd_main, MainID %in% unique(gpdd_data3$MainID))
idx_tr <- which(gpdd_main2$SourceTransform == "None" | 
                  gpdd_main2$SourceTransformReference == "Log (e)")
gpdd_main3 <- gpdd_main2[idx_tr,]

gpdd_data4 <- filter(gpdd_data3, MainID %in% gpdd_main3$MainID)
n2 <- table(gpdd_data4$MainID)


# Calculate WPE for each time series --------------------------------------

ts <- unique(gpdd_data4$MainID)

gpdd_wpe <- foreach(i = 1:length(ts), .combine = "c") %do% {
  
  x <- filter(gpdd_data4, MainID == ts[i])$Population
  x2 <- filter(gpdd_main3, MainID == ts[i])
  
  if (x2$SourceTransform == "Log") x <- exp(x)
  
  calculate_pe(x)[1]
  
}

gpdd_pe <- foreach(i = 1:length(ts), .combine = "c") %do% {
  
  x <- filter(gpdd_data4, MainID == ts[i])$Population
  x2 <- filter(gpdd_main3, MainID == ts[i])
  
  if (x2$SourceTransform == "Log") x <- exp(x)
  
  calculate_pe(x)[2]
  
}


# Apply WPE randomization test to time series  ----------------------------

pb <- txtProgressBar(0, length(ts), style = 3)

p_pe_wh <- foreach(i = 1:length(ts), .combine = "c") %do% {
  
  x <- filter(gpdd_data4, MainID == ts[i])$Population
  x2 <- filter(gpdd_main3, MainID == ts[i])
  
  if (x2$SourceTransform == "Log") x <- exp(x)
  
  y <- test_pe(x, n_random = 1000)
  
  setTxtProgressBar(pb, i)
  
  return(y[3])
  
}

p_wpe_wh <- foreach(i = 1:length(ts), .combine = "c") %do% {
  
  x <- filter(gpdd_data4, MainID == ts[i])$Population
  x2 <- filter(gpdd_main3, MainID == ts[i])
  
  if (x2$SourceTransform == "Log") x <- exp(x)
  
  y <- test_wpe(x, n_random = 1000)
  
  setTxtProgressBar(pb, i)
  
  return(y[3])
  
}

res_gpdd<- list(p_wpe_wh = p_wpe_wh,
                gpdd_wpe = gpdd_wpe,
                gpdd_data = gpdd_data4)

saveRDS(res_gpdd, "results/res_gpdd.rds")


# Summarize predictability across taxa ------------------------------------

sig_wpe_wh <- as.factor(ifelse(p_wpe_wh > 0.05 & n2 < 30, 1, 
                               ifelse(p_wpe_wh > 0.05 & n2 > 30, 2, 0)))

gpdd_main4 <- filter(gpdd_main, MainID %in% gpdd_data4$MainID) %>%
  left_join(gpdd_taxa, by = "TaxonID") %>%
  select(MainID, TaxonID, TaxonomicClass, SourceTransform, Reliability) %>%
  add_column(p_wpe_wh = p_wpe_wh) %>%
  add_column(sig_wpe_wh = sig_wpe_wh) %>%
  add_column(wpe = gpdd_wpe) %>%
  filter(TaxonomicClass %in% c("Aves", "Insecta", "Mammalia", "Osteichthyes"))

df1 <- gpdd_main4 %>%
  mutate(TaxonomicClass = factor(TaxonomicClass, 
                                 levels = c("Mammalia", "Aves", "Insecta", 
                                            "Osteichthyes"))) %>%
  group_by(TaxonomicClass, sig_wpe_wh) %>% 
  summarise(n = n()) %>%
  mutate(per = n/sum(n)) %>%
  ungroup()

theme_set(theme_bw())

g_gpdd1 <- ggplot(data = df1, 
                  mapping = 
                    aes(x = factor(TaxonomicClass, levels = levels(TaxonomicClass)[4:1]), 
                        y = per, 
                        group = factor(sig_wpe_wh, levels = levels(sig_wpe_wh)[3:1]),
                        fill = factor(sig_wpe_wh, levels = levels(sig_wpe_wh)[3:1]))) + 
  geom_col(width = 0.75, position = 'dodge', size = 0.5, color = "black", 
           alpha = 0.75) +
  geom_label(aes(y = per + 0.05,  label = round(per, 2)), size = 3, 
             position = position_dodge(0.75), colour = "white", 
             show.legend = F) +
  labs(y = "Frequency") +
  theme(legend.position = "bottom",
        legend.key.size = unit(3, "mm"),
        title = element_text(size = 10),
        legend.text = element_text(size = 10),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank()) +
  scale_fill_manual(name = NULL,
                    labels = c("1" = "p > 0.05, t < 30",
                               "2" = "p > 0.05, t > 30",
                               "0" ="p < 0.05"),
                    values = c("0" = "blue4", 
                               "1" = "darkred",
                               "2" = "darkorchid3")) +
  scale_y_continuous(limits = c(0, 0.6), breaks = seq(0, 0.5, 0.1)) +
  coord_flip()

g_gpdd2 <- ggplot() +
  #geom_jitter(aes(y = gpdd_wpe, x = sig_wpe_wh), alpha = 0.5, 
  #width = 0.15) +
  geom_violin(aes(y = gpdd_wpe, x = sig_wpe_wh, fill = sig_wpe_wh),
              alpha = 0.75, linewidth = 0.5) +
  labs(y = "WPE") +
  scale_fill_manual(values = c("0" = "blue4", 
                               "1" = "darkred",
                               "2" = "darkorchid3")) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  scale_x_discrete(labels = c("1" = "p > 0.05,\n t < 30",
                              "2" = "p > 0.05,\n t > 30",
                              "0" ="p < 0.05")) +
  theme(legend.position = "none",
        title = element_text(size = 10),
        legend.text = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        panel.border = element_blank(),
        panel.grid.minor = element_blank())

g_gpdd2 + g_gpdd1 +
  plot_annotation(tag_levels = "a") +
  plot_layout(widths = c(0.7, 1))
ggsave("figures/plot_gpdd.pdf", width = 180, height = 120, units = "mm",
       dpi = 600)

# Statistics for manuscript results
r1 <- filter(gpdd_main4, TaxonomicClass == "Mammalia")$sig_wpe_wh %>%
  table()
r1/sum(r1)  

r2 <- filter(gpdd_main4, TaxonomicClass == "Aves")$sig_wpe_wh %>%
  table()
r2/sum(r2)

r3 <- filter(gpdd_main4, TaxonomicClass == "Insecta")$sig_wpe_wh %>%
  table()
r3/sum(r3)

r4 <- filter(gpdd_main4, TaxonomicClass == "Osteichthyes")$sig_wpe_wh %>%
  table()
r4/sum(r4)


# Plot for SI
gpdd_main4_v2 <- filter(gpdd_main, MainID %in% gpdd_data4$MainID) %>%
  left_join(gpdd_taxa, by = "TaxonID") %>%
  select(MainID, TaxonID, TaxonomicClass, SourceTransform, Reliability) %>%
  add_column(p_wpe_wh = p_wpe_wh) %>%
  add_column(sig_wpe_wh = sig_wpe_wh) %>%
  add_column(wpe = gpdd_wpe) %>%
  filter(!TaxonomicClass %in% c("Aves", "Insecta", "Mammalia", "Osteichthyes"))

df1_v2 <- gpdd_main4_v2 %>%
  mutate(TaxonomicClass = factor(TaxonomicClass)) %>%
  group_by(TaxonomicClass, sig_wpe_wh) %>% 
  summarise(n = n()) %>%
  mutate(per = n/sum(n)) %>%
  ungroup()

ggplot(data = df1_v2, 
       mapping = 
       aes(x = TaxonomicClass, 
           y = per, 
           group = factor(sig_wpe_wh, levels = levels(sig_wpe_wh)[3:1]),
           fill = factor(sig_wpe_wh, levels = levels(sig_wpe_wh)[3:1]))) + 
  geom_col(width = 0.75, position = 'dodge', size = 1.1, color = "black", 
           alpha = 0.75) +
  geom_label(aes(y = per + 0.050,  label = paste("n", "=", n)), 
             position = position_dodge(0.75), colour = "white", 
             show.legend = F, size = 3) +
  labs(y = "Frequency") +
  theme(legend.position = "bottom",
        title = element_text(size = 18),
        legend.text = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank()) +
  scale_fill_manual(name = NULL,
                    labels = c("1" = "p > 0.05, t < 30",
                               "2" = "p > 0.05, t > 30",
                               "0" ="p < 0.05"),
                    values = c("0" = "blue4", 
                               "1" = "darkred",
                               "2" = "darkorchid3")) +
  coord_flip()

ggsave("plot_gpdd_SI.pdf", width = 12, height = 8, units = "in")


# Reliability plots -------------------------------------------------------

g_r1 <- ggplot(data = gpdd_main4) +
  geom_boxplot(aes(x = Reliability, group = Reliability, y = wpe))

dfr <- gpdd_main4 %>%
  group_by(Reliability, sig_wpe_wh) %>% 
  summarise(n = n()) %>%
  mutate(per = n/sum(n)) %>%
  ungroup()

g_r2 <- ggplot(data = dfr, 
       mapping = 
         aes(x = factor(Reliability), 
             y = per, 
             group = factor(sig_wpe_wh, levels = levels(sig_wpe_wh)[3:1]),
             fill = factor(sig_wpe_wh, levels = levels(sig_wpe_wh)[3:1]))) + 
  geom_col(width = 0.75, position = 'dodge', size = 1.1, color = "black", 
           alpha = 0.75) +
  geom_label(aes(y = per + 0.025,  label = round(per, 2)), 
             position = position_dodge(0.75), colour = "white", 
             show.legend = F) +
  labs(y = "Frequency", x = "Reliability") +
  theme(legend.position = "bottom",
        title = element_text(size = 18),
        legend.text = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 14),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()) +
  scale_fill_manual(name = NULL,
                    labels = c("1" = "p > 0.05, t < 30",
                               "2" = "p > 0.05, t > 30",
                               "0" ="p < 0.05"),
                    values = c("0" = "blue4", 
                               "1" = "darkred",
                               "2" = "darkorchid3")) +
  scale_y_continuous(limits = c(0, 0.7), breaks = seq(0, 0.7, 0.1))

g_r1 + g_r2 +
  plot_annotation(tag_levels = "a")
ggsave("plot_reliability.pdf", width = 17, height = 8, units = "in")


# Calculate trends --------------------------------------------------------

ts2 <- unique(gpdd_main4$MainID)

trends <- foreach(i = 1:length(ts2), .combine = "c") %do% {
  
  x <- filter(gpdd_data4, MainID == ts2[i])$Population
  if (i == 743) x <- filter(gpdd_data4, MainID == ts2[i])$PopulationUntransformed
  x2 <- filter(gpdd_main4, MainID == ts2[i])
  
  if (x2$SourceTransform == "Log") x <- exp(x)
  
  r <- log(x[2:length(x)]) - log(x[1:(length(x)-1)])
  
  r <- r[!is.infinite(r)]
  
  mean(r, na.rm = T)
  
}

gpdd_main5 <- gpdd_main4 %>%
  add_column(trends = trends)

g1 <- ggplot(data = filter(gpdd_main5, TaxonomicClass ==  "Mammalia")) +
  geom_density(mapping = aes(x = trends, fill = sig_wpe_wh), alpha = 0.6) +
  labs(y = "Density", x = "Average Growth (log)", title = "Mammalia") +
  scale_fill_manual(name = NULL,
                    labels = c("1" = "p > 0.05, t < 30",
                               "2" = "p > 0.05, t > 30",
                               "0" ="p < 0.05"),
                    values = c("0" = "blue4", 
                               "1" = "darkred",
                               "2" = "darkorchid3")) +
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  scale_y_continuous(breaks = c(0,6,12), limits = c(0,15)) +
  scale_x_continuous(breaks = c(-0.2,0,0.2,0.4,0.6))

g2 <- ggplot(data = filter(gpdd_main5, TaxonomicClass ==  "Aves")) +
  geom_density(mapping = aes(x = trends, fill = sig_wpe_wh), alpha = 0.6) +
  labs(y = "Density", x = "Average Growth (log)", title = "Aves") +
  scale_fill_manual(name = NULL,
                    labels = c("1" = "p > 0.05, t < 30",
                               "2" = "p > 0.05, t > 30",
                               "0" ="p < 0.05"),
                    values = c("0" = "blue4", 
                               "1" = "darkred",
                               "2" = "darkorchid3")) +
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  scale_y_continuous(breaks = c(0,6,12), limits = c(0,15))

g3 <- ggplot(data = filter(gpdd_main5, TaxonomicClass ==  "Insecta")) +
  geom_density(mapping = aes(x = trends, fill = sig_wpe_wh), alpha = 0.6) +
  labs(y = "Density", x = "Average Growth (log)", title = "Insecta") +
  scale_fill_manual(name = NULL,
                    labels = c("1" = "p > 0.05, t < 30",
                               "2" = "p > 0.05, t > 30",
                               "0" ="p < 0.05"),
                    values = c("0" = "blue4", 
                               "1" = "darkred",
                               "2" = "darkorchid3")) +
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  scale_x_continuous(breaks = c(-0.4, -0.2,0,0.2,0.4,0.6))

g4 <- ggplot(data = filter(gpdd_main5, TaxonomicClass ==  "Osteichthyes")) +
  geom_density(mapping = aes(x = trends, fill = sig_wpe_wh), alpha = 0.6) +
  labs(y = "Density", x = "Average Growth (log)", title = "Osteichthyes") +
  scale_fill_manual(name = NULL,
                    labels = c("1" = "p > 0.05, t < 30",
                               "2" = "p > 0.05, t > 30",
                               "0" ="p < 0.05"),
                    values = c("0" = "blue4", 
                               "1" = "darkred",
                               "2" = "darkorchid3")) +
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none") +
  scale_y_continuous(breaks = c(0,3,6))

ggplot(data = filter(gpdd_main5, TaxonomicClass ==  "Mammalia")) +
  geom_density(mapping = aes(x = trends, fill = sig_wpe_wh), alpha = 0.75) +
  labs(y = "Density", x = "Avereage Growth (log)") +
  scale_fill_manual(name = NULL,
                    labels = c("1" = "p > 0.05, t < 30",
                               "2" = "p > 0.05, t > 30",
                               "0" ="p < 0.05"),
                    values = c("0" = "blue4", 
                               "1" = "darkred",
                               "2" = "darkorchid3")) +
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")
ggsave("lgd.jpeg", width = 10, height = 8, units = "in")


(g1 + g2) / (g3 + g4)
ggsave("plot_gpdd_trends.pdf", width = 10, height = 8, units = "in")

# r vs WPE
ggplot() +
  geom_point(aes(y = gpdd_main5$wpe, x = gpdd_main5$trends, 
                 col = gpdd_main5$sig_wpe_wh), 
             alpha = 0.8, size = 3) +
  scale_color_manual(name = NULL,
                     labels = c("1" = "p > 0.05, t < 30",
                                "2" = "p > 0.05, t > 30",
                                "0" ="p < 0.05"),
                     values = c("0" = "blue4", 
                                "1" = "darkred",
                                "2" = "darkorchid3")) +
  labs(y = "WPE", x = "r") +
  theme(panel.border = element_blank(),
       panel.grid.minor = element_blank(),
       legend.position = "bottom")
ggsave("plot_WPEvsR.jpeg", width = 10, height = 8, units = "in")

