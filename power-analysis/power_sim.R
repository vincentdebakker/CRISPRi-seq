#### Naïve power simulations sgRNA enrichment calling for CRISPRi-seq ####
# Author: Vincent de Bakker
# 2020 NOV 11
# Veening Lab, DMF, FBM, University of Lausanne, Switzerland
# vincent.debakker@unil.ch
####

# PACKAGES ----
library(DESeq2)
library(parallel)
library(reshape2)
library(ggplot2)
source("/makeData.R")

# SETTINGS ----
# effect sizes (Log2FC)
E <- seq(1, 4, by = 0.5)
# nr of replicates per treatment group
N <- c(3:5, 8, 10)
# nr of sgRNAs per sample
G <- 1499
# nr of sgRNAs DE 
DEG <- c(1, 5, 10, seq(50, G / 2, by = 50))
# nr of reads per sgRNA
R <- c(10, 50, 100, 200, 300, 400, 500, 750, 1000)
# nr of MC simulation rounds per combination
M <- 10

# dispersion-mean counts relationship for simulation
## based very roughly on own data
dispMeanRel_emp <- function(mean_count){4 / mean_count + 0.1}
# Hypothesis test settings
alpha <- 0.05
LFCthr <- 1 # so FC of 2^1 = 2: doubling/halving
# cluster
ncores <- detectCores() - 2
# seed for reproducibility of "random" draws
seed <- 1992

# SIMULATIONS ----
message(paste(length(E) * length(N) * length(G) * length(DEG) * length(R) * M, 
              "combinations to run"))
cl <- makeCluster(ncores)
clusterExport(cl, varlist = c("N", "G", "R", "E", "M", "alpha", "LFCthr", "dispMeanRel_emp"))
clusterEvalQ(cl, {
  library(DESeq2)
  source("C:/Users/vince/Documents/PhD/Projects/CRISPRi-seq_NatProtoc/simulations_reads_per_sgRNA/makeData.R")
})
start <- Sys.time()
set.seed(seed)
sim_ls <- parLapply(cl, DEG, function(d){
  do.call(rbind, lapply(N, function(n){
    do.call(rbind, lapply(G, function(g){
      do.call(rbind, lapply(R, function(r){
        do.call(rbind, lapply(E, function(e){
          do.call(rbind, lapply(1:M, function(m){
            dds <- makeExampleDESeqDataSet_vdb(n = g, m = 2 * n, 
                                               dispMeanRel = dispMeanRel_emp, 
                                               interceptMean = log2(r), 
                                               interceptSD = 0, 
                                               betaMU = c(rep(-e, d), 
                                                          rep(0, g - d)), 
                                               betaSD = 0)
            dds <- DESeq(dds)
            res <- results(dds, contrast = c("condition", "B", "A"), 
                           alpha = alpha, lfcThreshold = LFCthr)
            # change NA p-values to 1 for summing
            res$padj <- ifelse(is.na(res$padj), 1, res$padj)
            # fraction of correctly rejected H0
            pow <- sum(res$padj[1:d] < alpha) / d
            # fraction of incorrectly rejected H0
            typeIerr <- sum(res$padj[(d + 1):g] < alpha) / (g - d)
            # return
            cbind("power" = pow, "typeIerr" = typeIerr, 
                  "trueDEG" = d, 
                  "Nreplicates" = n, 
                  "NsgRNAs" = g, 
                  "Nreads" = r, 
                  "effectSize" = e, 
                  "NsimMC" = m)
          }))
        }))
      }))
    }))
  }))
})
end <- Sys.time()
end - start
save(sim_ls, file = paste0(Sys.Date(), "_sim-ls_LFCthr-", LFCthr, "_alpha-", alpha, ".RData"))
stopCluster(cl)

# TIDY UP ----
sim_df <- as.data.frame(do.call(rbind, sim_ls))
# mean / SEM df
SEM <- function(x){
  sd(x) / sqrt(length(x))
}
sim_sum_df <- melt(tapply(sim_df$power, 
                          list(sim_df$trueDEG, 
                               sim_df$Nreplicates, 
                               sim_df$NsgRNAs, 
                               sim_df$Nreads, 
                               sim_df$effectSize),
                          mean))
colnames(sim_sum_df) <- c("trueDE", "Nreplicates", "NsgRNAs", "Nreads", "trueLFC", "meanPower")
sim_sum_df$semPower <- melt(tapply(sim_df$power, 
                              list(sim_df$trueDEG, 
                                   sim_df$Nreplicates, 
                                   sim_df$NsgRNAs, 
                                   sim_df$Nreads, 
                                   sim_df$effectSize),
                              SEM))$value
# add type I error rate
sim_sum_df$meanTypeIerr <- melt(tapply(sim_df$typeIerr, 
                                       list(sim_df$trueDEG, 
                                            sim_df$Nreplicates, 
                                            sim_df$NsgRNAs, 
                                            sim_df$Nreads, 
                                            sim_df$effectSize),
                                       mean))$value
sim_sum_df$semTypeIerr <- melt(tapply(sim_df$typeIerr, 
                                      list(sim_df$trueDEG, 
                                           sim_df$Nreplicates, 
                                           sim_df$NsgRNAs, 
                                           sim_df$Nreads, 
                                           sim_df$effectSize),
                                      SEM))$value

# PLOT ----
# mean / SEM df
ggplot(sim_sum_df, 
       aes(Nreads, meanPower, 
           col = as.factor(Nreplicates))) + 
  geom_hline(yintercept = 0.9, col = "grey", lty = 2, size = 1) + 
  geom_point(size = 1) + 
  geom_line(size = 1) + 
  geom_errorbar(aes(ymin = meanPower - semPower, 
                    ymax = meanPower + semPower), 
                size = 1) + 
  facet_grid(trueLFC ~ trueDE, labeller = "label_both") + 
  scale_color_brewer(palette = "Dark2") + 
  #scale_x_continuous(breaks = R) + 
  labs(title = paste0("tested with alpha: ", alpha, " and LFC threshold: ", LFCthr), 
       subtitle = paste0("Monte Carlo simulation rounds per combination with ", G, " sgRNAs: ", M), 
       y = "Average power", 
       x = "Average nr. of reads (counts) per non-induced / unaffected sgRNA", 
       col = "Biological \nreplicates \nper condition", 
       caption = "Averages with Standard Error of the Mean (SEM) as error bars. \nGrey dotted line at power of 0.9.") +
  theme_bw(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 40, 
                                   vjust = 1, 
                                   hjust = 1), 
        panel.grid = element_blank())


# Type I error rate
ggplot(sim_sum_df, 
       aes(Nreads, meanTypeIerr, 
           col = as.factor(Nreplicates))) + 
  geom_hline(yintercept = 0.05, col = "grey", lty = 2, size = 1) + 
  geom_point(size = 1) + 
  geom_line(size = 1) + 
  geom_errorbar(aes(ymin = meanTypeIerr - semTypeIerr, 
                    ymax = meanTypeIerr + semTypeIerr), 
                size = 1) + 
  facet_grid(trueLFC ~ trueDE, labeller = "label_both") + 
  scale_color_brewer(palette = "Dark2") + 
  #scale_x_continuous(breaks = R) + 
  labs(title = paste0("tested with alpha: ", alpha, " and LFC threshold: ", LFCthr), 
       subtitle = paste0("Monte Carlo simulation rounds per combination with ", G, " sgRNAs: ", M), 
       y = "Average type I error rate", 
       x = "Average nr. of reads (counts) per non-induced / unaffected sgRNA", 
       col = "Biological \nreplicates \nper condition", 
       caption = "Averages with Standard Error of the Mean (SEM) as error bars. \nGrey dotted line at power of 0.05.") +
  theme_bw(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 40, 
                                   vjust = 1, 
                                   hjust = 1), 
        panel.grid = element_blank())

# save
ggpower <- ggplot(subset(sim_sum_df, trueDE > 5 & trueDE < 550), 
                  aes(Nreads, meanPower, 
                      col = as.factor(Nreplicates))) + 
  geom_hline(yintercept = 0.9, col = "grey", lty = 2, size = 1) + 
  geom_point(size = 1) + 
  geom_line(size = 1) + 
  geom_errorbar(aes(ymin = meanPower - semPower, 
                    ymax = meanPower + semPower), 
                size = 1) + 
  facet_grid(trueLFC ~ trueDE, labeller = "label_both") + 
  scale_color_brewer(palette = "Dark2") + 
  #scale_x_continuous(breaks = R) + 
  labs(y = "Average power", 
       x = "Average nr. of reads (counts) per non-induced / unaffected sgRNA", 
       col = "Biological \nreplicates \nper condition") +
  theme_bw(base_size = 14) + 
  theme(axis.text.x = element_text(angle = 40, 
                                   vjust = 1, 
                                   hjust = 1), 
        panel.grid = element_blank())
ggpower
ggsave("Fig_Sx_naive-power-analysis.jpg", 
       width = 16, height = 10, units = "cm", scale = 3, 
       dpi = 300)
