# ============================================================================
# 05_plot_data.R
# Plaque-size and figure plots: generate the main figures showing the
# relationship between extent of recoding and plaque size / fitness.
#
# Data dependencies (relative to repo root):
#   * results/*.csv           (compiled fitness data)
#   * data/design/construct_table.csv
#
# NOTE: as-written plotting script; not guaranteed turnkey. Depends on
# outputs from 01-04. Paths are relative to the code/ directory.
# ============================================================================

##script for visualizing phiX deoptimization data associated with R01
##relies on
####R_generate_deopt_seqs.R
####R_run_mfold.R

#libraries
{
  library(stringr)
  library(ggplot2)
  library(Biostrings)
  library(seqinr)
  library(ggsci)
  library(reshape2)
  library(viridis)
  library(cowplot)
}

##load up WT phiX data
phix.vals <- read.csv("data/phix_seq_vals.csv")

##plot distributions from random codon changes
per25 <- read.csv("data/G_rand_0.25per_vals.csv", header = T, row.names = 1)
per50 <- read.csv("data/G_rand_0.5per_vals.csv", header = T, row.names = 1)
per75 <- read.csv("data/G_rand_0.75per_vals.csv", header = T, row.names = 1)
per100 <- read.csv("data/G_rand_1per_vals.csv", header = T, row.names = 1)
max <- read.csv("data/G_rand_max_vals.csv", header = T, row.names = 1)
max$gene <- str_replace(max$gene,"_1_","_max_")



#look at new sequences
diffcount <- read.csv(file="results/diffcounts.csv")
diffcount.l <- melt(diffcount)
ggplot(data=diffcount.l, aes(x=variable, y=X, fill=value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme(text = element_text(size=8))
#ggsave(filename = "plots/AAdiffCount.pdf", width = 6, height= 20, units = "in")

tab.g <- read.csv("results/G_new_diffs.csv", header = T, stringsAsFactors = F)
names(tab.g) <- c(1:176, "id")
plot.g <- melt(tab.g, id.vars = c("id"))
index.g <- aggregate(. ~ id, plot.g[,c("id","value")], sum)
index.g <- index.g[order(index.g$value),]
plot.g$id <- factor(plot.g$id, levels=index.g$id)
ggplot(data=plot.g, aes(x=variable, y=id, fill=value)) +
  geom_tile() +
  scale_fill_gradient(low = "steelblue", high = "white") +
  theme(text = element_text(size=10), legend.position = "none", axis.text.x=element_text(angle=-45)) +
  xlab("") +
  ylab("")
#ggsave(filename = "plots/G_diffs.pdf", width = 16, height= 10, units = "in")

tab.f <- read.csv("results/F_new_diffs.csv", header = T, stringsAsFactors = F)
names(tab.f) <- c(1:428, "id")
plot.f <- melt(tab.f, id.vars = c("id"))
index.f <- aggregate(. ~ id, plot.f[,c("id","value")], sum)
index.f <- index.f[order(index.f$value),]
plot.f$id <- factor(plot.f$id, levels=index.f$id)
ggplot(data=plot.f, aes(x=variable, y=id, fill=value)) +
  geom_tile() +
  scale_fill_gradient(low = "#31a354", high = "white") +
  theme(text = element_text(size=10), legend.position = "none", axis.text.x=element_text(angle=-45)) +
  xlab("") +
  ylab("")
#ggsave(filename = "plots/F_diffs.pdf", width = 16, height= 10, units = "in")

tab.h <- read.csv("results/H_new_diffs.csv", header = T, stringsAsFactors = F)
names(tab.h) <- c(1:329, "id")
plot.h <- melt(tab.h, id.vars = c("id"))
index.h <- aggregate(. ~ id, plot.h[,c("id","value")], sum)
index.h <- index.h[order(index.h$value),]
plot.h$id <- factor(plot.h$id, levels=index.h$id)
ggplot(data=plot.h, aes(x=variable, y=id, fill=value)) +
  geom_tile() +
  scale_fill_gradient(low = "#756bb1", high = "white") +
  theme(text = element_text(size=10), legend.position = "none", axis.text.x=element_text(angle=-45)) +
  xlab("") +
  ylab("")
#ggsave(filename = "plots/H_diffs.pdf", width = 16, height= 10, units = "in")

cdat <- read.csv("Craig modeling/Fitness_effects_table.csv")
plq <- read.csv("data/PlqSize_my_results_summary.csv", row.names = 1)
plq <- plq[!plq$plq=='wtanc_1',]
##names need decoded again
names <- readxl::read_xlsx("data/Combined Fitness Assay Planning.xlsx", sheet=1)
plq[!str_detect(plq$plq, 'R'),]$plq <- paste('R', plq[!str_detect(plq$plq, 'R'),]$plq, sep='')
#how many problems are there going to be
table(plq$plq %in% names$`label of repick`)  ##good deal. only 2 to fix
plq[!plq$plq %in% names$`label of repick`,]  ##good deal. only 2 to fix
nrow(names)
nrow(plq)
nrow(cdat)
names.s <- names[,c('Construct', 'avg fitness', 'plq', 'label of repick')]
merge <- merge(plq, names.s, by.x='plq', by.y='label of repick')
merge <- merge[,c('Construct', 'plq', 'avg fitness', 'radius.getquant', 'radius.st.err', 'radius.length')]
cdat$X <- str_remove(cdat$X, 'redo')

#fix names
merge$Construct <- str_remove(merge$Construct, '%')
merge$Construct <- str_replace_all(merge$Construct, ' ', '_')
merge$Construct <- str_replace(merge$Construct, '0.5', '50')
merge$Construct <- str_replace(merge$Construct, '0.25', '25')
merge$Construct <- str_replace(merge$Construct, '0.75', '75')
merge$Construct <- str_replace(merge$Construct, '1.0', '100')
merge$Construct <- str_replace(merge$Construct, 'shuffle', 'shuf')
merge$Construct <- str_replace(merge$Construct, 'random', 'rand')
merge$Construct <- str_replace(merge$Construct, 'rep_', 'rep')


amount <- c('max','100','75','50','25')
type <- c('shuf', 'rand', 'worst','stab')

merge$treat <- merge$Construct
for(i in 1:nrow(merge)){
  con <- merge[i,]$Construct
  if(!str_detect(con, 'stab')){
    gene <- str_split(con, '_', simplify = T)[,1]
    lev <- as.character(na.omit(str_extract(con, amount)))
    mod <- as.character(na.omit(str_extract(con, type)))
    rep <- str_extract(con, 'rep\\d')
    merge[i,]$treat <- paste(gene, lev, mod, rep, sep='_')
  }
}
table(merge$treat %in% cdat$X)
merge[!merge$treat %in% cdat$X,]  ##a few G shuffles are missing. not in Craigs output csv. figure out why

merge2 <- merge(merge, cdat, by.x='treat', by.y='X')

p1 <- ggplot(merge2, aes(x=Flask, y=EN.Reg, size=radius.getquant, color=radius.getquant)) +
  geom_point() +
  scale_color_viridis() +
  theme_bw() +
  labs(x='Fitness Effect (Flask)', y='Fitness Effect (PacBio)')
p1
ggsave(filename = "plots/plq_size.png", width = 6.5, height= 5, units = "in")
ggsave(filename = "plots/plq_size.pdf", width = 6.5, height= 5, units = "in")


p2 <- ggplot(merge2, aes(x=Flask, y=radius.getquant, size=radius.getquant, color=radius.getquant)) +
  geom_point() +
  scale_color_viridis() +
  theme_bw() +
  labs(x='Fitness Effect (Flask)', y='Plaque size')

p3 <- ggplot(merge2, aes(x=EN.Reg, y=radius.getquant, size=radius.getquant, color=radius.getquant)) +
  geom_point() +
  scale_color_viridis() +
  theme_bw() +
  labs(x='Fitness Effect (PacBio)', y='Plaque size')

plot_grid(p1+theme(legend.position = 'none'),p2+theme(legend.position = 'none'),p3, nrow = 1, rel_widths = c(0.75,0.75,1))
ggsave(filename = "plots/plq_size_all.png", width = 14, height= 4, units = "in")
ggsave(filename = "plots/plq_size_all.pdf", width = 14, height= 4, units = "in")





