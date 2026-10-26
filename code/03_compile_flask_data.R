# ============================================================================
# 03_compile_flask_data.R
# Flask assay data assembly: read the four Excel assay workbooks, normalise
# construct identifiers, and compile the per-assay fitness table.
#
# Data dependencies (relative to repo root):
#   * data/*.xlsx           (4 assay workbooks)
#   * data/design/construct_table.csv
#
# NOTE: as-written pipeline script; not guaranteed turnkey. The
# Construct Failures.xlsx read is commented out but the `fails` variable
# is still referenced downstream, so a full run will error at that point.
# Paths are relative to the code/ directory.
# ============================================================================

#libraries
{
library(seqinr)
library(stringr)
library(reshape2)
library(ggplot2)
library(dplyr)
library(readxl)
library(ggsci)
library(cowplot)
library(ggthemes)
library(tidyr)
library(ggrepel)
}

lu2 <- read_xlsx("../data/assay/Constructs_Transcript Stability_Worst.xlsx", sheet = 5, skip = 3) 
lu2 <- lu2[which(!is.na(lu2[,1])),]
LuComments <- rbind(lu2[,c(1,38)])
LuComments <- LuComments[!is.na(LuComments[,2]),]
lu2 <- lu2[,1:36]  ##drop comments and average fit
lu.l <- melt(lu2, id.vars = c('...1', 'sample #'))
lu.l <- lu.l[which(!str_detect(lu.l$value, '[a-zA-Z]')),]  ##drop cells with comments and stuff
names(lu.l) <- c('construct', 'sample', 'date', 'fitness')
lu.l$date <- openxlsx::convertToDate(as.character(lu.l$date))
lu.l$fitness <- round(as.numeric(lu.l$fitness), digits=2)

emma <- read_xlsx("../data/assay/Constructs_Shuffles_Randoms.xlsx", sheet = 4, skip = 7)
emma.s <- na.omit(emma[,c("Construct", "Date", "flask #", "fitness")])

emma2 <- read_xlsx("../data/assay/Constructs_Shuffles_Randoms.xlsx", sheet = 5) 
emma2 <- emma2[which(!is.na(emma2$Construct)),]
EmmaComments <- emma2[,c('Construct', 'Photo?', 'times assembled', 'notes', 'not enough/too many plqs on plate to count')]
emma2 <- emma2[,1:40]  ##drop comments and average fit
emma.l <- melt(emma2, id.vars = c('Sample #', 'Construct'))
emma.l <- emma.l[which(!str_detect(emma.l$value, '[a-zA-Z]')),]  ##drop cells with comments and stuff
names(emma.l) <- c( 'sample', 'construct', 'date', 'fitness')
emma.l$date <- openxlsx::convertToDate(as.character(emma.l$date))
emma.l$fitness <- round(as.numeric(emma.l$fitness), digits=2)
emma.l[emma.l$construct=='wt Sanger2mut (10^3)',]$construct <- 'wt Sanger2mut'


lu.cnt <- nrow(lu.l[str_detect(lu.l$construct, "wt"),])
emma.cnt <- nrow(emma.l[str_detect(emma.l$construct, "wt"),])

wt <- rbind(emma.l[str_detect(emma.l$construct, "wt"),], lu.l[str_detect(lu.l$construct, "wt"),])
wt$sci <- c(rep("emma",emma.cnt), rep("lu", lu.cnt))
wt$construct <- "wt"
ggplot(wt, aes(x=sci, y=fitness, color=sci, fill=sci)) +
  geom_violin(alpha=0.2, show.legend = FALSE) +
  geom_jitter(size=3, width = 0.1, show.legend = FALSE) +
  scale_color_aaas() +
  scale_fill_aaas() +
  theme_classic()
#ggsave(filename = "plots/wt_fitness.pdf", width = 5.5, height= 5.5, units = "in")
#ggsave(filename = "plots/wt_fitness.png", width = 5.5, height= 5.5, units = "in")
  

emma.s <- emma.l
##more formatting
emma.s$sci <- rep("emma",nrow(emma.s))
emma.s$gene <- str_split(emma.s$construct, " ", simplify = T)[,1]
emma.s$type <- str_split(emma.s$construct, " ", simplify = T)[,2]
emma.s$treat <- str_split(emma.s$construct, " ", simplify = T)[,3]
emma.s$rep <- str_split(emma.s$construct, " ", simplify = T)[,5]
emma.s[emma.s$gene=="wt",]$treat <- "0%"
emma.s$treat_numeric <- str_replace(emma.s$treat, "%", "")
emma.s$treat_numeric <- str_replace(emma.s$treat_numeric, "max", "125")
emma.s$treat_numeric <- as.numeric(emma.s$treat_numeric)/100
emma.s[emma.s$rep=="",]$rep <- "1"
emma.s$type <- str_replace(emma.s$type, "Sanger2mut", "wt")

lu.s <- lu.l
lu.s$sci <- rep("lu", nrow(lu.s))
lu.s$gene <- str_replace(str_extract(lu.s$construct, "[FGH]_"), "_", "")
lu.s[str_detect(lu.s$construct, "wt"),]$gene <- "wt"
lu.s$type <- "NA"
lu.s[str_detect(lu.s$construct, "shuf"),]$type <- "shuffle"
lu.s[str_detect(lu.s$construct, "worst"),]$type <- "worst"
lu.s[str_detect(lu.s$construct, "Hi"),]$type <- "shuffle_Hi"
lu.s[str_detect(lu.s$construct, "Lo"),]$type <- "shuffle_Lo"
lu.s[str_detect(lu.s$construct, "wt"),]$type <- "wt"
lu.s$treat <- str_extract(lu.s$construct, "[\\d\\.]{2,}")
lu.s[lu.s$gene=="wt",]$treat <- "0"
lu.s[str_detect(lu.s$type, "shuffle"),]$treat <- "0.5"
lu.s$treat <- paste(as.character(as.numeric(lu.s$treat)*100), "%", sep="")
lu.s$rep <- str_replace(str_extract(lu.s$construct, "rep\\d"), "rep", "")
lu.s[str_detect(lu.s$construct, "stab"),]$rep <- str_split(lu.s[str_detect(lu.s$construct, "stab"),]$construct,"_",simplify = T)[,3]
lu.s$rep <- str_replace(lu.s$rep,"Lo","")
lu.s$rep <- str_replace(lu.s$rep,"Hi","")
lu.s[lu.s$gene=="wt",]$rep <- "1"
lu.s$treat_numeric <- as.numeric(str_replace(lu.s$treat, "%", ""))/100


## Construct Failures.xlsx is not distributed with this repo. Provide it locally to rebuild the 11 construction-failure rows.
## fails <- read_excel("../data/assay/Construct Failures.xlsx", sheet = 1, col_names = F)
fails <- fails[which(!is.na(fails$...2)),]
emma.fails <- which(fails$...2 == 'Emma')
lu.fails <- which(fails$...2 == 'Lu')
names(fails) <- c('date', 'sci', 'construct')
fails <- fails[fails$sci %in% c('Emma','Lu'), 1:3]
fails$construct <- str_replace(fails$construct, 'transcript stability_Lo1', 'shuffle-Lo_50%_rep1')
fails$construct <- str_replace(fails$construct, 'H_0.25worst_rep3_redo3_O/N2', 'H_worst_25%_rep3')
fails$construct <- str_replace(fails$construct, 'H_0.25worst_rep1', 'H_worst_25%_rep1')
fails$sample <- 1
fails$fitness <- 0
fails$gene <- str_split(fails$construct, '_', simplify = T)[,1]
fails$type <-str_split(fails$construct, '_', simplify = T)[,2]
fails$treat <- str_split(fails$construct, '_', simplify = T)[,3]
fails$rep <- str_split(fails$construct, '_', simplify = T)[,4]
fails$treat_numeric <- str_split(fails$construct, '_', simplify = T)[,3]
fails$treat <- str_replace(fails$treat, 'max', '125%') 
fails$treat_numeric <- as.numeric(str_replace(fails$treat, "%", ""))/100
fails$treat <- str_replace(fails$treat, '125%', 'max') 
fails$date <- openxlsx::convertToDate(as.character(fails$date))
fails$sci <- tolower(fails$sci)
fails$type <- str_replace(fails$type, '-', '_')

emma.s <- emma.s[,names(lu.s)]
fails <- fails[,names(lu.s)]
head(emma.s,1)
head(lu.s,1)
head(fails,1)
fit.all <- rbind(emma.s, lu.s, fails)
fit.all <- fit.all[!fit.all$fitness < 0,]
fit.all$rep <- str_extract(fit.all$rep, '\\d+')
fit.s <- fit.all[,c("fitness", "date", "sci", "gene", "type", "treat", "rep")]
fit.s$date <- as.character(fit.s$date)

#write.csv(fit.all, 'results/fitness_all.csv', quote = F, row.names = F)

#explore data before modeling
ggplot(fit.all, aes(x=treat_numeric, y=fitness, color=type, shape=gene)) +
  geom_jitter(width=0.05, size=2, alpha=0.75) +
  theme_classic() +
  scale_colour_colorblind() +
  scale_fill_colorblind()
#ggsave(filename = "plots/all_oneplot.pdf", width = 8, height= 6, units = "in")
#ggsave(filename = "plots/all_oneplot.png", width = 8, height= 6, units = "in")

fit.s.nostab <- fit.s[!str_detect(fit.s$type, "Lo"),]
fit.s.nostab <- fit.s.nostab[!str_detect(fit.s.nostab$type, "Hi"),]
#fit.s.nostab[fit.s.nostab$type == "wt",]$type <- NA


size=3
p.shuf <- ggplot(fit.s.nostab[fit.s.nostab$type %in% c("wt","shuffle"),], aes(x=treat, y=fitness, shape=gene, color=gene, group=gene))+
  geom_jitter(position = position_jitterdodge(dodge.width=0.4, jitter.width=0.1, jitter.height = 0), size=size, alpha=0.75) +
  #facet_grid(type ~., margins = T, scales = "fixed", drop = T) +
  scale_x_discrete(limits=c("0%", "25%", "50%", "75%", "100%", "max")) +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  ggtitle("Shuffle recoding")
p.rand <- ggplot(fit.s.nostab[fit.s.nostab$type %in% c("wt","random"),], aes(x=treat, y=fitness, shape=gene, color=gene, group=gene))+
  geom_jitter(position = position_jitterdodge(dodge.width=0.4, jitter.width=0.1, jitter.height = 0), size=size, alpha=0.75) +
  #facet_grid(type ~., margins = T, scales = "fixed", drop = T) +
  scale_x_discrete(limits=c("0%", "25%", "50%", "75%", "100%", "max")) +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) +
  ggtitle("Random recoding")
p.worst <- ggplot(fit.s.nostab[fit.s.nostab$type %in% c("wt","worst"),], aes(x=treat, y=fitness, shape=gene, color=gene, group=gene))+
  geom_jitter(position = position_jitterdodge(dodge.width=0.4, jitter.width=0.1, jitter.height = 0), size=size, alpha=0.75) +
  #facet_grid(type ~., margins = T, scales = "fixed", drop = T) +
  scale_x_discrete(limits=c("0%", "25%", "50%", "75%", "100%", "max")) +
  theme_bw() +
  ylim(0,25) +
  ggtitle("Worst codon")
plot_grid( p.shuf, p.rand, p.worst, ncol = 1)


fit.s.stab <- fit.s[fit.s$type %in% c("wt", "shuffle", "shuffle_Hi", "shuffle_Lo"),]
fit.s.stab <- fit.s.stab[fit.s.stab$treat %in% c("0%","50%"),]
p.shuff <- ggplot(fit.s.stab, aes(x=type, y=fitness, shape=gene, color=gene, group=gene))+
  geom_jitter(position = position_jitterdodge(dodge.width=0.4, jitter.width=0.1, jitter.height = 0), size=size, alpha=0.75) +
  #facet_grid(type ~., margins = T, scales = "fixed", drop = T) +
  scale_x_discrete(limits=c("wt", "shuffle_Lo", "shuffle", "shuffle_Hi")) +
  theme_bw() +
  ylim(0,23) +
  ggtitle("Stability")

plot_grid(p.shuf, p.rand, p.worst, p.shuff, ncol = 1)
#ggsave(filename = "plots/all.pdf", width = 8.5, height= 11, units = "in")
#ggsave(filename = "plots/all.png", width = 8.5, height= 11, units = "in")

fit.nofail <- fit.all[which(!fit.all$fitness==0),] ##dropping constructs that couldn't be made 
###calculate means and stuff
fit.summary <-  fit.nofail %>% 
  group_by(gene, type, treat_numeric) %>%
  summarise(mean = mean(fitness),
            n = n(),
            sd = sd(fitness),
            se = sd/sqrt(n))
#write.csv(fit.summary, 'results/fit_summary.csv', quote = F, row.names = F)

fit.summary.tech <-  fit.nofail %>% 
  group_by(gene, type, treat_numeric, rep) %>%
  summarise(mean = mean(fitness),
            n = n(),
            sd = sd(fitness),
            se = sd/sqrt(n))
#write.csv(fit.summary.tech, 'results/fit_summary_techreps.csv', quote = F, row.names = F)

ggplot(fit.summary, aes(x=mean, y=sd, fill=gene, color=gene, shape=gene)) +
  geom_point(alpha=0.75) +
  theme_bw() +
  scale_color_jama() +
  xlab('Mean fitness') +
  ylab('Fitness SD')
#ggsave(filename = "plots/errorByfitness.pdf", width = 4.5, height= 4, units = "in")
#ggsave(filename = "plots/errorByfitness.png", width = 4.5, height= 4, units = "in")

add <- fit.all[which(fit.all$fitness == 0),
               c('gene', 'type', 'treat_numeric', 'fitness')]
add$n <- 1 
add$sd <- 0  
add$se <- 0


ggplot() +
  geom_pointrange(data=fit.summary, 
                  aes(x=treat_numeric, 
                      y=mean, 
                      color=type, 
                      shape=gene, 
                      ymin=mean-sd, 
                      ymax=mean+sd), 
                  alpha=0.75, 
                  position=position_jitter(width=0.07)) +
  geom_point(data=add, aes(x=treat_numeric, y=fitness, color=type, shape=gene),
             alpha=0.75, 
             position=position_jitter(width=0.07),
             size=2) +
  theme_classic() +
  scale_colour_colorblind() +
  scale_fill_colorblind() +
  scale_shape_manual(values = c(15,16,17,18)) +
  xlab("Proportion of codons changed") +
  ylab('Fitness (doublings/hr)')
#ggsave(filename = "plots/all_summary_sd.pdf", width = 6, height= 4.5, units = "in")
#ggsave(filename = "plots/all_summary_sd.png", width = 6, height= 4.5, units = "in")


#figure out how to look at assumptions for model
ggplot(fit.nofail[!fit.nofail$gene=="wt",], aes(fitness, fill=treat)) +
  geom_histogram(bins=35) +
  facet_grid(treat ~ ., margins = T, scales = "free") +
  theme_bw() +
  scale_colour_colorblind() +
  scale_fill_colorblind()
#ggsave(filename = "plots/check_dists.pdf", width = 6.5, height= 6.5, units = "in")
#ggsave(filename = "plots/check_dists.png", width = 6.5, height= 6.5, units = "in")

with(fit.nofail[!fit.nofail$gene=="wt",], tapply(fitness, treat, function(x) {
  sprintf("M (SD) = %1.2f (%1.2f)", mean(x), sd(x))
}))
#means are greater than errors

fit.all$date <- as.character(fit.all$date)
fit.all$treatment <- paste(fit.all$gene, fit.all$type, fit.all$treat, fit.all$rep, sep="_")
#fit.all$treatment <- relevel(fit.all$treatment, ref='wt_wt_0%_1')
mod.1 <- glm(fitness ~ treatment, data=fit.all)
mod.2 <- glm(fitness ~ treatment * sci, data=fit.all)
AIC(mod.1)
AIC(mod.2)
logLik(mod.1)
logLik(mod.2)
anova(mod.1, mod.2)
###don't worry about person performing experiment

#mod.date <- glm(fitness ~ )
#date.coeff <- as.matrix(mod.2$coefficients)



###calculate means and stuff for technical replicates
fit.summary.techs <-  fit.nofail %>% 
  group_by(gene, type, treat_numeric, rep) %>%
  summarise(mean = mean(fitness),
            n = n(),
            sd = sd(fitness),
            se = sd/sqrt(n))

#fit.s.nostab[fit.s.nostab$type %in% c("wt","shuffle"),]

size=3
fit.summary.techs$treat <- paste(fit.summary.techs$treat_numeric*100, '%', sep = '')
fit.summary.techs$treat <- str_replace(fit.summary.techs$treat, '125%', 'max')
add$treat <- paste(add$treat_numeric*100, '%', sep = '')
add$treat <- str_replace(add$treat, '125%', 'max')

#panel plots
{
p.shuf <- ggplot() +
  geom_pointrange(data=fit.summary.techs[fit.summary.techs$type %in% c('wt', 'shuffle'),], 
                  aes(x=treat, 
                      y=mean, 
                      color=gene, 
                      shape=gene, 
                      ymin=mean-sd, 
                      ymax=mean+sd), 
                  alpha=0.75, 
                  position=position_jitter(width=0.07)) +
  geom_point(data=add[add$type %in% c('wt','shuffle'),], aes(x=treat, y=fitness, color=gene, shape=gene),
             alpha=0.75, 
             position=position_jitter(width=0.07, height=0),
             size=2) +
  theme_classic() +
  scale_colour_colorblind() +
  scale_fill_colorblind() +
  scale_shape_manual(values = c(15,16,17,18)) +
  scale_x_discrete(limits=c("0%", "25%", "50%", "75%", "100%")) +
  theme(panel.grid.major.y = element_line(color = 'grey', size=0.1)) +
  ggtitle("Shuffle recoding") +
  labs(y='Fitness', x='Level of recoding')

p.rand <- ggplot() +
    geom_pointrange(data=fit.summary.techs[fit.summary.techs$type %in% c('wt', 'random'),], 
                    aes(x=treat, 
                        y=mean, 
                        color=gene, 
                        shape=gene, 
                        ymin=mean-sd, 
                        ymax=mean+sd), 
                    alpha=0.75, 
                    position=position_jitter(width=0.07)) +
    geom_point(data=add[add$type %in% c('wt','random'),], aes(x=treat, y=fitness, color=gene, shape=gene),
               alpha=0.75, 
               position=position_jitter(width=0.07, height=0),
               size=2) +
    theme_classic() +
    scale_colour_colorblind() +
    scale_fill_colorblind() +
    scale_shape_manual(values = c(15,16,17,18)) +
    scale_x_discrete(limits=c("0%", "25%", "50%", "75%", "100%", 'max'), labels=c("0%", "25%", "50%", "75%", "100%", 'all codons')) +
    theme(panel.grid.major.y = element_line(color = 'grey', size = 0.1)) +
    ggtitle("Random recoding") +
    labs(y='Fitness', x='Level of recoding')
  
p.worst <- ggplot() +
  geom_pointrange(data=fit.summary.techs[fit.summary.techs$type %in% c('wt', 'worst'),], 
                  aes(x=treat, 
                      y=mean, 
                      color=gene, 
                      shape=gene, 
                      ymin=mean-sd, 
                      ymax=mean+sd), 
                  alpha=0.75, 
                  position=position_jitter(width=0.07)) +
  geom_point(data=add[add$type %in% c('wt','worst'),], aes(x=treat, y=fitness, color=gene, shape=gene),
             alpha=0.75, 
             position=position_jitter(width=0.07, height=0),
             size=2) +
  theme_classic() +
  scale_colour_colorblind() +
  scale_fill_colorblind() +
  scale_shape_manual(values = c(15,16,17,18)) +
  scale_x_discrete(limits=c("0%", "25%", "50%", "75%", "100%")) +
  theme(panel.grid.major.y = element_line(color = 'grey', size=0.1)) +
  ggtitle("Worst recoding") +
  labs(y='Fitness', x='Level of recoding')

fit.sum.stab <- fit.summary.techs[fit.summary.techs$type %in% c("wt", "shuffle", "shuffle_Hi", "shuffle_Lo"),]
fit.sum.stab <- fit.sum.stab[fit.sum.stab$treat %in% c("0%","50%"),]
add.stab <- add[add$type %in% c("wt", "shuffle", "shuffle_Hi", "shuffle_Lo"),]
add.stab <- add.stab[add.stab$treat %in% c("0%","50%"),]

p.stab <- ggplot() +
  geom_pointrange(data=fit.sum.stab, 
                  aes(x=type, 
                      y=mean, 
                      color=gene, 
                      shape=gene, 
                      ymin=mean-sd, 
                      ymax=mean+sd), 
                  alpha=0.75, 
                  position=position_jitter(width=0.07)) +
  geom_point(data=add.stab, aes(x=type, y=fitness, color=gene, shape=gene),
             alpha=0.75, 
             position=position_jitter(width=0.07, height=0),
             size=2) +
  theme_classic() +
  scale_colour_colorblind() +
  scale_fill_colorblind() +
  scale_shape_manual(values = c(15,16,17,18)) +
  scale_x_discrete(limits=c("shuffle_Lo", "shuffle", "shuffle_Hi"), labels=c("Low stability", "50% shuffle", "Hi stability")) +
  theme(panel.grid.major.y = element_line(color = 'grey', size=0.1)) +
  ggtitle("Stability") +
  labs(y='Fitness', x='Recoding type')

plot_grid(p.shuf, p.rand, p.worst, p.stab, ncol = 2)
#ggsave(filename = "plots/all_techreps.pdf",  width = 7.5, height= 6.5, units = "in")
#ggsave(filename = "plots/all_techreps.png", width = 7.5, height= 6.5, units = "in")


ggplot() +
  geom_hline(yintercept = 20.35) +
  geom_pointrange(data=fit.summary.techs[fit.summary.techs$type %in% c('shuffle', 'random', 'worst', 'shuffle_Hi', 'shuffle_Lo'),], 
                  aes(x=treat, 
                      y=mean, 
                      color=gene, 
                      shape=gene, 
                      ymin=mean-sd, 
                      ymax=mean+sd,
                      group=type), 
                  alpha=0.75, 
                  position=position_jitterdodge(jitter.width=0.05, dodge.width=0.7)) +
  geom_point(data=add[add$type %in% c('shuffle', 'random','worst', 'shuffle_Hi', 'shuffle_Lo'),], aes(x=treat, y=fitness, color=gene, shape=gene, group=type),
             alpha=0.75, 
             position=position_jitterdodge(jitter.width=0.05, dodge.width=0.7),
             size=2) +
  theme_classic() +
  scale_colour_colorblind() +
  scale_fill_colorblind() +
  scale_shape_manual(values = c(15,16,17,18)) +
  scale_x_discrete(limits=c("25%", "50%", "75%", "100%")) +   #####error gets thrown bc max not in here
  theme(panel.grid.major.y = element_line(color = 'grey', size=0.1)) +
  #ggtitle("Worst recoding") +
  labs(y='Fitness', x='Level of recoding') 

ggsave(filename = "plots/all_techreps_oneplot.pdf",  width = 7, height= 4.5, units = "in")
ggsave(filename = "plots/all_techreps_oneplot.png", width = 7, height= 4.5, units = "in")

ggplot(fit.s.nostab, 
       aes(y=fitness, 
           fill=gene, 
           group=gene,
           alpha=treat,)) + 
  geom_density(color=NA)
}




##################################################process combo data

############################copied data files over into Rproject data folder!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
############################update if files are changed!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
combo_sum <- read_xlsx("data/FGH combo titer_fitness data.xlsx", sheet = 5, skip = 1, col_names = T, col_types = "text")
names(combo_sum)[6:24] <- c("14-May-21",	"20-May-21",	"2-Jun-21",	"4-Jun-21",	"8-Jun-21",	"9-Jun-21",	"10-Jun-21",	"11-Jun-21",	"16-Jun-21",	"17-Jun-21",	"18-Jun-21",	"22-Jun-21",	"23-Jun-21",	"24-Jun-21",	"30-Jun-21",	"1-Jul-21",	"2-Jul-21",	"7-Jul-21",	"8-Jul-21")
names(combo_sum)[29] <- "notes"
names(combo_sum)[1] <- "sample"
names(combo_sum) <- str_remove_all(names(combo_sum), "\\.")
combo_sum <- combo_sum[,-25]
combo_long <- melt(combo_sum, id.vars = c("sample", "construct", "F", "G", "H", "variant avg", "S2M avg", "predicted fitness", "notes"))
#combo_long[which(str_detect(combo_long$notes, "construction failed")),]$value <- "0"
combo_long <- combo_long[!is.na(combo_long$value),]

###have to fix all the formatting issues
##more formatting
types <- c("random", "shuffle", "worst", "S2M")
amounts <- c("100%", "75%", "50%", "25%", "0.75", "max", "0\\.5", "0\\.25", "1\\.0")

combo_long$`G` <- str_replace(combo_long$`G`, 'shuf_', 'shuffle_')
combo_long$`F` <- str_replace_all(combo_long$`F`, " ", "_")
combo_long$`G` <- str_replace_all(combo_long$`G`, " ", "_")
combo_long$`H` <- str_replace_all(combo_long$`H`, " ", "_")
combo_long$`F` <- str_replace(combo_long$`F`, "rep_", "rep")
combo_long$`G` <- str_replace(combo_long$`G`, "rep_", "rep")
combo_long$`H` <- str_replace(combo_long$`H`, "rep_", "rep")
combo_long$F_name <- NA
combo_long$G_name <- NA
combo_long$H_name <- NA

for(i in 1:nrow(combo_long)){
  F <- combo_long[i,]$`F`
  G <- combo_long[i,]$`G`
  H <- combo_long[i,]$`H`
  ftype <- as.character(na.omit(str_extract(F, types)))
  gtype <- as.character(na.omit(str_extract(G, types)))
  htype <- as.character(na.omit(str_extract(H, types)))
  flev <- as.character(na.omit(str_extract(F, amounts)))
  glev <- as.character(na.omit(str_extract(G, amounts)))
  hlev <- as.character(na.omit(str_extract(H, amounts)))
  frep <- as.character(na.omit(str_extract(F, 'rep\\d')))
  grep <- as.character(na.omit(str_extract(G, 'rep\\d')))
  hrep <- as.character(na.omit(str_extract(H, 'rep\\d')))
  combo_long[i,]$F_name <- paste('F', ftype, flev, frep, sep="_")
  combo_long[i,]$G_name <- paste('G', gtype, glev, grep, sep="_")
  combo_long[i,]$H_name <- paste("H", htype, hlev, hrep, sep="_")
}

combo_long$F_name <- str_replace(combo_long$F_name, '25%', '0.25') 
combo_long$F_name <- str_replace(combo_long$F_name, '50%', '0.5') 
combo_long$F_name <- str_replace(combo_long$F_name, '75%', '0.75') 
combo_long$F_name <- str_replace(combo_long$F_name, '100%', '1.0') 

combo_long$G_name <- str_replace(combo_long$G_name, '25%', '0.25') 
combo_long$G_name <- str_replace(combo_long$G_name, '50%', '0.5') 
combo_long$G_name <- str_replace(combo_long$G_name, '75%', '0.75') 
combo_long$G_name <- str_replace(combo_long$G_name, '100%', '1.0') 

combo_long$H_name <- str_replace(combo_long$H_name, '25%', '0.25') 
combo_long$H_name <- str_replace(combo_long$H_name, '50%', '0.5') 
combo_long$H_name <- str_replace(combo_long$H_name, '75%', '0.75') 
combo_long$H_name <- str_replace(combo_long$H_name, '100%', '1.0') 

combo_long[which(str_detect(combo_long$F_name, "S2M")),]$F_name <- 'wt'
combo_long[which(str_detect(combo_long$G_name, "S2M")),]$G_name <- 'wt'
combo_long[which(str_detect(combo_long$H_name, "S2M")),]$H_name <- 'wt'

names(combo_long) <- str_replace_all(names(combo_long), " ", "_")
combo_long$value <- as.numeric(combo_long$value)
combo_long$predicted_fitness <- as.numeric(combo_long$predicted_fitness)
combo_long$S2M_avg <- as.numeric(combo_long$S2M_avg)
names(combo_long) <- str_replace_all(names(combo_long), " ", "_")
fit.summary.combo.tech <-  combo_long %>%
  select(c('construct','S2M_avg','predicted_fitness', 'value', 'F_name', 'G_name', 'H_name')) %>%
  group_by(construct) %>%
  summarise(mean = mean(value),
            expect = mean(predicted_fitness),
            F_name = unique(F_name),
            G_name = unique(G_name),
            H_name = unique(H_name),
            S2M = mean(S2M_avg),
            n = n(),
            sd = sd(value),
            se = sd/sqrt(n))
fit.summary.combo.tech$set <- str_split(fit.summary.combo.tech$construct, " ", simplify = T)[,1]
fit.summary.combo.tech$ngenes <- str_length(str_split(fit.summary.combo.tech$construct, " ", simplify = T)[,2])-2

fit.summary.tech <- read.csv('results/fit_summary_techreps.csv')
fit.summary.tech$treat_numeric <- as.character(fit.summary.tech$treat_numeric)
fit.summary.tech[which(fit.summary.tech$treat_numeric == 1), ]$treat_numeric <- "1.0" 
row.names(fit.summary.tech) <- paste(fit.summary.tech$gene, "_",
                                     fit.summary.tech$type, "_",
                                     fit.summary.tech$treat_numeric, "_",
                                     "rep_",
                                     fit.summary.tech$rep, sep="")
row.names(fit.summary.tech) <- str_replace(row.names(fit.summary.tech), 'rep_', 'rep')
row.names(fit.summary.tech) <- str_replace(row.names(fit.summary.tech), '1\\.25', 'max')
fit.summary.tech <- fit.summary.tech[!fit.summary.tech$gene=='wt',]

##add the failed triples
failed <- data.frame(construct=c("Rand-1 FGH", "Rand-2 FGH", "Worst-2 FGH", "Worst-3 FGH"), 
                     expect=c(6.04, 7.49, 6.07, 4.35),
                     gene=c("FGH", "FGH", "FGH", "FGH"),
                     ngenes=c(3,3,3,3),
                     set=c("Rand-1", "Rand-2", "Worst-2", "Worst-3"),
                     n=c(3,3,3,3),
                     mean=c(0,0,0,0),
                     sd=c(0,0,0,0),
                     se=c(0,0,0,0))
row.names(failed) <- failed$set

sets <- unique(fit.summary.combo.tech$set)
sets <- sets[-10]  ##drop wt
plots <- list()
for(j in 1:length(sets)){
  i <- sets[j]
  tmp <- fit.summary.combo.tech[fit.summary.combo.tech$set == i,]
  singles <- na.omit(fit.summary.tech[unique(unlist(tmp[,c('F_name', 'G_name', 'H_name')])),])
  tmp$gene <- str_split(tmp$construct, " ", simplify = T)[,2]
  tmp$gene <- str_split(tmp$gene, "-", simplify = T)[,1]
  singles$ngenes <- 1
  singles$construct <- row.names(singles)
  singles$expect <- 50
  singles$set <- i
  triple <- failed[i,]
  triple$ngenes <- 3
  tmp2 <- tmp[,c("construct", "set", "gene", "ngenes", "expect", "mean", "sd", "se", "n")] 
  singles2 <- singles[,c("construct", "set", "gene", "ngenes","expect", "mean", "sd", "se", "n")]
  triple <- triple[,c("construct", "set", "gene", "ngenes","expect", "mean", "sd", "se", "n")]
  join <- rbind(singles2, tmp2, triple)
  join$ngenes <- as.character(join$ngenes)
  p <- ggplot() +
    geom_pointrange(data=join, 
                    aes(x=ngenes, 
                        y=mean, 
                        ymin=mean-sd, 
                        ymax=mean+sd), 
                    alpha=0.75, 
                    position=position_jitter(width=0.1)) +
    geom_point(data=join, aes(x=ngenes, y=expect, color='red', label=gene), 
               shape=17,
               alpha=0.75, 
               position=position_jitter(width=0.07, height=0),
               size=3) +
    geom_text(data=join, aes(label=gene, x=ngenes, y=mean), hjust=-1/as.numeric(join$ngenes), size=3) +
    geom_text(data=join, aes(label=gene, x=ngenes, y=expect), hjust=-1/as.numeric(join$ngenes), size=3) +
    theme_classic() +
    theme(panel.grid.major.y = element_line(color = 'grey', size=0.1), legend.position = 'none', plot.margin = unit(c(5.5,5.5,5.5,5.5),"pt")) +
    ggtitle(i) +
    labs(y='Fitness', x='Number of recoded genes') +
    ylim(c(-0.1,21))
  print(j)
  print(p)
  plots[[j]] = p  
}

pfinal <- plot_grid(plots[[1]] + labs(x='', y=''), 
                    plots[[2]] + labs(x='', y=''), 
                    plots[[3]] + labs(x='', y=''), 
                    plots[[4]] + labs(x=''), 
                    plots[[5]] + labs(x='', y=''), 
                    plots[[6]] + labs(x='', y=''), 
                    plots[[7]] + labs(x='', y=''), 
                    plots[[8]] + labs(y=''), 
                    plots[[9]] + labs(x='', y=''), 
                    nrow = 3)
pfinal
#ggsave(pfinal, filename = "plots/combos_crm2.pdf", width = 10, height= 10.5, units = "in")
#ggsave(pfinal, filename = "plots/combos.png", width = 10, height= 10.5, units = "in")



###################################################################################add recovery numbers
##format:
#gene_type_amount_rep
#G_random_0.5_rep1

rec_lu <- read_xlsx("data/FGH Passaging Data Lu 8NOv21.xlsx", sheet = "SUMMARY")
names(rec_lu) <- str_replace_all(names(rec_lu), " ", "_")
names(rec_lu) <- str_replace_all(names(rec_lu), "\\#", "num")
names(rec_lu) <- str_replace_all(names(rec_lu), "_\\(Y\\/N\\)", "")
names(rec_lu) <- str_replace_all(names(rec_lu), "_\\(hrs\\)", "")

rec_lu <- rec_lu[!is.na(rec_lu$T40_Titer),] 
rec_lu <- rec_lu[-which(rec_lu$Construct=="Construct"),] 

for(i in 1:nrow(rec_lu)){
  if(is.na(rec_lu[i,]$Sample_num)){
    rec_lu[i,1:3] <- fix
  }else{
    fix <- rec_lu[i,1:3]
  }
}
rec_lu$Date <- openxlsx::convertToDate(rec_lu$Date)
rec_lu$Fitness_Lrg_plqs <- NA
rec_lu$Fitness_Sm_plqs <- NA
fix <- which(str_detect(rec_lu$Construct, 'Worst.*rep'))
rec_lu[fix,]$Construct <- str_split(rec_lu[fix,]$Construct, '_', simplify = T)[,2]
fix <- which(str_detect(rec_lu$Construct, 'Worst.*Rep'))
rec_lu[fix,]$Construct <- str_split(rec_lu[fix,]$Construct, '_', simplify = T)[,2]
fix <- which(str_detect(rec_lu$Construct, 'Wosrt.*rep'))
rec_lu[fix,]$Construct <- str_split(rec_lu[fix,]$Construct, '_', simplify = T)[,2]

rec_emma <- read_xlsx("data/FGH Passaging Data Emma.xlsx", sheet='Summary')
names(rec_emma) <- str_replace_all(names(rec_emma), " ", "_")
names(rec_emma) <- str_replace_all(names(rec_emma), "\\?", "")
names(rec_emma) <- str_replace_all(names(rec_emma), "\\#", "num")
names(rec_emma) <- str_replace_all(names(rec_emma), "_\\(Y\\/N\\)", "")
names(rec_emma) <- str_replace_all(names(rec_emma), "_\\(hrs\\)", "")
rec_emma <- rec_emma[!is.na(rec_emma$Passage_Fitness),]
names(rec_emma) <- str_replace(names(rec_emma), "Large_Plqs", "Mutants")

#important columns
get <- c('Sample_num', 'Construct', "Date", 'Passage_num', 'Original_Fitness', 'Passage_Fitness', 'Passage_Time', "Mutants", 'Notes_on_Mutants', 'Fitness_Lrg_plqs', 'Fitness_Sm_plqs')

rec_tessa <- read_xlsx("data/FGH Passaging Data Tessa.xlsx", sheet='Organized passages', skip = 8)
names(rec_tessa) <- str_replace_all(names(rec_tessa), " ", "_")
names(rec_tessa) <- str_replace_all(names(rec_tessa), "\\?", "")
names(rec_tessa) <- str_replace_all(names(rec_tessa), "\\#", "num")
names(rec_tessa) <- str_replace_all(names(rec_tessa), "_\\(Y\\/N\\)", "")
names(rec_tessa) <- str_replace_all(names(rec_tessa), "_\\(hrs\\)", "")
names(rec_tessa) <- str_replace(names(rec_tessa), "sample", "Sample")
names(rec_tessa) <- str_replace(names(rec_tessa), "comment", "Notes_on_Mutants")
names(rec_tessa) <- str_replace(names(rec_tessa), "Time\\(Hr.Min\\)", "Passage_Time")
names(rec_tessa) <- str_replace(names(rec_tessa), "sample", "Sample")
names(rec_tessa) <- str_replace(names(rec_tessa), "fitness_comb", "Passage_Fitness")
rec_tessa <- rec_tessa[!is.na(rec_tessa$Passage_Fitness),]
get[!get %in% names(rec_tessa)]
tmp <- round(as.numeric(str_remove(rec_tessa$Sample_num, 'R')), digits = 4)
rec_tessa$Passage_num <- as.numeric(str_split(tmp, '\\.', simplify = T)[,2])
rec_tessa[is.na(rec_tessa$Passage_num),]$Passage_num <- 0
rec_tessa$Passage_num <- rec_tessa$Passage_num + 1
rec_tessa$Original_Fitness <- NA
rec_tessa$Mutants <- NA
rec_tessa$Fitness_Lrg_plqs <- NA
rec_tessa$Fitness_Sm_plqs <- NA
rec_tessa <- rec_tessa[,get]
rec_tessa$sci <- 'tessa'

rec_jacob <- read_xlsx("data/FGH Passaging Data Jacob.xlsx", sheet='Summary')
names(rec_jacob) <- str_replace_all(names(rec_jacob), " ", "_")
names(rec_jacob) <- str_replace_all(names(rec_jacob), "\\.", "")
names(rec_jacob) <- str_replace_all(names(rec_jacob), "\\#", "num")
names(rec_jacob) <- str_replace_all(names(rec_jacob), "_\\(Y\\/N\\)", "")
names(rec_jacob) <- str_replace_all(names(rec_jacob), "_\\(hrs\\)", "")
names(rec_jacob) <- str_replace(names(rec_jacob), "Number", "Sample_num")
names(rec_jacob) <- str_replace(names(rec_jacob), "Passage_Data_Fitness", "Passage_Fitness")
names(rec_jacob) <- str_replace(names(rec_jacob), "Large_Plaque_Fitness_\\(Best_Plate\\)", "Fitness_Lrg_plqs")
names(rec_tessa) <- str_replace(names(rec_tessa), "sample", "Sample")
get[!get %in% names(rec_jacob)]
rec_jacob$Original_Fitness <- NA
rec_jacob$Fitness_Sm_plqs <- NA
rec_jacob <- rec_jacob[,get]
rec_jacob <- rec_jacob[!is.na(rec_jacob$Construct),]
rec_jacob$sci <- 'jacob'


rec_lu <- rec_lu[,get]
rec_lu$sci <- 'lu'
rec_emma <- rec_emma[,get]
rec_emma$sci <- 'emma'

recov_all <- rbind(rec_lu, rec_emma, rec_tessa)  ###droped Jacobs
##lu added the set name to the singles. this would have been really helpful if everyone would have
recov_all[recov_all$ngenes==1,]$Construct <- str_replace(recov_all[recov_all$ngenes==1,]$Construct, "Worst-2_", "")
recov_all[recov_all$ngenes==1,]$Construct <- str_replace(recov_all[recov_all$ngenes==1,]$Construct, "Worst-1_", "")
recov_all[recov_all$ngenes==1,]$Construct <- str_replace(recov_all[recov_all$ngenes==1,]$Construct, "Wosrt-1_", "")
recov_all$Construct <- str_replace(recov_all$Construct, "Worst-1_", "Worst-1 ") ##lu changed naming convention
recov_all$Construct <- str_replace(recov_all$Construct, "Worst-2_", "Worst-2 ") ##lu changed naming convention
recov_all$Construct <- str_replace(recov_all$Construct, "Worst-3_", "Worst-3 ") ##lu changed naming convention
recov_all$Construct <- str_replace(recov_all$Construct, "Shuf-3-", "Shuf-3 ")  ##lu changed naming convention
recov_all$Construct <- str_replace(recov_all$Construct, "G_0.5.Shuf", "G 50% shuffle ") ##lu changed naming convention
recov_all$Construct <- str_replace(recov_all$Construct, "Shuffle 100%", "100% shuffle")  ##lu changed naming convention
recov_all$Construct <- str_replace(recov_all$Construct, "F_", "F ")  ##lu changed naming convention
recov_all$Construct <- str_replace(recov_all$Construct, "G_", "G ")  ##lu changed naming convention
recov_all$Construct <- str_replace(recov_all$Construct, "Frandom", "F random") 
recov_all$Construct <- str_replace(recov_all$Construct, " _rep", " rep") 
recov_all$Construct <- str_replace(recov_all$Construct, "t_rep", "t rep")
recov_all$Construct <- str_replace(recov_all$Construct, "Rep", "rep")
recov_all$Construct <- str_replace(recov_all$Construct, "rep ", " rep") 
recov_all$Construct <- str_replace(recov_all$Construct, "100% shuffle", "shuffle 100%")
recov_all$Construct <- str_replace(recov_all$Construct, "100% worst", "worst 100%")
recov_all$Construct <- str_replace(recov_all$Construct, "75% worst", "worst 75%")
recov_all$Construct <- str_replace(recov_all$Construct, "50% worst", "worst 50%")
recov_all$Construct <- str_replace(recov_all$Construct, "50% shuffle", "shuffle 50%")
recov_all$Construct <- str_replace(recov_all$Construct, "  ", " ")
############need to separate single recoveries from combos, or at least label them
recov_all$ngenes <- str_detect(recov_all$Construct, "rep")
recov_all$ngenes <- as.numeric(recov_all$ngenes)


rec_singles <- recov_all[recov_all$ngenes==1,]
rec_combos <- recov_all[!recov_all$ngenes==1,]
rec_combos$set <- str_split(rec_combos$Construct, " ", simplify = T)[,1]
rec_combos$gene <- str_split(rec_combos$Construct, " ", simplify = T)[,2]
rec_combos$gene <- str_split(rec_combos$gene, "-", simplify = T)[,1]
rec_combos$ngenes <- str_length(rec_combos$gene)
rec_combos$Passage_num <- as.numeric(rec_combos$Passage_num) + 3
rec_combos$Passage_Fitness <- as.numeric(rec_combos$Passage_Fitness)
rec_combos <- rec_combos[!is.na(rec_combos$Passage_Fitness),]

fit.summary.rec.combos <-  rec_combos %>%
  select(c('Construct','Passage_num','Passage_Fitness', 'ngenes', 'set', 'gene')) %>%
  group_by(Construct, Passage_num) %>%
  summarise(mean = mean(Passage_Fitness),
            Construct = unique(Construct),
            Passage_num = unique(Passage_num),
            set = unique(set),
            gene = unique(gene),
            n = n(),
            sd = sd(Passage_Fitness),
            se = sd/sqrt(n))


###have to fix all the formatting issues
##more formatting
types <- c("random", "shuffle", "worst", "S2M")
amounts <- c("100%", "75%", "50%", "25%", "0.75", "max", "0\\.5", "0\\.25", "1\\.0")

rec_singles$Construct <- str_replace(combo_long$rec_singles$Construct, 'shuf_', 'shuffle_')
rec_singles$Construct <- str_replace_all(rec_singles$Construct, " ", "_")
rec_singles$Construct <- str_replace(rec_singles$Construct, "rep_", "rep")

for(i in 1:nrow(rec_singles)){
  tmp <- rec_singles[i,]$Construct
  gene <- str_split(tmp, "_", simplify = T)[1]
  type <- as.character(na.omit(str_extract(tmp, types)))
  lev <- as.character(na.omit(str_extract(tmp, amounts)))
  rep <- as.character(na.omit(str_extract(tmp, 'rep\\d')))
  rec_singles[i,]$Construct <- paste(gene, type, lev, rep, sep="_")
}

rec_singles$Construct <- str_replace(rec_singles$Construct, '25%', '0.25') 
rec_singles$Construct <- str_replace(rec_singles$Construct, '50%', '0.5') 
rec_singles$Construct <- str_replace(rec_singles$Construct, '75%', '0.75') 
rec_singles$Construct <- str_replace(rec_singles$Construct, '100%', '1.0') 
##fixing a couple of things
rec_singles$Construct <- str_replace(rec_singles$Construct, 'Frandom', 'F')
rec_singles$Construct <- str_replace(rec_singles$Construct, '__', '_shuffle_')
rec_singles$set <- "NA"
rec_singles$gene <- str_split(rec_singles$Construct, "_", simplify = T)[,1]

###added some singles that are not part of a set
##must pull them out
rec_singles_no_set <- rec_singles[!rec_singles$Construct %in% unique(c(fit.summary.combo.tech$F_name,fit.summary.combo.tech$G_name,fit.summary.combo.tech$H_name)),]
rec_singles <- rec_singles[rec_singles$Construct %in% unique(c(fit.summary.combo.tech$F_name,fit.summary.combo.tech$G_name,fit.summary.combo.tech$H_name)),]


for(i in 1:nrow(rec_singles)){
  gene <- str_split(rec_singles[i,]$Construct, "_", simplify = T)[,1]
  find <- rec_singles[i,]$Construct
  search <- fit.summary.combo.tech[,which(names(fit.summary.combo.tech) == paste(gene, "_name", sep=""))]
  set <- unique(fit.summary.combo.tech[which(search==find),]$set)
  if(length(set) > 1){
    print("Problem. Found more than one set for this fragment")
  }else{
    rec_singles[i,]$set <- set
  }
}
rec_singles$Passage_Fitness <- as.numeric(rec_singles$Passage_Fitness)
rec_singles <- rec_singles[!is.na(rec_singles$Passage_Fitness),]
fit.summary.rec.singles <-  rec_singles %>%
  select(c('Construct','Passage_num','Passage_Fitness', 'ngenes', 'set', 'gene')) %>%
  group_by(Construct, Passage_num) %>%
  summarise(mean = mean(Passage_Fitness),
            Construct = unique(Construct),
            Passage_num = unique(Passage_num),
            set = unique(set),
            gene = unique(gene),
            n = n(),
            sd = sd(Passage_Fitness),
            se = sd/sqrt(n))

#write.csv(rec_combos, file="results/recovery_gene_combos.csv", quote = FALSE, row.names = FALSE)
#write.csv(rec_singles, file="results/recovery_gene_singles.csv", quote = FALSE, row.names = FALSE)




plots2 <- list()
for(j in 1:length(sets)){
  i <- sets[j]
  tmp <- fit.summary.combo.tech[fit.summary.combo.tech$set == i,]
  singles <- na.omit(fit.summary.tech[unique(unlist(tmp[,c('F_name', 'G_name', 'H_name')])),])
  tmp$gene <- str_split(tmp$construct, " ", simplify = T)[,2]
  tmp$gene <- str_split(tmp$gene, "-", simplify = T)[,1]
  singles$ngenes <- 1
  singles$construct <- row.names(singles)
  singles$expect <- 50
  singles$set <- i
  triple <- failed[i,]
  triple$ngenes <- 3
  tmp2 <- tmp[,c("construct", "set", "gene", "ngenes", "expect", "mean", "sd", "se", "n")] 
  singles2 <- singles[,c("construct", "set", "gene", "ngenes","expect", "mean", "sd", "se", "n")]
  triple <- triple[,c("construct", "set", "gene", "ngenes","expect", "mean", "sd", "se", "n")]
  join <- rbind(singles2, tmp2, triple)
  #join$ngenes <- as.numeric(join$ngenes)
  join$ngenes <- as.character(join$ngenes)
  
  ###add recovery stuff
  #tmp_rec <- rec_combos[rec_combos$set==i,]
  #tmp_rec <- tmp_rec[!tmp_rec$sci=="jacob",]  ###does not work with summary data
  tmp_rec <- fit.summary.rec.combos[fit.summary.rec.combos$set==i,]
  tmp_rec$Passage_num <- paste("p", as.character(tmp_rec$Passage_num-3), sep="")
  
  tmp_rec_singles <- fit.summary.rec.singles[fit.summary.rec.singles$set==i,]
  #tmp_rec_singles <- rec_singles[rec_singles$set==i,]
  #tmp_rec_singles <- tmp_rec_singles[!tmp_rec_singles$sci=="jacob",]
  tmp_rec_singles$Passage_num <- paste("p", as.character(tmp_rec_singles$Passage_num), sep="")
  #tmp_rec_singles$Passage_Fitness <- as.numeric(tmp_rec_singles$Passage_Fitness)
  tmp_rec_singles$mean <- as.numeric(tmp_rec_singles$mean)
  
  p <- ggplot() +
  geom_vline(xintercept = "3", color="grey") +
  geom_pointrange(data=join, 
                  aes(x=ngenes, 
                      y=mean, 
                      ymin=mean-sd, 
                      ymax=mean+sd), 
                  alpha=0.75, 
                  position=position_jitter(width=0.05)) +
  geom_point(data=join, aes(x=ngenes, y=expect),
             color='red',
             shape=17,
             alpha=0.75, 
             position=position_jitter(width=0.05, height=0),
             size=3) +
  geom_point(data=tmp_rec, aes(x=Passage_num, y=mean), 
            size=3,
            color="#2b8cbe",
            position=position_jitter(width=0.05, height=0)) +
  geom_line(data=tmp_rec, aes(x=Passage_num, y=mean, group=gene), 
            color="#2b8cbe") +
  geom_point(data=tmp_rec_singles, aes(x=Passage_num, y=mean), 
            size=3,
            color="#8856a7",
            shape=15,
            position=position_jitter(width=0.05, height=0)) +
  geom_line(data=tmp_rec_singles, aes(x=Passage_num, y=mean, group=gene), 
              color="#8856a7") +
  geom_label_repel(data=join, aes(label=gene, x=ngenes, y=mean), fill = NA, label.size=NA) +
  geom_label_repel(data=join, aes(label=gene, x=ngenes, y=expect), fill = NA, label.size=NA) +
  geom_label_repel(data=tmp_rec, aes(label=gene, group=gene, x=Passage_num, y=mean), fill = NA, label.size=NA) +
  geom_label_repel(data=tmp_rec_singles, aes(label=gene, group=gene, x=Passage_num, y=mean), fill = NA, label.size=NA) +
  theme_classic() +
  theme(panel.grid.major.y = element_line(color = 'grey', size=0.1),
        legend.position = 'none',
        plot.margin=unit(c(0.25,0.75,0.25,0.25),"cm")) +
  ggtitle(i) +
  labs(y='Fitness', x='Number of recoded genes / Passage number') +
  ylim(c(-0.1,21))
  print(j)
  print(p)
  plots2[[j]] = p  
}

plot_grid(plots2[[1]], plots2[[2]], plots2[[3]], plots2[[4]], plots2[[5]], plots2[[6]], plots2[[7]], plots2[[8]], plots2[[9]], nrow = 3)
#ggsave(filename = "plots/combos_rec.pdf", width = 16, height= 14, units = "in")
#ggsave(filename = "plots/combos_rec.png", width = 16, height= 14, units = "in")

