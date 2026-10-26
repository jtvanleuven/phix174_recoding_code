# ============================================================================
# 04_calc_divergence.R
# Codon-divergence tables: for each recoded construct, compute the number
# of codons changed relative to WT and the resulting divergence metrics.
#
# Data dependencies (relative to repo root):
#   * data/fastas/all.fasta       (WT + recoded sequences)
#   * data/design/construct_table.csv
#
# NOTE: as-written pipeline script; not guaranteed turnkey. One F_old2
# Dropbox path is commented out. Paths are relative to the code/ directory.
# ============================================================================

###load phix sequences to calculate distance between variants
library(sangerseqR)
library(stringr)
library(Biostrings)
library(seqinr)
#need to fix Biostrings::translate GENETIC_CODE, otherwise TTG and CTG get translated to methionine
PHIX_GENETIC_CODE <- GENETIC_CODE
attr(PHIX_GENETIC_CODE, "alt_init_codons") <- "ATG"

##NOTE: uses genome from old deoptimized work
vars <- readDNAStringSet("data/phiX_genomes.fa")

pa <- pairwiseAlignment(vars$phiXanc,vars$phiXlow,type="global-local")
writePairwiseAlignments(pa)
missum <- mismatchSummary(pa)
snp.tbl <- missum$pattern$position
View(snp.tbl)
#need to rearrange so that position 1 is start of A (3981)  offset by 1407
snp.tbl.fix <- rbind(snp.tbl[3981:nrow(snp.tbl),],snp.tbl[1:3980,])
#write.csv(snp.tbl.fix, file="phiX_attenuation_snps.csv",row.names = F)

####need to create table 20 rows. cells contain number of codons changed for each amino acid
genes_orig <- readDNAStringSet("data/all.fasta")   ##old sequences from GBE paper
G_new <- readDNAStringSet("data/__G_FullSetToOrder.fasta")
F_new <- readDNAStringSet("data/F_FullSetToOrder_fixSangerBsmb1.fasta")
F_old <- readDNAStringSet("data/__F_FullSetToOrder.fasta")
# F_old2 <- readDNAStringSet("../data/fastas/__F_FullSetToOrder.fasta")  # not distributed; only used in one commented-out line
H_new <- readDNAStringSet("data/__H_FullSetToOrder.fasta")
genes <- c(genes_orig, G_new, F_new, F_old, H_new)
#genes <- c(genes_orig, G_new, F_new, F_old, H_new, F_old2)
tab <- as.data.frame(matrix(ncol=20, nrow=length(genes)))
row.names(tab) <- names(genes)
names(tab) <- sort(unique(PHIX_GENETIC_CODE))[2:21]
for(i in 1:nrow(tab)){
  name <- row.names(tab)[i]
  if(!str_detect(name,"phix")){
    name.s <- str_split(name, "", simplify = T)[1]
  }else{
    name.s <- str_split(name, "", simplify = T)[6]
  }
  wt <- which(str_detect(row.names(tab), paste(name.s,"wt",sep="")))
  seq1 <- genes[wt]
  seq2 <- genes[i]
  cods1 <- substring(seq1, seq(1,nchar(seq1),3), seq(3,nchar(seq1),3))
  cods2 <- substring(seq2, seq(1,nchar(seq2),3), seq(3,nchar(seq2),3))
  aas1 <- c(str_split(as.character(Biostrings::translate(seq1)), "", simplify = T))
  diffs <- table(aas1[!cods1==cods2])
  tab[i,] <- diffs[names(tab[i,])]
}

tab2 <- tab[!str_detect(row.names(tab),"wt"),]
tab2[is.na(tab2)] <- 0
tab3 <- tab2[,colSums(tab2) > 0]
row.names(tab3) <- str_replace(row.names(tab3),"phix_","")
#write.csv(tab3,file="results/diffcounts.csv", quote = F)



#######look at new constructs
###want to quantify unique and shared mutations
tab.g <- as.data.frame(matrix(ncol=176, nrow=length(G_new)))
row.names(tab.g) <- names(G_new)
names(tab.g) <- as.character(1:176)
tab.g.diffs <- tab.g
for(i in 1:nrow(tab.g)){
  name <- row.names(tab.g)[i]
  seq1 <- genes["phix_Gwt"]
  seq2 <- genes[name]
  cods1 <- substring(seq1, seq(1,nchar(seq1),3), seq(3,nchar(seq1),3))
  cods2 <- substring(seq2, seq(1,nchar(seq2),3), seq(3,nchar(seq2),3))
  tab.g[i,] <- cods1==cods2
  tab.g.diffs[i,] <- paste(cods1, cods2, sep = "->")
}
tab.g <- tab.g*1
tab.g$id <- row.names(tab.g)
#write.csv(tab.g, file="results/G_new_diffs.csv", quote = F, row.names = F)

#F_all <- c(F_new,F_old2)
F_all <- c(F_new,F_old[which(!str_detect(names(F_old),"shuffle"))])
tab.f <- as.data.frame(matrix(ncol=428, nrow=length(F_all)))
row.names(tab.f) <- names(F_all)
names(tab.f) <- as.character(1:428)
tab.f.diffs <- tab.f
for(i in 1:nrow(tab.f)){
  name <- row.names(tab.f)[i]
  seq1 <- genes["phix_Fwt"]
  seq2 <- genes[name]
  cods1 <- substring(seq1, seq(1,nchar(seq1),3), seq(3,nchar(seq1),3))
  cods2 <- substring(seq2, seq(1,nchar(seq2),3), seq(3,nchar(seq2),3))
  tab.f[i,] <- cods1==cods2
  tab.f.diffs[i,] <- paste(cods1, cods2, sep = "->")
}
tab.f <- tab.f*1
tab.f$id <- row.names(tab.f)
#write.csv(tab.f, file="results/F_new_diffs.csv", quote = F, row.names = F)


tab.h <- as.data.frame(matrix(ncol=329, nrow=length(H_new)))
row.names(tab.h) <- names(H_new)
names(tab.h) <- as.character(1:329)
tab.h.diffs <- tab.h
for(i in 1:nrow(tab.h)){
  name <- row.names(tab.h)[i]
  seq1 <- genes["phix_Hwt"]
  seq2 <- genes[name]
  cods1 <- substring(seq1, seq(1,nchar(seq1),3), seq(3,nchar(seq1),3))
  cods2 <- substring(seq2, seq(1,nchar(seq2),3), seq(3,nchar(seq2),3))
  tab.h[i,] <- cods1==cods2
  tab.h.diffs[i,] <- paste(cods1, cods2, sep = "->")
}
tab.h <- tab.h*1
tab.h$id <- row.names(tab.h)
#write.csv(tab.h, file="results/H_new_diffs.csv", quote = F, row.names = F)


tab.g <- read.csv("results/G_new_diffs.csv", header = T, stringsAsFactors = F)
names(tab.g) <- c(1:176, "id")
tab.g.l <- melt(tab.g, id.vars = c("id"))
##read.csv.....
names(tab.g.diffs) <- as.character(1:176)
tab.g.diffs$id <- row.names(tab.g.diffs)
tab.g.diffs.l <- melt(tab.g.diffs, id.vars = c("id"))
tab.g.diffs.l <- tab.g.diffs.l[which(tab.g.l$value==0),]
tab.g.diffs.l$mut <- paste(tab.g.diffs.l$variable,tab.g.diffs.l$value,sep=":")
g.mut.table <- table(tab.g.diffs.l$mut)
g.mut.table[g.mut.table==1]
length(g.mut.table[g.mut.table==1])
length(g.mut.table)

g.low.1 <- tab.g.diffs.l[tab.g.diffs.l$id=="G_stability_low_rep1_seq2272",]
others <- tab.g.diffs.l[!tab.g.diffs.l$id=="G_stability_low_rep4_seq2272",]
table(g.low.1$mut %in% others$mut)
