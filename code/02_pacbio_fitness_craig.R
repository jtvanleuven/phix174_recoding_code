# ============================================================================
# 02_pacbio_fitness_craig.R
# Craig's PacBio fitness calculation: derive per-construct fitness from
# PacBio read counts (relative to WT).
#
# Data dependencies (relative to repo root):
#   * results/pacbio_may_counts.csv   (from 01)
#   * data/craig/                     (Craig's modeling inputs)
#
# NOTE: as-written pipeline script; not guaranteed turnkey. Paths are
# relative to the code/ directory.
# ============================================================================

setwd("../results/")   # NOTE: run from code/ OR adjust this to your layout
counts <- read.csv(file= "pacbio_may_counts.csv")
rownames(counts) <- counts$X
counts <- counts[,-1]
counts[c(1,2),] <- counts[c(2,1),]
rownames(counts)[c(1,2)] <- rownames(counts)[c(2,1)]
freqs <- read.csv(file="pacbio_may_freqs.csv")
rownames(freqs) <- freqs$X
freqs <- freqs[,-1]
freqs[c(1,2),] <- freqs[c(2,1),]
rownames(freqs)[c(1,2)] <- rownames(freqs)[c(2,1)]
split <- sapply(colnames(freqs), function(x) strsplit(x, "_"))
n.samps <- length(split)
no.data.dex <-which(apply(counts, 1, function(x) length(which(is.na(x))))==n.samps)
no.data.names <- rownames(counts)[no.data.dex]
counts <- counts[-no.data.dex,]
no.data.dex <-which(apply(freqs, 1, function(x) length(which(is.na(x))))==n.samps)
freqs <- freqs[-no.data.dex,]
EV.col.names<- c("EV_r1_mix", "EV_r1_T40", "EV_r1_T80", "EV_r2_mix", "EV_r2_T40",
                  "EV_r2_T80", "EV_r3_mix", "EV_r3_T40", "EV_r3_T80")
EV.col.dex <- sapply(EV.col.names, function(x) which(colnames(counts)==x))
EN.col.names<- c("EN_r1_mix", "EN_r1_T40", "EN_r1_T80", "EN_r2_mix", "EN_r2_T40",
                 "EN_r2_T80", "EN_r3_mix", "EN_r3_T40", "EN_r3_T80")
EN.col.dex <- sapply(EN.col.names, function(x) which(colnames(counts)==x))

EV.nas <- apply(counts[,EV.col.dex], 1, function(x) length(which(is.na(x))))
EN.nas <- apply(counts[,EN.col.dex], 1, function(x) length(which(is.na(x))))
Na.by.var <- as.data.frame(cbind(EN=EN.nas, EV=EV.nas))

n.vars <- dim(freqs)[1]

f.long<- as.data.frame(matrix(nrow=prod(dim(freqs)), ncol=5))
colnames(f.long) <- c("Variant", "Treat", "Rep", "Time", "Freq")
times <- c(0, 40, 80)
names(times) <- c("mix", "T40", "T80")
times.hr <- times/60
vars.v <- rownames(freqs)
treats.v <- c("EN", "EV")
c.long <- f.long
colnames(c.long) <- c("Variant", "Treat", "Rep", "Time", "Count")


r<- 1
for (i in 1:n.vars){
  Var <- rownames(freqs)[i]
  for (j in 1:n.samps){
    f.long$Variant[r] <- Var
    f.long$Treat[r] <- split[[j]][1]
    f.long$Rep[r] <- unlist(strsplit(split[[j]][2], split=""))[2]
    time.pt <- split[[j]][3]
    f.long$Time[r] <- times[which(names(times)==time.pt)]
    f.long$Freq[r] <- freqs[i,j]
    
    c.long$Variant[r] <- Var
    c.long$Treat[r] <- split[[j]][1]
    c.long$Rep[r] <- unlist(strsplit(split[[j]][2], split=""))[2]
    time.pt <- split[[j]][3]
    c.long$Time[r] <- times[which(names(times)==time.pt)]
    c.long$Count[r] <- counts[i,j]
    r <- r + 1
  }
}

wt.name <- "NC_001422"
t <- 40/60

fit.est <- as.data.frame(matrix(nrow=prod(dim(freqs)),ncol=5))
colnames(fit.est) <- c("Variant", "Treat", "Rep", "Times", "Fit")
fit.reg <- as.data.frame(matrix(nrow=n.vars*2*3, ncol=4))
colnames(fit.reg) <- c("Variant", "Treat", "Rep", "Fit")

r <- 1
r.reg <- 1

for (var.i in 1:n.vars){
  var <- vars.v[var.i]
  print(var.i)
  for (treat.i in 1:2){
    treat <- treats.v[treat.i]
    for (rep.i in 1:3){
      rows <- which(c.long$Variant==var & c.long$Treat==treat & c.long$Rep==rep.i)
      r.wt.0 <- which(c.long$Variant==wt.name & c.long$Treat==treat & 
                        c.long$Rep==rep.i & c.long$Time==0)
      r.wt.40 <- which(c.long$Variant==wt.name & c.long$Treat==treat & 
                        c.long$Rep==rep.i & c.long$Time==40)
      r.wt.80 <- which(c.long$Variant==wt.name & c.long$Treat==treat & 
                        c.long$Rep==rep.i & c.long$Time==80)
      r.0 <- rows[which(c.long$Time[rows]==0)]
      r.40 <- rows[which(c.long$Time[rows]==40)]
      r.80 <- rows[which(c.long$Time[rows]==80)]
      rat.0 <- c.long$Count[r.0]/c.long$Count[r.wt.0]
      rat.40 <- c.long$Count[r.40]/c.long$Count[r.wt.40]
      rat.80 <- c.long$Count[r.80]/c.long$Count[r.wt.80]
      fit.0.40 <- log(rat.40/rat.0,2)/t
      fit.0.80 <- log(rat.80/rat.0,2)/(2*t)
      fit.40.80 <- log(rat.80/rat.40,2)/t
      
      fit.est$Variant[c(r, r+1, r+2)] <- var
      fit.est$Treat[c(r, r+1, r+2)] <- treat
      fit.est$Rep[c(r, r+1, r+2)] <- rep.i
      fit.est$Times[r] <- "0.40"
      fit.est$Fit[r] <- fit.0.40
      fit.est$Times[r+1] <- "0.80"
      fit.est$Fit[r+1] <- fit.0.80
      fit.est$Times[r+2] <- "40.80"
      fit.est$Fit[r+2] <- fit.40.80
      r <- r+3
      
      reg.rat.40 <- log((c.long$Count[r.40]/c.long$Count[r.wt.40])/
                    (c.long$Count[r.0]/c.long$Count[r.wt.0]),2)
      reg.rat.80 <- log((c.long$Count[r.80]/c.long$Count[r.wt.80])/
                          (c.long$Count[r.0]/c.long$Count[r.wt.0]),2)
        
      if (is.na(reg.rat.40)==FALSE & is.na(reg.rat.80)==FALSE){
        slope <- lm(c(reg.rat.40, reg.rat.80) ~ times.hr[2:3] - 1)
        plot(times.hr[2:3], c(reg.rat.40, reg.rat.80), xlim=c(0,1.5), ylim=c(-10,2))
        abline(0,slope$coefficients)
        abline(v=0)
        fit.reg$Variant[r.reg] <- var
        fit.reg$Treat[r.reg] <- treat
        fit.reg$Rep[r.reg] <- rep.i
        fit.reg$Fit[r.reg] <- slope$coefficients
        r.reg <- r.reg+1
      }
    }
  }
}

fit.reg.means <- as.data.frame(matrix(nrow=n.vars, ncol=11))
colnames(fit.reg.means) <- c("Fit.mean", "Mean.EN", "Mean.EV", "SE.EN", "SE.EV","EN.1", "EN.2", "EN.3", "EV.1", "EV.2", "EV.3")
rownames(fit.reg.means) <- vars.v
EV.cols <- which(colnames(fit.reg.means) %in% c("EV.1", "EV.2", "EV.3"))
EN.cols <- which(colnames(fit.reg.means) %in% c("EN.1", "EN.2", "EN.3"))

for (i in 1:n.vars){
  var <- vars.v[i]
  dex1 <- which(fit.reg$Variant==var &
                  fit.reg$Treat=="EN" &
                  fit.reg$Rep==1)
  if (length(dex1)==1){
    fit.reg.means$EN.1[i] <- fit.reg$Fit[dex1]
  }
  dex2 <- which(fit.reg$Variant==var & 
                  fit.reg$Treat=="EN" & 
                  fit.reg$Rep==2)
  if (length(dex2)==1){
    fit.reg.means$EN.2[i] <- fit.reg$Fit[dex2]
  }
  dex3 <- which(fit.reg$Variant==var &
                  fit.reg$Treat=="EN" &
                  fit.reg$Rep==3)
  if (length(dex3)==1){
    fit.reg.means$EN.3[i] <- fit.reg$Fit[dex3]
  }
  dex4 <- which(fit.reg$Variant==var &
                  fit.reg$Treat=="EV" &
                  fit.reg$Rep==1)
  if (length(dex4)==1){
    fit.reg.means$EV.1[i] <- fit.reg$Fit[dex4]
  }
  dex5 <- which(fit.reg$Variant==var &
                  fit.reg$Treat=="EV" &
                  fit.reg$Rep==2)
  if (length(dex5)==1){
    fit.reg.means$EV.2[i] <- fit.reg$Fit[dex5]
  }
  dex6 <- which(fit.reg$Variant==var &
                  fit.reg$Treat=="EV" &
                  fit.reg$Rep==3)
  if (length(dex6)==1){
    fit.reg.means$EV.3[i] <- fit.reg$Fit[dex6]
  }
  
  fit.reg.means$Fit.mean[i] <- mean(as.numeric(fit.reg.means[i,c(EN.cols, EV.cols)]), na.rm=T)
  fit.reg.means$Mean.EN[i] <- mean(as.numeric(fit.reg.means[i,EN.cols]), na.rm=T)
  fit.reg.means$Mean.EV[i] <- mean(as.numeric(fit.reg.means[i, EV.cols]), na.rm=T)
  fit.reg.means$SE.EN[i] <- sd(as.numeric(fit.reg.means[i, EN.cols]), na.rm=T)/sqrt(3-length(which(is.na(fit.reg.means[i, EN.cols]))))
  fit.reg.means$SE.EV[i] <- sd(as.numeric(fit.reg.means[i, EV.cols]), na.rm=T)/sqrt(3-length(which(is.na(fit.reg.means[i, EV.cols]))))
}


bg.EN.EV <- c("yellow", "deepskyblue")

plot(1,1,type="n", xlim=c(0, 130), ylim=c(-20,5), xaxt="n", ylab="Fitness effect", xlab="")
for (i in 1:n.vars){
  var <- vars.v[i]
  for (j in 1:2){
    treat <- treats.v[j]
    for (rep.i in 1:3){
      fits <- fit.est$Fit[which(fit.est$Variant==var & fit.est$Treat==treat & fit.est$Rep==rep.i)]
      points(rep(i, length(fits)), fits, pch=21, bg="white", cex=0.8)
    }
  }
}






setwd("../data/craig/")   # NOTE: run from code/ OR adjust this to your layout
flask.fit <- read.csv(file="fitness of construct long form.csv")
flask.names.lookup <- apply(flask.fit, 1, function(x) paste(c(x[2], x[1], x[4], x[5]), collapse="_"))
flask.fit <- cbind(Lookup=flask.names.lookup, flask.fit)

Pacbio.names <- read.csv(file="Pacbio var name lookup.csv")
to.remove <- sapply(no.data.names, function(x) which(Pacbio.names$Name == x))
Pacbio.names <- Pacbio.names[-to.remove,]
lookup.names <- apply(Pacbio.names[,2:5], 1, function(x) paste(x, collapse="_"))
unqs <- unique(lookup.names)
u.list <- lapply(unqs, function(x) which(lookup.names==x))
names(u.list) <- unqs
Pacbio.names <- cbind(Pacbio.names, Lookup=lookup.names, Flask.fit=rep(NA, length(lookup.names)))
for (i in 1:length(Pacbio.names$Lookup)){
  PB.name <- Pacbio.names$Lookup[i]
  index <- which(flask.fit$Lookup == PB.name)
  if (length(index)==1){
    Pacbio.names$Flask.fit[i] <- flask.fit$Fitness[index]
  }
}
wt.index <- which(Pacbio.names$Lookup == "NC_001422_NC_001422_NC_001422_NC_001422")
wt.fit <- Pacbio.names$Flask.fit[wt.index]
Pacbio.names <- cbind(Pacbio.names, Fit.diff=Pacbio.names$Flask.fit-wt.fit)
Pacbio.names <- cbind(Pacbio.names, Fit.Pacbio.reg.mean=rep(NA, length(Pacbio.names$Name)))
for (i in 1:length(Pacbio.names$Name)){
  var <- Pacbio.names$Name[i]
  fit <- fit.reg.means$Fit.mean[which(rownames(fit.reg.means)==var)]
  Pacbio.names$Fit.Pacbio.reg.mean[i] <- fit
}

fit.reg.conf <- apply(fit.reg.means[,2:7], 1, function(x) length(which(is.na(x))))
Pacbio.names <- cbind(Pacbio.names, na.fit.reg=fit.reg.conf)
point.col <- c("white", "grey20", "grey40", "grey60", "grey80", "black", "black")
plot(Pacbio.names$Fit.Pacbio.reg.mean, Pacbio.names$Fit.diff, pch=21, bg=point.col[Pacbio.names$na.fit.reg+1],
     xlab="Pacbio Fitness based on 0-40-80 regression (all flasks)", ylab="Flask based fitness",
     xlim=c(-15,2), ylim=c(-15,2))
lm <- lm(Pacbio.names$Fit.diff ~Pacbio.names$Fit.Pacbio.reg.mean)
abline(lm$coefficients)
abline(0,1, lty="dashed")
abline(0,0, col="grey")
abline(v=0, col="grey")
mtext(side=3, text=paste("R-squared =", round(summary(lm)$r.squared,2)))

rank.table <- as.data.frame(matrix(nrow=n.vars, ncol=9))
colnames(rank.table) <- c("flask", "EN.Reg", "EV.Reg", "EN.0.40", "EV.0.40", "EN.0.80", "EV.0.80",
                          "EN.40.80", "EV.40.80")
rownames(rank.table) <- vars.v
var.dex <- which(is.na(Pacbio.names$Fit.diff)==FALSE)
var.dex <- var.dex[-1]
rank.table$flask[var.dex] <- rank(Pacbio.names$Fit.diff[var.dex], na.last="keep")
rank.table$EN.Reg[var.dex] <- rank(fit.reg.means$Mean.EN[var.dex], na.last="keep")
rank.table$EV.Reg[var.dex] <- rank(fit.reg.means$Mean.EV[var.dex], na.last="keep")




# Do we need the 80 minute sample?
fit.0.40.df <- as.data.frame(matrix(nrow=n.vars, ncol=11))
colnames(fit.0.40.df) <- c("All.mean", "EN.mean", "EV.mean", 
                           "EN.1", "EN.2", "EN.3", "EV.1", "EV.2", "EV.3",
                           "SE.EN", "SE.EV")
rownames(fit.0.40.df) <- vars.v
fit.40.80.df <- fit.0.40.df
fit.0.80.df <- fit.0.40.df

for (i in 1:n.vars){
  var <- vars.v[i]
  dex.EN1 <- which(fit.est$Variant==var & fit.est$Treat=="EN" & fit.est$Rep==1)
  dex.EN2 <- which(fit.est$Variant==var & fit.est$Treat=="EN" & fit.est$Rep==2)
  dex.EN3 <- which(fit.est$Variant==var & fit.est$Treat=="EN" & fit.est$Rep==3)
  dex.EV1 <- which(fit.est$Variant==var & fit.est$Treat=="EV" & fit.est$Rep==1)
  dex.EV2 <- which(fit.est$Variant==var & fit.est$Treat=="EV" & fit.est$Rep==2)
  dex.EV3 <- which(fit.est$Variant==var & fit.est$Treat=="EV" & fit.est$Rep==3)
  fit.0.40.df$EN.1[i] <- fit.est$Fit[dex.EN1[which(fit.est$Time[dex.EN1]=="0.40")]] 
  fit.0.40.df$EN.2[i] <- fit.est$Fit[dex.EN2[which(fit.est$Time[dex.EN2]=="0.40")]]
  fit.0.40.df$EN.3[i] <- fit.est$Fit[dex.EN3[which(fit.est$Time[dex.EN3]=="0.40")]] 
  fit.0.80.df$EN.1[i] <- fit.est$Fit[dex.EN1[which(fit.est$Time[dex.EN1]=="0.80")]] 
  fit.0.80.df$EN.2[i] <- fit.est$Fit[dex.EN2[which(fit.est$Time[dex.EN2]=="0.80")]]
  fit.0.80.df$EN.3[i] <- fit.est$Fit[dex.EN3[which(fit.est$Time[dex.EN3]=="0.80")]] 
  fit.40.80.df$EN.1[i] <- fit.est$Fit[dex.EN1[which(fit.est$Time[dex.EN1]=="40.80")]] 
  fit.40.80.df$EN.2[i] <- fit.est$Fit[dex.EN2[which(fit.est$Time[dex.EN2]=="40.80")]]
  fit.40.80.df$EN.3[i] <- fit.est$Fit[dex.EN3[which(fit.est$Time[dex.EN3]=="40.80")]]
  
  fit.0.40.df$EV.1[i] <- fit.est$Fit[dex.EV1[which(fit.est$Time[dex.EV1]=="0.40")]] 
  fit.0.40.df$EV.2[i] <- fit.est$Fit[dex.EV2[which(fit.est$Time[dex.EV2]=="0.40")]]
  fit.0.40.df$EV.3[i] <- fit.est$Fit[dex.EV3[which(fit.est$Time[dex.EV3]=="0.40")]] 
  fit.0.80.df$EV.1[i] <- fit.est$Fit[dex.EV1[which(fit.est$Time[dex.EV1]=="0.80")]] 
  fit.0.80.df$EV.2[i] <- fit.est$Fit[dex.EV2[which(fit.est$Time[dex.EV2]=="0.80")]]
  fit.0.80.df$EV.3[i] <- fit.est$Fit[dex.EV3[which(fit.est$Time[dex.EV3]=="0.80")]] 
  fit.40.80.df$EV.1[i] <- fit.est$Fit[dex.EV1[which(fit.est$Time[dex.EV1]=="40.80")]] 
  fit.40.80.df$EV.2[i] <- fit.est$Fit[dex.EV2[which(fit.est$Time[dex.EV2]=="40.80")]]
  fit.40.80.df$EV.3[i] <- fit.est$Fit[dex.EV3[which(fit.est$Time[dex.EV3]=="40.80")]]
  
  fit.0.40.df$All.mean[i] <- mean(as.numeric(fit.0.40.df[i,4:9]), na.rm=T)
  fit.0.40.df$EN.mean[i] <- mean(as.numeric(fit.0.40.df[i,4:6]), na.rm=T)
  fit.0.40.df$EV.mean[i] <- mean(as.numeric(fit.0.40.df[i,7:9]), na.rm=T)
  fit.0.40.df$SE.EN[i] <- sd(as.numeric(fit.0.40.df[i,4:6]))/sqrt(length(which(is.na(fit.0.40.df[i,4:6])==F)))
  fit.0.40.df$SE.EV[i] <- sd(as.numeric(fit.0.40.df[i,7:9]))/sqrt(length(which(is.na(fit.0.40.df[i,7:9])==F)))
  
  fit.0.80.df$All.mean[i] <- mean(as.numeric(fit.0.80.df[i,4:9]), na.rm=T)
  fit.0.80.df$EN.mean[i] <- mean(as.numeric(fit.0.80.df[i,4:6]), na.rm=T)
  fit.0.80.df$EV.mean[i] <- mean(as.numeric(fit.0.80.df[i,7:9]), na.rm=T)
  fit.0.80.df$SE.EN[i] <- sd(as.numeric(fit.0.80.df[i,4:6]))/sqrt(length(which(is.na(fit.0.80.df[i,4:6])==F)))
  fit.0.80.df$SE.EV[i] <- sd(as.numeric(fit.0.80.df[i,7:9]))/sqrt(length(which(is.na(fit.0.80.df[i,7:9])==F)))
  
  fit.40.80.df$All.mean[i] <- mean(as.numeric(fit.40.80.df[i,4:9]), na.rm=T)
  fit.40.80.df$EN.mean[i] <- mean(as.numeric(fit.40.80.df[i,4:6]), na.rm=T)
  fit.40.80.df$EV.mean[i] <- mean(as.numeric(fit.40.80.df[i,7:9]), na.rm=T)
  fit.40.80.df$SE.EN[i] <- sd(as.numeric(fit.40.80.df[i,4:6]))/sqrt(length(which(is.na(fit.40.80.df[i,4:6])==F)))
  fit.40.80.df$SE.EV[i] <- sd(as.numeric(fit.40.80.df[i,7:9]))/sqrt(length(which(is.na(fit.40.80.df[i,7:9])==F)))
}

# Which treatment (EN or EV) and which comparision (0-40, 0-80, 40-80)
# has the most consistent results--as quantified by lowest SE?
SE.comp <- as.data.frame(matrix(nrow=4, ncol=4))
colnames(SE.comp) <- c("Pairwise.EN", "Pairwise.EV", "Regress.EN", "Regress.EV")
rownames(SE.comp) <- c("0.40", "40.80", "0.80", "All")
n.var.est <- SE.comp

SE.comp$Pairwise.EN[1] <- round(mean(fit.0.40.df$SE.EN[2:n.vars], na.rm=T),2)
SE.comp$Pairwise.EV[1] <- round(mean(fit.0.40.df$SE.EV[2:n.vars], na.rm=T),2)
SE.comp$Regress.EN[4] <- round(mean(fit.reg.means$SE.EN[2:n.vars], na.rm=T),2)
SE.comp$Regress.EV[4] <- round(mean(fit.reg.means$SE.EV[2:n.vars], na.rm=T),2)
SE.comp$Pairwise.EN[2] <- round(mean(fit.40.80.df$SE.EN[2:n.vars], na.rm=T),2)
SE.comp$Pairwise.EV[2] <- round(mean(fit.40.80.df$SE.EV[2:n.vars], na.rm=T),2)
SE.comp$Pairwise.EN[3] <- round(mean(fit.0.80.df$SE.EN[2:n.vars], na.rm=T),2)
SE.comp$Pairwise.EV[3] <- round(mean(fit.0.80.df$SE.EV[2:n.vars], na.rm=T),2)

tech.reps <- read.csv(file="fit_summary_techreps_CRM.csv")
Names <- apply(tech.reps, 1, function(x) paste(c(x[1], x[2], x[3], x[4]), collapse="_"))
wt.dex <- which(Names=="NC_001422_NC_001422_NC_001422_NC_001422")
Effect <- tech.reps$Mean[wt.dex] - tech.reps$Mean
tech.reps <- cbind(Name=Names, tech.reps, Effect=Effect)
tech.reps.cons <- tech.reps[-wt.dex,]


SE.table <- as.data.frame(matrix(nrow=9, ncol=1))
colnames(SE.table) <- "SE"
rownames(SE.table) <- c("flask", "EN.0.40", "EV.0.40", "EN.40.80", "EV.40.80", "EN.0.80", 
                        "EV.0.80", "EN.Reg", "EV.Reg")
SE.table[1,1] <- round(mean(tech.reps.cons$SE, na.rm=T),2)
SE.table["EN.0.40",1] <- round(mean(fit.0.40.df$SE.EN[2:n.vars], na.rm=T),2)
SE.table["EV.0.40",1] <- round(mean(fit.0.40.df$SE.EV[2:n.vars], na.rm=T),2)
SE.table["EN.0.80",1] <- round(mean(fit.0.80.df$SE.EN[2:n.vars], na.rm=T),2)
SE.table["EV.0.80",1] <- round(mean(fit.0.80.df$SE.EV[2:n.vars], na.rm=T),2)
SE.table["EN.40.80",1] <- round(mean(fit.40.80.df$SE.EN[2:n.vars], na.rm=T),2)
SE.table["EV.40.80",1] <- round(mean(fit.40.80.df$SE.EV[2:n.vars], na.rm=T),2)
SE.table["EN.Reg",1] <- round(mean(fit.reg.means$SE.EN[2:n.vars], na.rm=T),2)
SE.table["EV.Reg",1] <- round(mean(fit.reg.means$SE.EV[2:n.vars], na.rm=T),2)



# === Back to comparing Pacbio results with flask assay results
rank.table$EN.0.40[var.dex] <- rank(fit.0.40.df$EN.mean[var.dex], na.last="keep")
rank.table$EV.0.40[var.dex] <- rank(fit.0.40.df$EV.mean[var.dex], na.last="keep")
rank.table$EN.0.80[var.dex] <- rank(fit.0.80.df$EN.mean[var.dex], na.last="keep")
rank.table$EV.0.80[var.dex] <- rank(fit.0.80.df$EV.mean[var.dex], na.last="keep")
rank.table$EN.40.80[var.dex] <- rank(fit.40.80.df$EN.mean[var.dex], na.last="keep")
rank.table$EV.40.80[var.dex] <- rank(fit.40.80.df$EV.mean[var.dex], na.last="keep")

rank.cor.table <- as.data.frame(matrix(nrow=9, ncol=9))
rownames(rank.cor.table) <- colnames(rank.table)
colnames(rank.cor.table) <- colnames(rank.table)
for (i in 1:9){
  for (j in 1:9){
    rank.cor.table[i,j] <- round(cor.test(rank.table[,i], rank.table[,j], method="spearman")$estimate,2)
  }
}

fit.eff.table <- as.data.frame(cbind(Flask=Pacbio.names$Fit.diff, 
                       EN.Reg=fit.reg.means$Mean.EN,
                       EV.Reg=fit.reg.means$Mean.EV,
                       EN.0.40=fit.0.40.df$EN.mean,
                       EV.0.40=fit.0.40.df$EV.mean,
                       EN.0.80=fit.0.80.df$EN.mean,
                       EV.0.80=fit.0.80.df$EV.mean,
                       EN.40.80=fit.40.80.df$EN.mean,
                       EV.40.80=fit.40.80.df$EV.mean))
rownames(fit.eff.table) <- Pacbio.names$Name

write.csv(fit.eff.table, file="Fitness_effects_table.csv", row.names=T)
fit.cor.table <- as.data.frame(matrix(nrow=9, ncol=9))
rownames(fit.cor.table) <- colnames(fit.eff.table)
colnames(fit.cor.table) <- colnames(fit.eff.table)
for (i in 1:9){
  for (j in 1:9){
    fit.cor.table[i,j] <- round(cor.test(fit.eff.table[,i], fit.eff.table[,j], method="pearson")$estimate,2)
  }
}


n.var.est <- as.data.frame(matrix(nrow=8, ncol=4))
colnames(n.var.est) <- c(3,2,1,0)
rownames(n.var.est) <- c("EN.0.40", "EV.0.40", "EN.40.80", "EV.40.80", "EN.0.80", 
                         "EV.0.80", "EN.Reg", "EV.Reg")

x <- apply(fit.0.40.df[,EN.cols], 1, function(x) length(which(is.na(x)==F)))
n.var.est[1,] <- sapply(seq(3,0), function(y) length(which(x==y)))
x <- apply(fit.0.40.df[,EV.cols], 1, function(x) length(which(is.na(x)==F)))
n.var.est[2,] <- sapply(seq(3,0), function(y) length(which(x==y)))

x <- apply(fit.40.80.df[,EN.cols], 1, function(x) length(which(is.na(x)==F)))
n.var.est[3,] <- sapply(seq(3,0), function(y) length(which(x==y)))
x <- apply(fit.40.80.df[,EV.cols], 1, function(x) length(which(is.na(x)==F)))
n.var.est[4,] <- sapply(seq(3,0), function(y) length(which(x==y)))

x <- apply(fit.0.80.df[,EN.cols], 1, function(x) length(which(is.na(x)==F)))
n.var.est[5,] <- sapply(seq(3,0), function(y) length(which(x==y)))
x <- apply(fit.0.80.df[,EV.cols], 1, function(x) length(which(is.na(x)==F)))
n.var.est[6,] <- sapply(seq(3,0), function(y) length(which(x==y)))

x <- apply(fit.reg.means[,EN.cols], 1, function(x) length(which(is.na(x)==F)))
n.var.est[7,] <- sapply(seq(3,0), function(y) length(which(x==y)))
x <- apply(fit.reg.means[,EV.cols], 1, function(x) length(which(is.na(x)==F)))
n.var.est[8,] <- sapply(seq(3,0), function(y) length(which(x==y)))

est.na.by.var <- as.data.frame(matrix(nrow=n.vars, ncol=8))
colnames(est.na.by.var) <- c("EN.0.40", "EV.0.40", "EN.40.80", "EV.40.80", "EN.0.80", 
                                "EV.0.80", "EN.Reg", "EV.Reg")
rownames(est.na.by.var) <- vars.v
est.na.by.var$EN.0.40 <- apply(fit.0.40.df[,EN.cols], 1, function(x) length(which(is.na(x))))
est.na.by.var$EV.0.40 <- apply(fit.0.40.df[,EV.cols], 1, function(x) length(which(is.na(x))))
est.na.by.var$EN.0.80 <- apply(fit.0.80.df[,EN.cols], 1, function(x) length(which(is.na(x))))
est.na.by.var$EV.0.80 <- apply(fit.0.80.df[,EV.cols], 1, function(x) length(which(is.na(x))))
est.na.by.var$EN.40.80 <- apply(fit.40.80.df[,EN.cols], 1, function(x) length(which(is.na(x))))
est.na.by.var$EV.40.80 <- apply(fit.40.80.df[,EV.cols], 1, function(x) length(which(is.na(x))))
est.na.by.var$EN.Reg <- apply(fit.reg.means[,EN.cols], 1, function(x) length(which(is.na(x))))
est.na.by.var$EV.Reg <- apply(fit.reg.means[,EV.cols], 1, function(x) length(which(is.na(x))))
  



# ===================================================
# ================ For Phage Meeting ================
# ===================================================

# ===== How correlated are Pacbio estimates of fitness to flask estiamtes? =====
fit.cor.table
rank.cor.table


apply(est.na.by.var, 2, function(x) length(which(x==0)))

bg.by.na <- c("white", "grey70", "grey50", "red")
lm <- lm(fit.eff.table$Flask ~ fit.eff.table$EV.Reg)
r <- fit.cor.table$Flask[3]
plot(fit.eff.table$EV.Reg, fit.eff.table$Flask,
     main=paste("R2=", round(summary(lm)$r.squared,2), ", r=", r),
     pch=21, bg=bg.by.na[est.na.by.var$EV.Reg+1],
     xlab="EV Regression", ylab="Flask",
     xlim=c(-15,2), ylim=c(-15,2))
abline(lm$coefficients)
abline(0,1, lty="dashed")
legend("topleft", legend=seq(3,0), bty="n", pch=c(21,21,21,22), pt.bg=bg.by.na, title="# reps observed (of 3)")

r <- rank.cor.table$flask[3]
plot(rank.table$EV.Reg, rank.table$flask,
     main=paste("r=", r),
     pch=21, bg=bg.by.na[est.na.by.var$EV.Reg+1],
     xlab="Rank from EV Regression", ylab="Rank from Flask",
     xlim=c(0,120), ylim=c(0,120))
dropouts <- which(is.na(fit.eff.table$EV.Reg))
points(rep(122, length(dropouts)), fit.eff.table$Flask[dropouts], pch=22, bg="red", cex=1.2)
legend("topleft", legend=seq(3,0), bty="n", pch=c(21,21,21,22), pt.bg=bg.by.na, title="# reps obs")


lm <- lm(fit.eff.table$Flask ~ fit.eff.table$EN.Reg)
r <- fit.cor.table$Flask[2]
plot(fit.eff.table$EN.Reg, fit.eff.table$Flask,
     main=paste("R2=", round(summary(lm)$r.squared,2), ", r=", r),
     pch=21, bg=bg.by.na[est.na.by.var$EN.Reg+1],
     xlab="EN Regression", ylab="Flask",
     xlim=c(-15,2), ylim=c(-15,2))
abline(lm$coefficients)
abline(0,1, lty="dashed")
dropouts <- which(is.na(fit.eff.table$EN.Reg))
points(rep(-15.5, length(dropouts)), fit.eff.table$Flask[dropouts], pch=22, bg="red", cex=1.2)
legend("topleft", legend=seq(3,0), bty="n", pch=c(21,21,21,22), pt.bg=bg.by.na, title="# reps observed (of 3)")

r <- rank.cor.table$flask[2]
plot(rank.table$EN.Reg, rank.table$flask,
     main=paste("r=", r),
     pch=21, bg=bg.by.na[est.na.by.var$EN.Reg+1],
     xlab="Rank from EN Regression", ylab="Rank from Flask",
     xlim=c(0,120), ylim=c(0,120))
dropouts <- which(is.na(fit.eff.table$EN.Reg))
points(rep(122, length(dropouts)), rank.table$flask[dropouts], pch=22, bg="red", cex=1.2)
legend("left", legend=seq(3,0), pch=c(21,21,21,22), pt.bg=bg.by.na, title="# reps obs")


lm <- lm(fit.eff.table$EN.Reg ~ fit.eff.table$EV.0.40)
r <- fit.cor.table$EN.Reg[5]
plot(fit.eff.table$EV.0.40, fit.eff.table$EN.Reg,
     main=paste("R2=", round(summary(lm)$r.squared,2), ", r=", r),
     pch=21, bg=bg.by.na[est.na.by.var$EN.Reg+1],
     xlab="EV 0 to 40", ylab="EN Regression",
     xlim=c(-15,2), ylim=c(-15,2))
abline(lm$coefficients)
abline(0,1, lty="dashed")
dropouts <- which(is.na(fit.eff.table$EN.Reg))
points(rep(-15.5, length(dropouts)), fit.eff.table$Flask[dropouts], pch=22, bg="red", cex=1.2)
legend("topleft", legend=seq(3,0), bty="n", pch=c(21,21,21,22), pt.bg=bg.by.na, title="# reps observed (of 3)")


r <- fit.cor.table$EN.0.40[1]
lm <- lm(fit.eff.table$Flask ~ fit.eff.table$EN.0.40)
plot(fit.eff.table$EN.0.40, fit.eff.table$Flask,
     main=paste("R2=", round(summary(lm)$r.squared,2), ", r=", r),
     pch=21, bg=bg.by.na[Na.by.var$EN+1],
     xlab="EN 0 to 40", ylab="Flask",
     xlim=c(-15,2), ylim=c(-15,2))
abline(lm$coefficients)
abline(0,1, lty="dashed")
dropouts <- which(is.na(fit.eff.table$EN.0.40))
points(rep(-15.5, length(dropouts)), fit.eff.table$Flask[dropouts], pch=22, bg="red", cex=1.2)
legend("topleft", legend=seq(3,0), bty="n", pch=c(21,21,21,22), pt.bg=bg.by.na, title="# reps observed (of 3)")


r <- fit.cor.table$EN.0.80[1]
lm <- lm(fit.eff.table$Flask ~ fit.eff.table$EN.0.80)
plot(fit.eff.table$EN.0.80, fit.eff.table$Flask,
     main=paste("R2=", round(summary(lm)$r.squared,2), ", r=", r),
     pch=21, bg=bg.by.na[Na.by.var$EN+1],
     xlab="EN 0 to 80", ylab="Flask",
     xlim=c(-15,2), ylim=c(-15,2))
abline(lm$coefficients)
abline(0,1, lty="dashed")
dropouts <- which(is.na(fit.eff.table$EN.0.80))
points(rep(-15.5, length(dropouts)), fit.eff.table$Flask[dropouts], pch=22, bg="red", cex=1.2)
legend("topleft", legend=seq(3,0), bty="n", pch=c(21,21,21,22), pt.bg=bg.by.na, title="# reps obs")



r <- fit.cor.table$EV.0.40[1]
lm <- lm(fit.eff.table$Flask ~ fit.eff.table$EV.0.40)
plot(fit.eff.table$EV.0.40, fit.eff.table$Flask,
     main=paste("R2=", round(summary(lm)$r.squared,2), ", r=", r),
     pch=21, bg=bg.by.na[Na.by.var$EV+1],
     xlab="EV 0 to 40", ylab="Flask",
     xlim=c(-15,2), ylim=c(-15,2))
abline(lm$coefficients)
abline(0,1, lty="dashed")
dropouts <- which(is.na(fit.eff.table$EV.0.40))
points(rep(-15.5, length(dropouts)), fit.eff.table$Flask[dropouts], pch=22, bg="red", cex=1.2)
legend("topleft", legend=seq(3,0), bty="n", pch=c(21,21,21,22), pt.bg=bg.by.na, title="# reps obs")


r <- fit.cor.table$EV.0.80[1]
lm <- lm(fit.eff.table$Flask ~ fit.eff.table$EV.0.80)
plot(fit.eff.table$EV.0.80, fit.eff.table$Flask,
     main=paste("R2=", round(summary(lm)$r.squared,2), ", r=", r),
     pch=21, bg=bg.by.na[Na.by.var$EV+1],
     xlab="EV 0 to 80", ylab="Flask",
     xlim=c(-15,2), ylim=c(-15,2))
abline(lm$coefficients)
abline(0,1, lty="dashed")
dropouts <- which(is.na(fit.eff.table$EV.0.80))
points(rep(-15.5, length(dropouts)), fit.eff.table$Flask[dropouts], pch=22, bg="red", cex=1.2)
legend("topleft", legend=seq(3,0), bty="n", pch=c(21,21,21,22), pt.bg=bg.by.na, title="# reps obs")



# ====== How similar are 0 to 40 and 0 to 80 samples?
# Answer: they are highly correlated. 
# 0.80 seems to hit a floor that 0.40 doesn't
lm <- lm(fit.eff.table$EN.0.40 ~ fit.eff.table$EN.0.80)
r <- fit.cor.table$EN.0.40[6]
plot(fit.eff.table$EN.0.40, fit.eff.table$EN.0.80,
     main=paste("R2=", round(summary(lm)$r.squared,2), ", r=", r),
     pch=21, bg=bg.by.na[Na.by.var$EN+1],
     xlab="EN 0 t0.40", ylab="EN 0 to 80",
     xlim=c(-15,2), ylim=c(-15,2))
abline(lm$coefficients)
abline(0,1, lty="dashed")

lm <- lm(fit.eff.table$EV.0.40 ~ fit.eff.table$EV.0.80)
r <- fit.cor.table$EV.0.40[7]
plot(fit.eff.table$EV.0.40, fit.eff.table$EV.0.80,
     main=paste("R2=", round(summary(lm)$r.squared,2), ", r=", r),
     pch=21, bg=bg.by.na[Na.by.var$EN+1],
     xlab="EV 0 t0.40", ylab="EV 0 to 80",
     xlim=c(-15,2), ylim=c(-15,2))
abline(lm$coefficients)
abline(0,1, lty="dashed")


# How correlated are the EV and EN treatmentws for the same time?
lm <- lm(fit.eff.table$EV.0.40 ~ fit.eff.table$EN.0.40)
r <- fit.cor.table$EV.0.40[4]
plot(fit.eff.table$EV.0.40, fit.eff.table$EN.0.40,
     main=paste("R2=", round(summary(lm)$r.squared,2), ", r=", r),
     pch=21, bg=bg.by.na[Na.by.var$EN+1],
     xlab="EV.0.40", ylab="EN.0.40",
     xlim=c(-15,2), ylim=c(-15,2))
abline(lm$coefficients)
abline(0,1, lty="dashed")

lm <- lm(fit.eff.table$EV.0.80 ~ fit.eff.table$EN.0.80)
r <- fit.cor.table$EV.0.80[4]
plot(fit.eff.table$EV.0.80, fit.eff.table$EN.0.80,
     main=paste("R2=", round(summary(lm)$r.squared,2), ", r=", r),
     pch=21, bg=bg.by.na[Na.by.var$EN+1],
     xlab="EV.0.80", ylab="EN.0.80",
     xlim=c(-15,2), ylim=c(-15,2))
abline(lm$coefficients)
abline(0,1, lty="dashed")


# What treatment is the most consistent?
# SE.table shows you that, EN outperforms EV
# and 0.80 outperforms 0.40. Reg similar to 0.80
# Notably, flask is not great with SE=0.37
SE.table


# Is SE related to fitness?
# For flask-based fitness, no:
wt.fit <- Pacbio.names$Flask.fit[1]
plot(tech.reps.cons$Mean-wt.fit, tech.reps.cons$SE,
     xlab="Flask fitness effect", ylab="Standard Error", main="Flask")
abline(0.35, 0, lty="dashed")

# For Pacbio-based fitness, yes:

plot(fit.0.80.df$EV.mean, fit.0.80.df$SE.EV,
     xlab="Pacbio 0 to 80 EV fitness effect", ylab="Standard Error",
     xlim=c(-10,2), main="0.80 EN and EV")
abline(0.35, 0, lty="dashed")
points(fit.0.80.df$EN.mean, fit.0.80.df$SE.EN, pch=21, bg="red")
legend("topright", legend=c("EV", "EN"), pch=21, pt.bg=c("white", "red"), bty="n")
# Adding the red dots shows you the value of even numbers treatment: 
# you get more consistent estimates of the low fitness variations


# Question: is low vs high SE a function of high vs low read count?
# What time point should we use?  Endpoint probably.
# This plot shows that high variance estimates come form low counts

EV.T80.count.dex <- sapply(c("EV_r1_T80", "EV_r2_T80", "EV_r3_T80"), function(x) which(colnames(counts)==x))
mean.EV.T80.counts <- apply(counts[,EV.T80.count.dex], 1, function(x) mean(x, na.rm=T))
plot(log(mean.EV.T80.counts, 10), fit.0.80.df$SE.EV)
abline(0.37, 0, lty="dashed")



# Thoughts/issues: 
# 1. Low fitness variants--are flask estimates biased low and Pacbio more accurate?
# 2. What fitness values should we use in the paper?
# 3. Pacbio can test if evolution has occurred within the growth assay. Need to check this.
# 4. JT will split wildtype--need to incorporate this in analysis









plot(tech.reps.cons$Mean, tech.reps.cons$SE)
plot(fit.0.80.df$EV.mean, fit.0.80.df$SE.EV)
plot(fit.0.80.df$EN.mean, fit.0.80.df$SE.EN)
plot(fit.0.40.df$EV.mean, fit.0.40.df$SE.EV)
plot(fit.reg.means$Mean.EV, fit.reg.means$SE.EV)

plot(fit.40.80.df$EN.mean, Pacbio.names$Fit.diff)
plot(fit.40.80.df$EV.mean, Pacbio.names$Fit.diff)

plot(fit.0.40.df$EN.mean, fit.40.80.df$EN.mean)
plot(fit.0.40.df$EV.mean, fit.40.80.df$EV.mean)

  





#  What level of correlation would we expect if we ran flask assays all over?
tech.reps <- read.csv(file="fit_summary_techreps_CRM.csv")
Names <- apply(tech.reps, 1, function(x) paste(c(x[1], x[2], x[3], x[4]), collapse="_"))
wt.dex <- which(Names=="wt_wt_wt_wt")
Effect <- tech.reps$Mean[wt.dex] - tech.reps$Mean
tech.reps <- cbind(Name=Names, tech.reps, Effect=Effect)
tech.reps.cons <- tech.reps[-wt.dex,]
n <- dim(tech.reps.cons)[1]
boot.set1 <- rnorm(n=n, mean=tech.reps.cons$Mean, sd=tech.reps.cons$SE)
boot.set2 <- rnorm(n=n, mean=tech.reps.cons$Mean, sd=tech.reps.cons$SE)
lm <- lm(boot.set2 ~ boot.set1)
plot(boot.set1, boot.set2,
     main=paste("R2=", round(summary(lm)$r.squared,3)))
abline(lm$coefficients)





