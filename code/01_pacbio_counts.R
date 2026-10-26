# ============================================================================
# 01_pacbio_counts.R
# PacBio genotype counts: quality-filter, trim, and count reads per construct.
#
# Data dependencies (relative to repo root):
#   * data/counts/          (18 PacBio count files)
#
# NOTE: as-written pipeline script; not guaranteed turnkey. Paths are
# relative to the code/ directory.
# ============================================================================

library(stringr)
library(Biostrings)
library(Rsamtools)

#quality filter and trim
#fastp -i INFILE.fq.gz -o OUTFILE.fq.gz --length_required 3000 --length_limit 3500

#########################run_fastp.sh#######################################

#for f in *.fastq.gz; do
  #name=$(echo $f | sed "s/.fastq.gz//g")
#  fread=$(echo -e "${name}.fastq.gz")
#  echo "Trimming $fread"
#  fout=$(echo -e "${name}_filt.fastq.gz")
#  report=$(echo -e "${name}_fastp.html")
#  cmdstring=$(echo "~/bin/fastp -i $fread -o $fout -h $report --thread 24 --length_required 3000 --length_limit 3500")
#  eval $cmdstring
#done


#########################run_bwa.sh#######################################

# for f in *.fastq.gz; do
# name=$(echo $f | sed "s/_filt.fastq.gz//g")
# fread=$(echo -e "${name}_filt.fastq.gz")
# sam=$(echo -e "${name}_bwa.sam")
# bam=$(echo -e "${name}_bwa_mapped.bam")
# sorted1=$(echo -e "${name}_bwa_mapped.sorted")
# sorted2=$(echo -e "${name}_bwa_mapped.sorted.bam")
# sorted3=$(echo -e "${name}_q20.bam")
# out=$(echo -e "${name}_bwa_mapped.mpileup")
# out2=$(echo -e "${name}_bwa_mapped.parsepile.txt")
# out3=$(echo -e "${name}_flagstat.txt")
# 
# cmdstring=$(echo "bwa mem -t 20 refs/all.fasta $fread > $sam")       
# eval $cmdstring
# cmdstring2=$(echo "samtools view -bS -F 4 $sam -o $bam")
# eval $cmdstring2
# cmdstring3=$(echo "samtools sort $bam -o $sorted2")
# eval $cmdstring3
# cmdstring4=$(echo "samtools view -b -q 20 $sorted2 -o $sorted3")
# eval $cmdstring4
# cmdstring5=$(echo "samtools depth $sorted3 > $out2")
# eval $cmdstring5
# cmdstring6=$(echo "samtools flagstat $sorted3 > $out3")


#########################run_seqtk.sh#######################################
# for f in *.bam; do
# name=$(echo $f | sed "s/.bam//g")
# bam=$(echo -e "${name}.bam")
# bam2=$(echo -e "${name}.20mis.bam")  #set to drop reads with more than 20 mismatches. This captures the majority of reads. Allows for up to 0.65% error rate
# report=$(echo -e "${name}_counts.txt")
# 
# cmdstring=$(echo "java -jar  ~/Software/jvarkit/dist/samjdk.jar -e 'return record.getIntegerAttribute("NM").intValue()<20;' $bam | samtools view -S -b - > $bam2")
# eval $cmdstring
# cmdstring2=$(echo "samtools index $bam2")
# eval $cmdstring2
# cmdstring3=$(echo "seqkit bam -C $bam &> $report")
# eval $cmdstring3
# done


bams <- list.files('data/counts')
counts <- bams[str_detect(bams, 'counts')]
names <- str_split(counts, '_q40', simplify = T)[,1]
cntTab <- data.frame(matrix(nrow=127, ncol=18))
names(cntTab) <- names

refs <- readDNAStringSet('data/fastas/all.fasta')
row.names(cntTab) <- names(refs)

for(i in counts){
  cnt <- read.csv(paste("data/counts/", i, sep=""), sep="\t")
  match <- str_split(i, '_q40', simplify = T)[,1]
  cntTab[cnt$Ref,match] <- cnt$Count
}
write.csv(cntTab, 'results/pacbio_may_counts.csv', quote = F, row.names = T)

freqTab <- cntTab
for(i in 1:ncol(freqTab)){
  sum <- sum(freqTab[,i], na.rm = T)
  freqTab[,i] <- freqTab[,i]/sum
}
write.csv(freqTab, 'results/pacbio_may_freqs.csv', quote = F, row.names = T)

freqL <- freqTab
freqL$geno <- row.names(freqL)

piles <- list.files('data/mpileup/')
varTab <- data.frame(matrix(nrow=2, ncol=18))
names(varTab) <- str_split(names, '_q20', simplify = T)[,1]
row.names(varTab) <- c( "anc", "NC_001422")
for(i in piles){
  pileup <- read.csv(paste("data/mpileup/",i, sep=""), sep = "\t")
  pileupSanger <- pileup[pileup$chr=="NC_001422",]
  pileupSanger2811 <- pileupSanger[pileupSanger$loc==2811,c(5,6)]
  samp <- str_split(i, '_parse', simplify = T)[,1]
  varTab[,samp] <- as.numeric(pileupSanger2811)
}

###redo tables
cntTab <- cntTab[-1,]
names(cntTab) <- str_split(names, '_q20', simplify = T)[,1]
cntTab2 <- rbind(varTab, cntTab)
write.csv(cntTab2, 'results/pacbio_may_counts.csv', quote = F, row.names = T)
freqTab <- cntTab2
for(i in 1:ncol(freqTab)){
  sum <- sum(freqTab[,i], na.rm = T)
  freqTab[,i] <- freqTab[,i]/sum
}
write.csv(freqTab, 'results/pacbio_may_freqs.csv', quote = F, row.names = T)





############################full length genome sequencing##################################

##amplicon spans begining and end of reference. Unfortunately, must reindex reference otherwise alignments are not working well.

refs <- readDNAStringSet('data/fastas/all_combos_singles.fasta')

