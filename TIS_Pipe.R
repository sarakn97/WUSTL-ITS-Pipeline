# TIS Pipeline

#initializations
library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")

# List demultiplexed sample fastq files
samples <- list.files(path = "/home/sarakn97/WUSTL/WangL")
# Sort Forward and Reverse Reads
R1s <- sort(list.files(path = "/home/sarakn97/WUSTL/WangL", pattern = "_R1_001.fastq.gz", full.names = TRUE))
R2s <- sort(list.files(path = "/home/sarakn97/WUSTL/WangL", pattern = "_R2_001.fastq.gz", full.names = TRUE))

# Delete files with less than 13 Reads 
for (x in R1s) {
  if ((countLines(x)/4) < 13) {
   unlink(x)
  }
}

for (x in R2s) {
  if ((countLines(x)/4) < 13) {
    unlink(x)
  }
}

# Remove Samples A1, A3, A4, A9, A11, A12, A14
exists <- file.exists(R1s) & file.exists(R2s)
R1s <- R1s[exists]
R2s <- R2s[exists]

# Identify Primers
FWD <- "GCATCGATGAAGAACGCAGC" 
REV <- "TCCTCCGCTTATTGATATGC" 

# Create all orientations of primers
allOrients <- function(primer) {
  require(Biostrings)
  dna <- DNAString(primer)
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), RevComp = reverseComplement(dna))
  return(sapply(orients, toString))
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

# filter sequences, remove those with "N"s
R1s.filtN <- file.path(path = "/home/sarakn97/WUSTL/WangL", "filtN", basename(R1s))
R2s.filtN <- file.path(path = "/home/sarakn97/WUSTL/WangL", "filtN", basename(R2s))
filterAndTrim(R1s, R1s.filtN, R2s, R2s.filtN, maxN = 0, multithread = TRUE)

# Count primer appearances in forward&reverse reads, considering all orientations
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits>0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = R1s.filtN), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = R2s.filtN), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = R1s.filtN), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = R2s.filtN))

# Remove Primers using cutadapt
cutadapt <- "/usr/bin/cutadapt"
# Check to see if cutadapt has been successfully identified
system2(cutadapt, args = "--version")

# Create Output files for cutadapt-ed files & define cutadapt parameters
path <- "/home/sarakn97/WUSTL/WangL"
path.cut <- file.path(path, "cutadapt")
if (!dir.exists(path.cut)) dir.create(path.cut)
R1s.cut <- file.path(path.cut, basename(R1s))
R2s.cut <- file.path(path.cut, basename(R2s))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD & RC of REV off of R1s
R1.flags <- paste("-g", FWD, "-a", REV.RC)
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt, "-n 2" required to remove both FWD & REV from reads
for(i in seq_along(R1s)) {
  system2(cutadapt, args=c(R1.flags, R2.flags, "-n", 2,
                           "-o", R1s.cut[i], "-p", R2s.cut[i], R1s.filtN[i], R2s.filtN[i]))
}

# 3 remaining REV primers of Forward Reads in RC orientation ??
# ^^ Removed C5 to fix issue.... How do I identify the sample with remaining primers after trimming?
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = R1s.cut), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = R2s.cut), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = R1s.cut), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = R2s.cut))

# Sort lists of forward and reverse cutadapt-ed fastq files
cutR1s <- sort(list.files(path.cut, pattern = "_R1_001.fastq.gz", full.names = TRUE))
cutR2s <- sort(list.files(path.cut, pattern = "_R2_001.fastq.gz", full.names = TRUE))


# Filter And Trim
postcut.path1 <- "/home/sarakn97/WUSTL/WangL/cutadapt/FcutR1s"
FiltR1s <- file.path(postcut.path1, "filtered", basename(R1s.Fcut))
postcut.path2 <- "/home/sarakn97/WUSTL/WangL/cutadapt/FcutR2s"
FiltR2s <- file.path(postcut.path2, "filtered", basename(R2s.Fcut))

#maxEE sets max # of expected errors allowed in read
out <- filterAndTrim(cutR1s, FiltR1s, cutR2s, FiltR2s, maxN = 0, maxEE = c(2,2), truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)

# Inspect Read Quality profiles
# Would not work before FilterAndTrime bc some reads in samples had 0 length
plotQualityProfile(FiltR1s[10:11])
plotQualityProfile(FiltR2s[10:11])

# Learn the Error Rates
errR1 <- learnErrors(FiltR1s, multithread = TRUE)
errR2 <- learnErrors(FiltR1s, multithread = TRUE)

# Visualize estimated error rates, all looks good
plotErrors(errR1, nominalQ = TRUE)

# Dereplicate identical reads
derepR1s <- derepFastq(FiltR1s, verbose = TRUE)
derepR2s <- derepFastq(FiltR2s, verbose = TRUE)

# Extract Sample Names
get.sample.name <- function(fname) (basename(fname))
remove.ext <- function(fname) (tools::file_path_sans_ext(fname, compression=TRUE))
sample.names <- unname(sapply(FiltR1s, get.sample.name))
sample.names <- unname(sapply(sample.names, remove.ext))
sample.names <- data.frame("data" = sample.names)
sample.names <- substr(sample.names$data,1, nchar(sample.names$data)-7)
sample.names <- as.character(sample.names)

# Name the derep-class objects by sample names
names(derepR1s) <- sample.names
names(derepR2s) <- sample.names

# Sample Inference, core sample inference algorithm
dadaR1s <- dada(derepR1s, err = errR1, multithread = TRUE)
dadaR2s <- dada(derepR2s, err = errR2, multithread = TRUE)

# Merge Paired Reads
mergers <- mergePairs(dadaR1s, derepR1s, dadaR2s, derepR2s, verbose=TRUE)

# Construct Sequence Table (ASV)
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Remove Chimeras
# Identified 191 bimeras out of 734 input sequences.
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
# inspect distribustion of seq lenghts
table(nchar(getSequences(seqtab.nochim)))

# Track Reads through Pipeline
getN <- function(x) sum(getUniques(x))
trak <- cbind(out, sapply(dadaR1s, getN), sapply(dadaR2s, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(trak) <- c("input", "filtered", "denoisedR1", "denoisedR2", "merged", "nochim")
rownames(trak) <- sample.names
head(trak)
