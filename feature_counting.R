args <- commandArgs(trailingOnly = TRUE)
ourdir <- args[1]
mode <- args[2]
X <- switch(mode, PE = FALSE, SE = TRUE)
setwd(ourdir)
list.files(ourdir)
csvfile <- file.path(ourdir, "sample.csv")
sampleTable <- read.csv(csvfile)
sampleTable


# import function
Bamfilenames <- file.path(ourdir, paste0(sampleTable$sample, ".bam"))
file.exists(Bamfilenames)

library("Rsamtools")
Bamfiles <- BamFileList(Bamfilenames)
seqinfo(Bamfiles)

library("GenomicFeatures")
gff3file <- file.path(ourdir, "genomeFile/NC000913.3.gff3")
txdb_gff <- makeTxDbFromGFF(gff3file, format = "gff3")

transc <- transcriptsBy(txdb_gff, by = "gene")

# Read counting step
library("GenomicAlignments")
library("BiocParallel")
if(X)
{
se_ecoli <- summarizeOverlaps(feature = transc, reads = Bamfiles,                       
                              mode="Union",
                              ignore.strand=FALSE,
                              singleEnd=X)
} else
{
se_ecoli <- summarizeOverlaps(feature = transc, reads = Bamfiles,                       
                              mode="Union",
                              ignore.strand=FALSE,
                              singleEnd=X,
                              preprocess.reads=invertStrand)
}


head(assay(se_ecoli))

write.csv(assay(se_ecoli),"rawcount.csv", row.names = TRUE)
