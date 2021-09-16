## Preprocessing NGS data (quality filter, adaptor removali and mapping) of a barcoded strand-specific RNA-seq via a direct RT step introducing barcode (after rRNA removal)
## The parameter setting here (demultiplexing step) is suitable for 7 nt barcode
## Global alignment is performed for single end read (R2)
## Written by Tianmin Wang
## Jintao Liu Lab, Tsinghua University
## Last updated Aug 10, 2021

# Accept arguments and create directory for log files ////////////////////////
sample=$1
mapping_version=$2
mkdir ${sample}_logs/

# Generate required scripts automatically ////////////////////////////////////
python generate_scripts.py adaptor_R1.fasta ${mapping_version} SE

# FastQC /////////////////////////////////////////////////////////////////////
mkdir fastqcLog
ls "${sample}"_*.fastq.gz | time parallel --bar --results fastqcLog -j8 fastqc {}

# Quality trimming ///////////////////////////////////////////////////////////
export PATH=$PATH:~/.local/bin/
cutadapt -q 20,20 --minimum-length 100:100 --max-n 3 --pair-filter=any -o "${sample}"_QF_R1.fastq.gz -p "${sample}"_QF_R2.fastq.gz "${sample}"_R1.fastq.gz "${sample}"_R2.fastq.gz > "${sample}"_logs/QF.log

# Demultiplexing /////////////////////////////////////////////////////////////
# when the library quality is high, say, each adaptor can be found strictly anchored at the 5' of R1
cutadapt --pair-filter=any --no-indels --minimum-length 20 --times 1 -e 0.15 --overlap 7 -g file:adaptor_R1.fasta -A file:adaptor_R2.fasta -o end5-{name}.R1.fastq.gz -p end5-{name}.R2.fastq.gz "${sample}"_QF_R1.fastq.gz "${sample}"_QF_R2.fastq.gz > "${sample}"_logs/demultiplexing.log

# Mapping to genome //////////////////////////////////////////////////////////
# get genome from NCBI RefSeq: NC000913.3: https://www.ncbi.nlm.nih.gov/nuccore/NC_000913.3 on Sep. 18, 2019
export PATH=$PATH:~/Documents/wtm/bowtie2-2.3.5.1/
export PATH=$PATH:/usr/local/bin/bin/
bash mapping_R2_"${mapping_version}".sh ${sample}

# Summarize the result of demultiplexing and mapping /////////////////////////
python result_summary.py "${sample}"_logs/

# Count the features from alignments /////////////////////////////////////////
Rscript feature_counting.R ./ SE > "${sample}"_logs/featureCounting.log

# Compress intermediate .bam files to save disk space ////////////////////////
bash compress_bam.sh

# Relocate intermediate files ////////////////////////////////////////////////
mv "${sample}"_*fastqc* fastqcLog/
mkdir QF_files
mv "${sample}"_QF*.fastq.gz QF_files/
mkdir demultiplexed_fastq/
mv end5*.fastq.gz demultiplexed_fastq/
mkdir bam_files/
mv lib*.bam.gz sample.csv bam_files/
mkdir rawcounts_tpm/
mv TPM*.csv rawcount.csv rawcounts_tpm/
mv "${sample}"_summary.txt QC_analysis/
mkdir generated_scripts/
mv mapping_R2*.sh compress_bam.sh generated_scripts/
rm -rf QF_files/
