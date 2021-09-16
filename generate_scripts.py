# this script is used to generate bash scripts for mapping (bowtie2) multiple demultiplexed .fastq file
# this script also generates another flat file "sample.csv", fed to the R scripts that perform feature counting
# finally this script generates a bash scipt to compress the intermediate .bam file to save disk space

import os
import sys

# arguments processing
# ///////////////////////////////////////////////////////////////////////////////////////////
barcodes_file = sys.argv[1]
mapping_version = sys.argv[2]
mapping_type = sys.argv[3]
if mapping_version not in ['v1', 'v2']:
    print('incorrect mapping version, v1 or v2 is expected!')
    sys.exit(-1)

if mapping_type not in ['SE', 'PE']:
    print('incorrect mapping type, SE or PE is expected!')
    sys.exit(-1)

# find all barcodes
barcodes = []
f = open(barcodes_file, 'r')
for line in f:
    if line[0] == '>':
        barcode = line.rstrip()[1:]
        if barcode[:7] == 'adaptor':
            barcodes.append(barcode[7:])
        else:
            print('incorrect barcode name for %s in %s'%(barcode, barcodes_file))
            sys.exit(-1)
f.close()

# generate the bash script used in mapping (bowtie2)
# ///////////////////////////////////////////////////////////////////////////////////////////
if mapping_type == 'SE':
    os.system('cat /dev/null > mapping_R2_%s.sh'%(mapping_version))
    g = open('mapping_R2_%s.sh'%(mapping_version), 'r+')
else:
    os.system('cat /dev/null > mapping_%s.sh'%(mapping_version))
    g = open('mapping_%s.sh'%(mapping_version), 'r+')

g.write('sample=$1\n')
g.write('mkdir "${sample}"_logs/bowtie2_logs/\n')
g.write('\n')

def mapping_command(mapping_version, mapping_type, barcode):
    # function to generate mapping command for each line
    this_command = ''
    if mapping_version == 'v1':
        if mapping_type == 'SE':
            this_command = this_command + '(bowtie2 -x MG1655_bowtie2Index/NC000913.3 -U end5-adaptor%s.R2.fastq.gz -S lib%s.sam)2>"${sample}"_logs/bowtie2_logs/adaptor%s_log.txt\n'%(barcode, barcode, barcode)            
        else: #PE
            this_command = this_command + '(bowtie2 -X 1000 -I 18 --no-mixed --no-discordant -x MG1655_bowtie2Index/NC000913.3 -1 end5-adaptor%s.R1.fastq.gz -2 end5-adaptor%s.R2.fastq.gz -S lib%s.sam)2>"${sample}"_logs/bowtie2_logs/adaptor%s_log.txt\n'%(barcode, barcode, barcode, barcode)
    elif mapping_version == 'v2':
        if mapping_type == 'SE':
            this_command = this_command + '(bowtie2 --very-sensitive-local -N 1 -L 16 --no-1mm-upfront --score-min G,9,8 -p 1 -x MG1655_bowtie2Index/NC000913.3 -U end5-adaptor%s.R2.fastq.gz -S lib%s.sam)2>"${sample}"_logs/bowtie2_logs/adaptor%s_log.txt\n'%(barcode, barcode, barcode)
        else: # PE
            this_command = this_command + '(bowtie2 --very-sensitive-local -N 1 -L 16 --no-1mm-upfront --score-min G,9,8 -p 1 -X 1000 -I 18 --no-mixed --no-discordant -x MG1655_bowtie2Index/NC000913.3 -1 end5-adaptor%s.R1.fastq.gz -2 end5-adaptor%s.R2.fastq.gz -S lib%s.sam)2>"${sample}"_logs/bowtie2_logs/adaptor%s_log.txt\n'%(barcode, barcode, barcode, barcode)
    this_command = this_command + 'samtools view -S -b lib%s.sam > lib%s.bam\n'%(barcode, barcode)
    this_command = this_command + 'rm -rf lib%s.sam\n'%(barcode)
    return this_command

for barcode in barcodes:
    this_command = mapping_command(mapping_version, mapping_type, barcode)
    g.write(this_command + '\n')
g.close()

# generate .csv file used in feature counting
# ///////////////////////////////////////////////////////////////////////////////////////////
os.system('cat /dev/null > sample.csv')
g = open('sample.csv', 'r+')
g.write('sample\n')
for barcode in barcodes: 
    g.write('lib%s\n'%(barcode))
g.close()

# generate bash script to compress .bam files
# ///////////////////////////////////////////////////////////////////////////////////////////
os.system('cat /dev/null > compress_bam.sh')
g = open('compress_bam.sh', 'r+')
for barcode in barcodes: 
    g.write('gzip lib%s.bam\n'%(barcode))
g.close()
