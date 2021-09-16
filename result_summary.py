# this sciprt is used to extract the essential information from all .log files
# it is expected to be executed under the working direcotry of barcoded_RNAseq_directRT where *_logs/ is
# run by python result_summary.py prefix_logs/
# written by Tianmin Wang from Liu Lab @ Tsinghua University
# last update, April 14, 2021

import os
import sys
import numpy as np

log_path = sys.argv[1]

## ////////////////////////////////////////////////
## Processing the log files

# QF and demultiplexing logs
QF_log = '%s/QF.log'%(log_path)
demultiplexing_log = '%s/demultiplexing.log'%(log_path)

# get all mapping logs
mapping_log_path = '%s/bowtie2_logs/'%(log_path)

os.system('ls %s > tmp.txt'%(mapping_log_path))
mapping_logs = []
f = open('tmp.txt', 'r')
for line in f:
    log_file = line.rstrip()
    mapping_logs.append('%s%s'%(mapping_log_path, log_file))
f.close()
os.system('rm -rf tmp.txt')

# get all files before mapping
os.system('ls |egrep "end5" |egrep "R2" |egrep -v "unknown" |egrep "adaptor" > tmp.txt')
dem_files = []
f = open('tmp.txt', 'r')
for line in f:
    dem_file = line.rstrip()
    dem_files.append(dem_file)
f.close()
os.system('rm -rf tmp.txt')

prior_dem_file = '%s_%s'%(log_path.split('_logs')[0], 'QF_R2.fastq.gz')

## ////////////////////////////////////////////////
## Summarize the result

# functions
def extract_QF(QF_log):
    f = open(QF_log, 'r')
    for line in f:
        if 'Total written (filtered):' in line:
            ratio = line.rstrip().split(' bp (')[1].split('%)')[0]
            break        
    f.close()
    return float(ratio)

def extract_demultiplexing(demultiplexing_log):
    f = open(demultiplexing_log, 'r')
    for line in f:
        if 'Read 1 with adapter' in line:
            R1_ratio = line.rstrip().split(' (')[1].split('%)')[0]
        if 'Read 2 with adapter' in line:
            R2_ratio = line.rstrip().split(' (')[1].split('%)')[0]
    f.close()
    return (float(R1_ratio), float(R2_ratio))

def extract_filter_short(dem_files, prior_dem_file):
    read_counts = []
    for dem_file in dem_files:
        os.system('zless %s |wc -l > tmp.txt'%(dem_file))    
        f = open('tmp.txt', 'r')
        for i,line in enumerate(f):
            if i == 0:
                reads = float(line.split(' ')[0])/4
                read_counts.append(reads)
                break
        f.close()
        os.system('rm tmp.txt')
    os.system('zless %s |wc -l > tmp.txt'%(prior_dem_file))
    f = open('tmp.txt', 'r')
    for i,line in enumerate(f):
        if i == 0:
            reads = float(line.split(' ')[0])/4
            prior_read_counts = reads
    f.close()
    os.system('rm tmp.txt')
    return (read_counts, prior_read_counts)

def extract_mapping(mapping_logs):
    mapping_ratios = []
    for mapping_log in mapping_logs:
        f = open(mapping_log, 'r')
        for i,line in enumerate(f):
            if '% overall alignment rate' in line:
                mapping_ratio = float(line.rstrip().split('% overall alignment rate')[0])
                mapping_ratios.append(mapping_ratio)
                break
        f.close()
    return mapping_ratios

# summarize one by one
os.system('cat /dev/null > %s_summary.txt'%(log_path.split('_logs')[0]))
g = open('%s_summary.txt'%(log_path.split('_logs')[0]), 'r+')

g.write('Quality filter (%)\n')
QF_ratio = extract_QF(QF_log)
g.write('%.2f\n\n'%(QF_ratio))

g.write('Demultiplexing (%)\n')
R1_ratio, R2_ratio = extract_demultiplexing(demultiplexing_log)
g.write('R1 with adaptor\t%.2f\nR2 with adaptor\t%.2f\n\n'%(R1_ratio, R2_ratio))

g.write('With adaptor and long enough R2 (%)\n')
writtenLine = 'All_adaptors\t'
read_counts, prior_read_counts = extract_filter_short(dem_files, prior_dem_file)
for read_count in read_counts:
    ratio = read_count/prior_read_counts*100
    writtenLine = writtenLine + '%.2f'%(ratio) + ','
writtenLine = writtenLine[:-1] + '\n'
g.write(writtenLine)
total = np.sum(np.array(read_counts))/prior_read_counts * 100
g.write('Total\t%.2f\n\n'%(total))

g.write('Alignment rate (%)\n')
writtenLine = 'All_adaptors\t'
mapping_ratios = extract_mapping(mapping_logs)
for mapping_ratio in mapping_ratios:
    writtenLine = writtenLine + '%.2f'%(mapping_ratio) + ','
writtenLine = writtenLine[:-1] + '\n'
g.write(writtenLine)
average = np.mean(np.array(mapping_ratios))
g.write('Average\t%.2f\n'%(average))
std = np.std(np.array(mapping_ratios))
g.write('STD\t%.2f\n'%(std))
g.close()
