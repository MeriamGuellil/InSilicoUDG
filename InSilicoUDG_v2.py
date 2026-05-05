import sys
import os
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description="""in silico UDG using FASTQ reads from a mapdamage recalibrated bam file""")
parser.add_argument('-f',metavar='FASTQ file', dest='FASTQ_R', required=True, type=str, help='FASTQ file containing recalibrated reads mapping to the reference sequence')
parser.add_argument('-o',metavar='Output File', dest='out_id', required=True, type=str, help='Output file')
parser.add_argument('-p5',metavar='ratio of read positions at the 5 prime end to be investigated', dest='p5_cutoff', required=False,  default="0.20", type=float, help='ratio of read positions at the 5 prime end to be investigated')
parser.add_argument('-p3',metavar='ratio of read positions at the 3 prime end to be investigated', dest='p3_cutoff', required=False, default="0.20", type=float, help='ratio of read positions at the 3 prime end to be investigated')
parser.add_argument('-z',metavar='maximum number of bases with BQ0 to consider read', dest='max_zero', required=False, default="3", type=int, help='maximum number of bases with BQ0 to consider read')
args= parser.parse_args()

cut=0
kept=0
total=0
discard=0

with open(args.out_id, "w") as output:

    FASTQ = list(SeqIO.parse(args.FASTQ_R, "fastq"))

    for record in SeqIO.parse(args.FASTQ_R, "fastq"): 
        VALUE_p5 = {}
        VALUE_p3 = {}
        READ = {}

        if 0 in record.letter_annotations["phred_quality"]: 
            qual = record.letter_annotations["phred_quality"]
            nb = qual.count(0) 
            LEN  = len(qual)

            if nb <= args.max_zero: 

                READ = {pos:qual[pos] for pos in range(len(qual))}

                p5_tresh = int(len(READ)*args.p5_cutoff) 
                p3_tresh = int(len(READ)*(1-args.p3_cutoff)) 

                VALUE_p5 = {p+1:bq for p,bq in READ.items() if bq == 0 and p <= p5_tresh} 
                VALUE_p3 = {p:bq for p,bq in READ.items() if bq == 0 and p >= p3_tresh} 

                if len(VALUE_p5.keys()) > 0:
                    p5_Cut = int(max(VALUE_p5.keys()))
                else:
                    p5_Cut = 0
                if len(VALUE_p3.keys()) > 0:
                    p3_Cut = int(min(VALUE_p3.keys()))
                else:
                    p3_Cut = LEN
            
                if len(VALUE_p5.keys()) > 0 or len(VALUE_p3.keys()) < LEN:
                    SeqIO.write(record[p5_Cut:p3_Cut], output, "fastq")
                    cut += 1
            else:
                discard += 1
        else:
            SeqIO.write(record, output, "fastq")
            kept += 1

perc=int(cut)/int(len(FASTQ))
print("Total Reads:  "+ str(len(FASTQ)))
print("Total Trimmed:  "+ str(cut))
print("Total Untrimmed:  "+ str(kept))
print("Total discarded:  "+ str(discard))
print("Percent Trimmed:  "+ str(perc))
