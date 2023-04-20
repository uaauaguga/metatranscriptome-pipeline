#!/usr/bin/env python
import argparse
import pysam

def main():
    parser = argparse.ArgumentParser(description='infer strandness of MTX dataset by metaphlen marker gene mapping')
    parser.add_argument('--input', '-i', type=str, required=True, help='bam file produced by metaphlen')
    parser.add_argument('--output', '-o', type=str, required=True, help='fraction of reads aligned to forward strand and reverse strand')
    parser.add_argument('--min-mapping-quality','-mapq', type=int, default=0, help= "minimum mapping quality for a read to be considered")
    args = parser.parse_args()
    samfile = pysam.AlignmentFile(args.input, "rb")
    n_reverse, n_forward = 0, 0
    for read in samfile:
        if read.is_paired and read.is_read2:
            continue
        if read.is_unmapped:
            continue
        if read.mapping_quality < args.min_mapping_quality:
            continue
        if read.is_reverse:
            n_reverse += 1
        else:
            n_forward += 1
    with open(args.output,"w") as f:
        f.write("forward\treverse\tfraction\n")
        f.write(f"{n_forward}\t{n_reverse}\t{n_forward/(n_forward+n_reverse)}\n")
           




if __name__ == "__main__":
    main()
