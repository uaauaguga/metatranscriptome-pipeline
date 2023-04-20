#!/usr/bin/env python
from collections import defaultdict

import argparse
import os
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('extract orphan transcript')

def main():
    parser = argparse.ArgumentParser(description='extract orpha  transcript and group by igr')
    parser.add_argument('--input','-i',type=str, required=True, help="input annotated intervals")
    parser.add_argument('--output','-o',type=str, required=True, help="output igr")
    parser.add_argument('--min-distance', '-md', type=int, default = 16, help="this length flanking CDS are considered as downsream/leader")
    args = parser.parse_args()
    

    fout = open(args.output, "w")
    igrs = defaultdict(list)
    with open(args.input) as f:
        for line in f:
            fields = line.strip().split("\t")
            #k119_7	47	248	STRG.1.1	21.950249	-	.|<	inf,27	concordant	downstream	.|None
            if "," not in fields[7]:
                continue
            if "inf" in fields[7]:
                continue
            ud, dd = fields[7].split(",")       
            #ud = -1 if ud == "inf" else ud
            #dd = -1 if dd == "inf" else dd
    
            ud, dd = int(ud), int(dd) 
            if min(ud, dd) < args.min_distance:
                continue
            seq_id, start, end, name, coverage, strand = fields[:6]
            start, end = int(start), int(end)
            length = end - start
            start, end = start - ud, end + dd
            iv = (seq_id, start, end, strand)
            s  = ud
            igrs[iv].append((name, s, length, float(coverage)))

    for iv in sorted(list(igrs.keys()),key=lambda x:(x[0],x[1])):
        seq_id, start, end, strand = iv
        for attrs in igrs[iv]:
            name, s, L, coverage = attrs
            print(seq_id, start, end, name, coverage, strand, s, L ,sep="\t",file=fout)
    fout.close()
        
if __name__ == "__main__":
    main()
