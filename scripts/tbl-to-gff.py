#!/usr/bin/env python
import argparse
import re

def attr_formatter(attrs):
    s = ""
    for k in attrs:
        v = attrs[k]
        if v == "-":
            if k == "ID" and "RNA_name" in attrs:
                v = attrs["RNA_name"]
            else:
                v = "None"
        s += f"{k}={v};"
    return s

def main():
    parser = argparse.ArgumentParser(description="Convert Infernal hits to gff format")
    parser.add_argument('--input','-i', help="Infernal tabular format", required=True)
    parser.add_argument('--output','-o', help="Output bed file",required=True)
    args = parser.parse_args()

    fout = open(args.output,"w")
    with open(args.input) as fin:
        for line in fin:
            if line.startswith("#"):
                continue
            line = line.strip()
            if len(line) == 0:
                continue
            # 0.  target seq name
            # 1.  target accession
            # 2.  query profile name 
            # 3.  query profile acession
            # 4.  model type hmm/cm
            # 5.  model from
            # 6.  model to
            # 7.  seq from 
            # 8.  seq to
            # 9.  strand -
            # 10. trunc  -
            # 11. pass 6
            # 12. gc 0.50
            # 13. bias 0.0 
            # 14. score 72.2
            # 15. E value 6.4e-18
            # 16. inc ! 
            # 17. target description flag=1 multi=2.0000 len=538
            #  k119_23629           -         5S_rRNA              RF00001   hmm        3      117      265      151      -     -    6 0.50   0.0   72.2   6.4e-18 !   flag=1 multi=2.0000 len=538
            fields = re.split("\s+",line)
            cm_name, cm_id, seq_id = fields[2], fields[3], fields[0]
            mdl_start, mdl_end = int(fields[5]),int(fields[6])
            seq_start, seq_end = int(fields[7]),int(fields[8])
            strand = fields[9]
            trunc = fields[10]       
            gc = fields[12]
            bias = fields[13]
            score = fields[14]
            e_value = fields[15]
            if seq_start == seq_end:
                continue
            if seq_start > seq_end:
                seq_start, seq_end = seq_end, seq_start
                assert strand == "-", line
            else:
                assert strand == "+", line
            attrs = {"ID":cm_id,
                     "RNA_name":cm_name,
                     "trunc":trunc,
                     "gc":gc,
                     "bias":bias,
                     "e_value":e_value}
            print(seq_id, 
                  "Infernal","ncRNA", 
                  seq_start,seq_end,score,
                  fields[9], "0",
                  attr_formatter(attrs), file=fout, sep="\t")

    fout.close()


if __name__ == "__main__":
    main()

