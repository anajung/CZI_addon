#!/usr/bin/env python

from Bio import SeqIO

def filter_fasta(combinedfadata):
    with open ("filtered_combined.fa", 'w') as out_handle:
        for record in SeqIO.parse(combinedfadata, "fasta"):
            if record.seq.count('N')<= 9000 and record.seq.count('-')<= 9000:
                out_handle.write(record.format("fasta"))

combinedfadata = "!{combinedfadata}"
filter_fasta(combinedfadata)