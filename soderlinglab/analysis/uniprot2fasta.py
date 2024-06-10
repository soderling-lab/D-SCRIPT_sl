import pandas as pd 
import numpy as np 
import os
import re
from Bio import SeqIO
import sys
import requests
from tqdm import tqdm
from io import StringIO

def get_fasta(accession):
    url = f"https://www.uniprot.org/uniprot/{accession}.fasta"
    response = requests.get(url)
    if response.status_code == 200:
        response = StringIO(response.text)
        rec = next(SeqIO.parse(response, "fasta"))
        return rec.seq
    else:
        return 'No data found'

def clean_uniprotfasta(record):
    recid = record.id.split('|')[1]
    fasta = str(record.seq)
    return recid, fasta

myfasta = os.path.expanduser('~/Desktop/D-SCRIPT/data/seqs/uniprots.fasta')
if os.path.exists(myfasta) is False:
    FASTA_FILE = f'data/uniprot_sprot.fasta'
    with open(FASTA_FILE, 'r') as f, open(myfasta, 'w') as n:
        for rec in SeqIO.parse(f, 'fasta'):
            rid, fasta = clean_uniprotfasta(rec)
            if len(fasta) > 2000:
                fasta = fasta[:2000]
            n.write(f'>{rid}\n{fasta}\n')

PAIRS = sys.argv[1] ## INPUT DATA
## read pairs
pairs = pd.read_csv(PAIRS, sep='\t', header=None, names = ['Item1', 'Item2'])

sequences = {}
uniqueacc = set(pairs['Item1'].values).union(set(pairs['Item2'].values))
with open(myfasta, 'r') as f:
    sequences = {rec.id: str(rec.seq) for rec in SeqIO.parse(f, 'fasta')}

# Prepare to fetch missing sequences
missing_sequences = []
for item in tqdm(uniqueacc):
    if item not in sequences:
        oseq = get_fasta(item)
        if oseq != 'No data found':
            missing_sequences.append((item, oseq))

# Write all missing sequences at once
with open(myfasta, 'a') as f:
    for item, seq in missing_sequences:
        if len(seq) > 2000:
            seq = seq[:2000]
        f.write(f'>{item}\n{seq}\n')

