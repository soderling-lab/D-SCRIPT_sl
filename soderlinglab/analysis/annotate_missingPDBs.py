import re
import time
import json
import zlib
from xml.etree import ElementTree
from urllib.parse import urlparse, parse_qs, urlencode
import requests
from requests.adapters import HTTPAdapter, Retry
import sys
import os 
import argparse
import pandas as pd
from Bio import SeqIO
import numpy as np 

parser = argparse.ArgumentParser(description="Pls input the pairs you will analyze.")
parser.add_argument('--pairs', help='Pairs to process (.tsv)')
args = parser.parse_args()

pairs = pd.read_csv(args.pairs, sep='\t', header=None, names = ['Item1', 'Item2'])
pairsl =np.unique(pairs.values.flatten())

my3difile = os.path.expanduser('~/Desktop/D-SCRIPT/data/3di/sodlab.fa')
fastafile = os.path.expanduser(f'~/Desktop/D-SCRIPT/data/seqs/uniprots.fasta')
sequences = {}
with open(my3difile, 'r') as f:
    sequences = {rec.id: str(rec.seq) for rec in SeqIO.parse(f, 'fasta')}

missingacc = [acc for acc in pairsl if acc not in sequences]

with open(fastafile, 'r') as fasta:
    for rec in SeqIO.parse(fastafile, 'fasta'):
        if rec.id in missingacc:
            length = len(rec.seq)
            with open(my3difile, 'a') as difile:
                difile.write(f'>{rec.id}\n{"X"*length}\n')



# for x in missingacc:
#     with open(fastafile, 'r') as fasta:
#         for rec in SeqIO.parse(fastafile, 'fasta'):
#             if x == rec.id:
#                 length = len(rec.seq)
#                 with open(my3difile, 'a') as difile:
#                     difile.write(f'>{x}\n{'X'*length}\n')