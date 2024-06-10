import pandas as pd 
import numpy as np 
import os
import re
from Bio import SeqIO
import sys
import argparse
import itertools
import math 

def str2bool(value):
    if isinstance(value, bool):
        return value
    if value.lower() in {'true', 't', 'yes', 'y', '1'}:
        return True
    elif value.lower() in {'false', 'f', 'no', 'n', '0'}:
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def main():
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    currdir = os.getcwd()
    parent = os.path.dirname(currdir)
    gparent = os.path.dirname(parent)
    g2parent = os.path.dirname(gparent)
    parser = argparse.ArgumentParser(description="Pls parse an Excel file and an optional "
                                    "protein regulation (up/down) parameter.")

    # Add arguments
    parser.add_argument('data_file', help='The path to the Excel file to process.')

    parser.add_argument('-r', '--regulation', help='Which regulation are you focused on?', default=None)
    parser.add_argument('--allvall', help='Run All proteins vs All? Or just POI vs found? T/F', type=str2bool, 
                       )

    # Parse arguments
    args = parser.parse_args()
    READ_FILE = args.data_file
    REGULATION = args.regulation
    ALLVALL = args.allvall
    sheetnames = pd.ExcelFile(READ_FILE).sheet_names
    # for i, sh in enumerate(sheetnames):
    #     if 'normalized' in sh.lower(
    #         NORMSHEET= sh
    #         break
    DATASHEET = sheetnames[-1]
    base = os.path.basename(READ_FILE)
    data = pd.read_excel(READ_FILE, sheet_name=DATASHEET)
    # stats = pd.read_excel(READ_FILE, sheet_name=STATSHEET)
    datastats = data[['Accession', 'Stats']]

    ## Only UP!
    if REGULATION:
        datastats = datastats[datastats['Stats'].str.lower() == REGULATION.lower()]
        print(f'\nFiltering to {REGULATION}: {datastats.shape} remains.', file=sys.stderr)
    else:
        REGULATION = ''
    POI, POIA, ORG = 'Cnksr2', 'Q80YA9', 'Mouse' 
    pairs = [(acc_num.split(';')[0], POIA) for acc_num in datastats['Accession'].values]
    if ALLVALL:
        foundprots = [acc.split(';')[0] for acc in datastats['Accession'].values]
        ## permutation 12 is diff than 21 (order matters). combination: order does not matter!
        print(f'Given {len(foundprots)} unique proteins, there will be {math.comb(len(foundprots), 2)} combinations added',
              file=sys.stderr)
        #### len(foundprots)! / 2!(3-2)!. binomial coefficient C(totalprotslength, combinationlength) = n!/ (r!(n-r)!)
        pairs.extend(list(itertools.combinations(foundprots, 2)))
        print(f'after addition: {len(pairs)}\n', file=sys.stderr)
        file = f'{gparent}/data/pairs/{base.split(".xlsx")[0]}{REGULATION.lower()}_ALLvALL.tsv'
        with open(file, 'w') as f:
            for p in pairs:
                f.write(p[0] + '\t' + p[1] + '\n')  # Correctly formats each line in the TSV file
    else:
        file = f'{gparent}/data/pairs/{base.split(".xlsx")[0]}{REGULATION.lower()}.tsv'
        with open(f'{gparent}/data/pairs/{base.split(".xlsx")[0]}{REGULATION.lower()}.tsv', 'w') as f:
            for p in pairs:
                f.write(p[0] + '\t' + p[1] + '\n')  # Correctly formats each line in the TSV file
    print(file)
    return file

if __name__ == "__main__":
    main()