#!/bin/bash

# activate env
source ~/micromamba/etc/profile.d/micromamba.sh
micromamba activate dscript

proteomic_data=~/Desktop/D-SCRIPT/soderlinglab/BioIDdata/Final_sltmm_Cnksr2_gad2.xlsx
regulation='Up'
allvall="True"

echo "Making pairs from your file: ${proteomic_data}"
pairsfile=$(python soderlinglab/analysis/make_pairs.py "$proteomic_data" --r "$regulation" \
--allvall "$allvall")

echo "Pairs created: Updating data/seqs/uniprots.fasta based on '$pairsfile'"
python soderlinglab/analysis/uniprot2fasta.py "$pairsfile"

echo "Extract 3di data"
out3di=$(dscript extract-3di data/pdbs/ data/3di/sodlab.fa 2>&1 | tr ' ' '\n' | tail -n 10)
echo "output from extract-3di '$out3di'"

echo "Renaming > line items to starting w accession"
sed -i '/z^>/ s/_/ /g' data/3di/sodlab.fa
sed -i '/^>/ s/-/ /g' data/3di/sodlab.fa

echo "adding missing accessions into 3di .fa file"
python soderlinglab/analysis/annotate_missingPDBs.py --pairs "$pairsfile"

echo "Running TT3D"
dscript predict --pairs "$pairsfile" --seqs data/seqs/uniprots.fasta \
--foldseek_fasta data/3di/sodlab.fa --device 0 --model \
pretrained_models/tt3d_v1.sav

# Find the most recent .tsv file and rename it
recent_tsv=$(ls -t *.tsv | head -n 2)
echo "Files created: $recent_tsv\n Now going to manipulate name!"

extractbasename=$(basename "$proteomic_data" .xlsx)
suf1="${extractbasename#*_}" ## splits after first underscore
suf2="${suf1#*_}"
if [[ "$allvall" == "True" ]]; then
    str2insert="${suf2}${regulation}_AllvAll"
else
    str2insert="${suf2}${regulation}"
fi

for tsv_file in $recent_tsv; do
    if [[ "$tsv_file" == *".predictions"* ]]; then
    ## syntax displays ${variable/pattern/replacement} replaces1st occurance of pattern
        new_name="${tsv_file/.predictions/_${str2insert}.predictions}"
        mv "$tsv_file" "$new_name"
        echo "Renamed $tsv_file to $new_name"
    else
        echo "File $tsv_file does not contain the pattern '.predictions'"
    fi
done

recent_tsv=$(ls -t | grep '\.predictions\.tsv$' | head -n 1)
if [ -n "$recent_tsv" ]; then
    proteins1=$(awk -F'\t' '{print $1}' "$recent_tsv" | tr '\n' ' ')
    proteins2=$(awk -F'\t' '{print $2}' "$recent_tsv" | tr '\n' ' ')
    echo "$proteins1" > accessions1.txt
    echo "$proteins2" > accessions2.txt

    col1Genes=$(python ~/Documents/UniprotIDmapping/write_tsv.py --idfile1 accessions1.txt \
    --idfile2 accessions2.txt --inputdb "UniProtKB_AC-ID" --outputdb "Gene_Name" --tsv_file "$recent_tsv")

    rm accessions1.txt
    rm accessions2.txt
    echo "Updated file: $recent_tsv with gene names and new values inserted."
else
    echo "No predictions.tsv file found."
fi