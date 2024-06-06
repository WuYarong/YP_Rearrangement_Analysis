# YP_Rearrangement_Analysis
This script is compatible with both Python 2 and Python 3, having been tested on Python 2.7+ and Python 3.6 & 3.7. 


## Software

Nucleotide-Nucleotide BLAST 2.13.0+

Python 2.7+ or Python 3.6 & 3.7

## Python Dependencies
```
scipy  
numpy  
statsmodels  
```

## Installation Guide
```
pip install scipy  
pip install numpy  
pip install statsmodels  
```

## Usage
### Identification of RBP (rearrangement-related breakpoint) hotspots
```
python identify_RBP_hotspots.py.py input_file [--significance_level SIGNIFICANCE_LEVEL] [--apply_fdr] [--prefix PREFIX]
```

`
python identify_RBP_hotspots.py.py --apply_fdr --prefix test.RBP_hotspots.out example/test.RBP_cout.txt
`

### Copy number count of insertion sequences (IS)


#### step1: blastn mapping (output format: 6)
`
makeblastdb -in example/test.IS.fa -dbtype nucl  
blastn -query example/test.strain1.chr.fasta -db example/test.IS.fa -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -out test.strain1.IS.m6.out  
blastn -query example/test.strain2.chr.fasta -db example/test.IS.fa -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" -out test.strain2.IS.m6.out 
` 

#### step2: distinguish and count different IS elements for multiple strains
`
grep IS100 test.strain*.m6.out >IS100.blast.m6.merge.txt  
grep IS1541 test.strain*.m6.out >IS1541.blast.m6.merge.txt  
python IS_count_merge.py IS100.blast.m6.merge.txt >IS100.blast.m6.merge.stat.xls  
python IS_count_merge.py IS1541.blast.m6.merge.txt >IS1541.blast.m6.merge.stat.xls  
`

