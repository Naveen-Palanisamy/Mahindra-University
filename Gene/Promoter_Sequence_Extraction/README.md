#  Promoter Region Extractor Using BioPython + NCBI Entrez

This script reads a list of gene names from a file (`geneList.txt`) and uses the **NCBI Entrez API** (via BioPython) to:

1. Retrieve each gene's genomic coordinates  
2. Extract the **promoter region** (2kb upstream of the transcription start site, or TSS)  
3. Save the coordinates to a `.tsv` file and promoter sequences to a `.fasta` file

---

## Input

### `geneList.txt`

- A plain text file with one human gene symbol per line
- Example:
TP53
BRCA1
EGFR

##Run the script
python filename.py

---

## Output

### 1. `geneCoordinates.tsv`

A tab-separated file with the following columns:

| Gene | Chromosome | Start | End | Strand |
|------|------------|-------|-----|--------|

This contains the parsed genomic information for each gene from NCBI.

---

### 2. `promoterSequences.fasta`

- A FASTA-formatted file containing the upstream **promoter regions** for each gene.
- By default, it fetches **2000 base pairs upstream** of the gene's start site.
- The FASTA headers are named like: `TP53_promoter`, `BRCA1_promoter`, etc.

---

## ⚙️ Requirements

- Python 3.x
- [BioPython](https://biopython.org/)

Install BioPython using pip:
```bash
pip install biopython
