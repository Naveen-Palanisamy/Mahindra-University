# Gene Coordinate Fetcher from NCBI using BioPython

This project contains a Python script that fetches **genomic coordinates** of human genes using the **NCBI Entrez API** and **BioPython**.

It takes a list of gene names from a text file (`geneList.txt`) and returns each gene's chromosome, start position, end position, and strand orientation. The result is saved in a TSV file: `geneCoordinates.tsv`.

---

## Features

- Automated querying of NCBI Gene database
- Parses chromosome location and strand from gene records
- Supports gene name lists from text file
- Handles both forward and reverse strand annotations
- Outputs results in a clean `.tsv` format

---

## Input

- A plain text file named `geneList.txt` in the same directory.
- Each line should contain a single gene symbol (e.g., `BRCA1`, `TP53`, `EGFR`).

**Example `geneList.txt`:**
  1. TBR1
  2. IL5
  3. CPPA
  4. SL3
