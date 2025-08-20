import sys
import time
from Bio import Entrez

# REQUIREMENTS:
# pip install biopython
# Ensure that geneList.txt is in the same directory.
# This script will generate a file called geneCoordinates.tsv with gene location data.

# Set your email for NCBI Entrez (this is required by NCBI to identify the user)
Entrez.email = 'yourmail@gmail.com'

#  Read gene list from a file
# The input file should have one gene name per line (e.g., TP53, BRCA1, EGFR)
genelist = []
with open(file='/home/usr/Downloads/geneList.txt', mode='r') as input_file:
    genelist = [x.strip() for x in input_file.readlines()]

# Set up the output file and write header
header = ["Gene", "Chromosome", "Start", "End", "Strand"]
with open(file='geneCoordinates.tsv', mode='w') as writer:
    writer.write('\t'.join(header) + "\n")

    #  Check if gene list is empty
    if len(genelist) < 1:
        print("No gene names provided in input list.\n")
        sys.exit(1)

    print("Input gene list is:", genelist)
    locationResult = {}

    #  Loop through each gene and query NCBI for genomic coordinates
    for gene in genelist:
        qterm = f'{gene}[Preferred Symbol] AND human[ORGN]'  # Entrez query
        print(f"Querying NCBI for: {gene}")
        
        # Step 1: Search for gene ID
        handle = Entrez.esearch(db='gene', term=qterm)
        record = Entrez.read(handle)
        handle.close()

        result = [gene]  # Prepare to store gene + coordinates

        # Step 2: Fetch data for each matching gene ID
        for gene_id in record['IdList']:
            handle = Entrez.efetch(db='gene', id=gene_id, retmode="gene_table")
            lines = handle.readlines()
            name = start = end = strand = None

            for line in lines:
                if 'Annotation' in line:
                    # 📍 Parse the annotation line to extract coordinates
                    strand = 1  # default is forward strand
                    words = line.strip().replace('Annotation: ', '').split(' ')
                    print("Parsed annotation words:", words)

                    # Check if on reverse strand
                    if "complement)" in words:
                        [chr, name, nucid, location, isComplement] = words
                        strand = -1
                    else:
                        [chr, name, nucid, location] = words

                    # Extract start and end positions
                    start, end = location.split('..')
                    start = int(start.replace('(', ''))
                    end = int(end.replace(',', '').replace(')', ''))

                    # Save results
                    locationResult[gene] = [name, nucid, start, end, strand]
                    result.extend([name, start, end, strand])
            handle.close()

        # Convert all values to string and write to file
        resultstr = [str(x) for x in result]
        writer.write('\t'.join(resultstr) + "\n")

        #  Wait before the next request to avoid overwhelming NCBI's servers
        time.sleep(5)
