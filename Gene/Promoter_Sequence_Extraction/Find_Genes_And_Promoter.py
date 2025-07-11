import sys
import time
from Bio import Entrez
from Bio import SeqIO

# REQUIREMENTS:
# pip install biopython
# Ensure that geneList.txt is in the same directory.
# This script generates:
#   - geneCoordinates.tsv: Tabular output of gene positions
#   - promoterSequences.fasta: FASTA file of upstream promoter sequences

# Set your email for NCBI Entrez (REQUIRED by NCBI)
Entrez.email = 'your_email@example.com'

# Number of bases upstream of TSS (Transcription Start Site)
upstream = 2000

# Read gene names from file
genelist = []
with open(file='/home/usr/Downloads/geneList.txt', mode='r') as input_file:
    genelist = [x.strip() for x in input_file.readlines()]

# Output header and file for gene coordinates
header = ["Gene", "Chromosome", "Start", "End", "Strand"]
writer = open(file='geneCoordinates.tsv', mode='w')
writer.write('\t'.join(header) + "\n")

# Exit if no genes provided
if len(genelist) < 1:
    print("No gene names provided in input list.\n")
    sys.exit(1)

print("Input gene list is:", genelist)

locationResult = {}

# Step 1: Fetch gene coordinates from NCBI
for gene in genelist:
    qterm = f'{gene}[Preferred Symbol] AND human[ORGN]'
    print("Querying NCBI for:", gene)

    # Search gene ID
    handle = Entrez.esearch('gene', qterm)
    record = Entrez.read(handle)
    handle.close()

    result = [gene]

    for gene_id in record['IdList']:
        handle = Entrez.efetch(db='gene', id=gene_id, retmode="gene_table")
        lines = handle.readlines()
        name = start = end = strand = None

        for line in lines:
            if 'Annotation' in line:
                # Parse the annotation line
                strand = 1  # default is forward
                words = line.strip().replace('Annotation: ', '').split(' ')

                if "complement" in words:
                    # Gene is on reverse strand
                    [chr, name, nucid, location, _] = words
                    strand = 2  # 2 = reverse strand in Entrez
                else:
                    [chr, name, nucid, location] = words

                start, end = location.split('..')
                start = int(start.replace('(', ''))
                end = int(end.replace(',', '').replace(')', ''))

                # Store gene metadata
                locationResult[gene] = [name, nucid, start, end, strand]
                result.extend([name, start, end, strand])
        handle.close()

    # Write result to file
    writer.write('\t'.join([str(x) for x in result]) + "\n")
    time.sleep(5)  # delay for NCBI API

writer.close()

# Step 2: Fetch promoter sequences upstream of gene start
fastaWriter = open("promoterSequences.fasta", "w")

for gene in locationResult.keys():
    [name, nucid, start, end, strand] = locationResult[gene]

    # Calculate promoter region (upstream of start)
    pstart = max(1, int(start) - upstream)  # avoid negative coordinates
    pend = int(start) - 1

    print(f"Pulling promoter sequence for {gene} ({nucid}:{pstart}-{pend})")

    # Fetch sequence using NCBI efetch
    try:
        handle = Entrez.efetch(
            db="nucleotide",
            id=nucid,
            rettype="fasta",
            retmode="text",
            strand=strand,
            seq_start=pstart,
            seq_stop=pend
        )

        for record in SeqIO.parse(handle, "fasta"):
            record.id = f"{gene}_promoter"
            record.description = ""
            SeqIO.write(record, fastaWriter, "fasta")

        handle.close()
        time.sleep(5)

    except Exception as e:
        print(f"Failed to fetch promoter for {gene}: {e}")

fastaWriter.close()
