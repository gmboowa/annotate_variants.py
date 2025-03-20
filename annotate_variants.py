import argparse
import os
from Bio import SeqIO

# Set up command-line argument parsing
parser = argparse.ArgumentParser(description="Annotate a VCF file with gene information from a GenBank reference.")
parser.add_argument("-i", "--input_vcf", required=True, help="Input VCF file")
parser.add_argument("--ref", "--reference_gb", required=True, help="Reference GenBank file")
parser.add_argument("-o", "--output_vcf", required=True, help="Output annotated VCF file")
args = parser.parse_args()

# Validate input files
if not os.path.isfile(args.input_vcf):
    raise FileNotFoundError(f"Error: VCF file '{args.input_vcf}' not found.")
if not os.path.isfile(args.ref):
    raise FileNotFoundError(f"Error: GenBank file '{args.ref}' not found.")

# Step 1: Parse the GenBank file to extract gene locations and IDs
gene_annotations = {}

for record in SeqIO.parse(args.ref, "genbank"):
    for feature in record.features:
        if feature.type == "gene":
            start = int(feature.location.start) + 1  # Convert to 1-based indexing
            end = int(feature.location.end)
            gene_id = feature.qualifiers.get("gene", ["Unknown"])[0]  # Extract gene name or assign "Unknown"
            
            # Store gene info in a dictionary
            for pos in range(start, end + 1):
                gene_annotations[pos] = gene_id

# Step 2: Read and annotate the VCF file
vcf_data = []
with open(args.input_vcf, "r") as vcf:
    for line in vcf:
        if line.startswith("#"):  # Preserve header lines
            vcf_data.append(line.strip())
        else:
            columns = line.strip().split("\t")
            chrom, pos, vid, ref, alt = columns[:5]
            pos = int(pos)

            # Get corresponding Gene ID or mark as "Intergenic"
            gene_id = gene_annotations.get(pos, "Intergenic")
            
            # Append the new column to the VCF entry
            annotated_line = line.strip() + f"\t{gene_id}"
            vcf_data.append(annotated_line)

# Step 3: Write the Annotated VCF File
with open(args.output_vcf, "w") as output:
    for entry in vcf_data:
        output.write(entry + "\n")

print(f"Annotated VCF file saved as: {args.output_vcf}")