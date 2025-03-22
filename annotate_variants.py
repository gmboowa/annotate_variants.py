import argparse
import os
from Bio import SeqIO, Seq
from Bio.Seq import MutableSeq

# Set up command-line argument parsing
parser = argparse.ArgumentParser(description="Annotate a VCF file with gene and variant effect from a GenBank reference.")
parser.add_argument("-i", "--input_vcf", required=True, help="Input VCF file")
parser.add_argument("--ref", "--reference_gb", required=True, help="Reference GenBank file")
parser.add_argument("-o", "--output_vcf", required=True, help="Output annotated VCF file")
args = parser.parse_args()

# Validate input files
if not os.path.isfile(args.input_vcf):
    raise FileNotFoundError(f"Error: VCF file '{args.input_vcf}' not found.")
if not os.path.isfile(args.ref):
    raise FileNotFoundError(f"Error: GenBank file '{args.ref}' not found.")

# Parse GenBank to extract gene & CDS regions
gene_annotations = {}
cds_features = []

record = SeqIO.read(args.ref, "genbank")
sequence = record.seq

for feature in record.features:
    if feature.type == "gene":
        gene_id = feature.qualifiers.get("gene", ["Unknown"])[0]
        for pos in range(int(feature.location.start) + 1, int(feature.location.end) + 1):
            gene_annotations[pos] = gene_id
    elif feature.type == "CDS":
        cds_features.append(feature)

# Function to check if variant causes missense or synonymous change
def get_variant_effect(pos, ref_allele, alt_allele):
    # Only handle SNPs (length 1 substitutions)
    if len(ref_allele) != 1 or len(alt_allele) != 1:
        return gene_annotations.get(pos, "Intergenic"), "Non-SNP"

    for cds in cds_features:
        start = int(cds.location.start) + 1
        end = int(cds.location.end)
        strand = cds.location.strand
        if start <= pos <= end:
            gene_id = cds.qualifiers.get("gene", ["Unknown"])[0]
            coding_seq = cds.extract(sequence)
            codon_pos = (pos - start) if strand == 1 else (end - pos)
            codon_index = codon_pos // 3
            base_index = codon_pos % 3
            if codon_index >= len(coding_seq) // 3:
                return gene_id, "Non-coding"

            codon = coding_seq[codon_index * 3: codon_index * 3 + 3]
            if len(codon) != 3:
                return gene_id, "Non-coding"

            ref_codon = MutableSeq(str(codon))
            if strand == 1:
                ref_codon[base_index] = alt_allele
            else:
                comp_alt = Seq.Seq(alt_allele).complement()
                ref_codon[2 - base_index] = str(comp_alt)[0]  # Use only first character

            new_codon = Seq.Seq(str(ref_codon))

            try:
                ref_aa = codon.translate()
                alt_aa = new_codon.translate()
            except Exception:
                return gene_id, "TranslationError"

            if ref_aa == alt_aa:
                return gene_id, "Synonymous"
            else:
                return gene_id, "Missense"

    return gene_annotations.get(pos, "Intergenic"), "Non-coding"

# Annotate the VCF file
vcf_data = []
with open(args.input_vcf, "r") as vcf:
    for line in vcf:
        if line.startswith("#"):
            if line.startswith("#CHROM"):
                vcf_data.append(line.strip() + "\tGene_ID\tVariant_Effect")
            else:
                vcf_data.append(line.strip())
        else:
            columns = line.strip().split("\t")
            if len(columns) < 5:
                continue
            chrom, pos, vid, ref, alt = columns[:5]
            pos = int(pos)
            gene_id, effect = get_variant_effect(pos, ref, alt)
            annotated_line = line.strip() + f"\t{gene_id}\t{effect}"
            vcf_data.append(annotated_line)

# Save the annotated VCF
with open(args.output_vcf, "w") as out:
    for entry in vcf_data:
        out.write(entry + "\n")

print(f"Annotated VCF with gene and variant effect saved to: {args.output_vcf}")
