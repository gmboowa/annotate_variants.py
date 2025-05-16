import argparse
import os
import sys
from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord
from intervaltree import Interval, IntervalTree
import gzip

def parse_args():
    parser = argparse.ArgumentParser(description="Annotate VCF with gene and variant effect from GenBank")
    parser.add_argument("-i", "--input", required=True, help="Input VCF file")
    parser.add_argument("-r", "--reference", required=True, help="Reference GenBank file")
    parser.add_argument("-o", "--output", required=True, help="Output annotated VCF file")
    parser.add_argument("--compress", action="store_true", help="Compress output with gzip")
    parser.add_argument("-v", "--verbose", action="store_true", help="Show progress messages")
    return parser.parse_args()

class GenbankParser:
    def __init__(self, gb_file):
        self.genes = IntervalTree()
        self.cds_features = []
        self.translation_tables = {}
        
        try:
            for record in SeqIO.parse(gb_file, "genbank"):
                self.sequence = record.seq
                for feat in record.features:
                    self._process_feature(feat)
        except Exception as e:
            raise ValueError(f"Invalid GenBank file: {str(e)}")

    def _process_feature(self, feat):
        start = feat.location.start.position + 1  # 1-based
        end = feat.location.end.position
        strand = feat.location.strand
        
        if feat.type == "gene":
            gene_id = feat.qualifiers.get("gene", ["Unknown"])[0]
            self.genes.addi(start, end, (gene_id, strand))
            
        elif feat.type == "CDS":
            transl_table = int(feat.qualifiers.get("transl_table", ["11"])[0])
            gene_id = feat.qualifiers.get("gene", ["Unknown"])[0]
            self.cds_features.append({
                'start': start,
                'end': end,
                'strand': strand,
                'gene': gene_id,
                'transl_table': transl_table
            })

class VCFAnnotator:
    def __init__(self, gb_parser):
        self.gb = gb_parser
        
    def _get_coding_effect(self, pos, ref, alt, chrom):
        # Handle only SNPs for amino acid changes
        if len(ref) != 1 or len(alt) != 1:
            return "Non-SNP", "Intergenic"
            
        # Find overlapping genes
        gene_hits = self.gb.genes.at(pos)
        if not gene_hits:
            return "Intergenic", "Intergenic"
            
        # Check CDS regions
        for cds in self.gb.cds_features:
            if cds['start'] <= pos <= cds['end']:
                return self._calculate_amino_acid_change(pos, ref, alt, cds)
                
        return "Non-coding", gene_hits.pop().data[0]

    def _calculate_amino_acid_change(self, pos, ref, alt, cds):
        try:
            # Get relative position in CDS
            cds_pos = pos - cds['start']
            if cds['strand'] == -1:
                cds_pos = cds['end'] - pos
                
            # Extract codon (0-based in 3bp chunks)
            codon_num, offset = divmod(cds_pos, 3)
            codon_start = cds['start'] + codon_num*3
            codon_end = codon_start + 3
            
            # Get original codon
            codon = self.gb.sequence[codon_start-1:codon_end-1]
            if cds['strand'] == -1:
                codon = codon.reverse_complement()
                
            # Create mutated codon
            mut_codon = list(str(codon))
            mut_codon[offset] = alt if cds['strand'] == 1 else Seq.Seq(alt).reverse_complement()[0]
            mut_codon = Seq.Seq(''.join(mut_codon))
            
            # Translate
            table = cds['transl_table']
            orig_aa = codon.translate(table=table)
            new_aa = mut_codon.translate(table=table)
            
            return ("Synonymous" if orig_aa == new_aa else "Missense"), cds['gene']
            
        except Exception:
            return "Translation-Error", cds['gene']

def main():
    args = parse_args()
    
    # Validate inputs
    for f in [args.input, args.reference]:
        if not os.path.exists(f):
            sys.exit(f"Error: Input file {f} not found")

    # Parse GenBank
    try:
        gb_parser = GenbankParser(args.reference)
    except Exception as e:
        sys.exit(f"GenBank parsing failed: {str(e)}")

    # Initialize annotator
    annotator = VCFAnnotator(gb_parser)
    
    # Handle compressed output
    opener = gzip.open if args.compress else open
    mode = 'wt' if args.compress else 'w'
    
    with opener(args.output, mode) as out:
        # Process VCF
        with gzip.open(args.input, 'rt') if args.input.endswith('.gz') else open(args.input) as vcf:
            for line in vcf:
                if line.startswith('##'):
                    if 'INFO' not in line:
                        out.write(line)
                    continue
                        
                if line.startswith('#CHROM'):
                    header = line.strip().split('\t')
                    header += ['##INFO=<ID=GENE,Number=1,Type=String,Description="Gene identifier">']
                    header += ['##INFO=<ID=EFFECT,Number=1,Type=String,Description="Variant effect">']
                    out.write('\t'.join(header) + '\n')
                    continue
                
                parts = line.strip().split('\t')
                if len(parts) < 5:
                    continue
                
                chrom, pos, = parts[0], int(parts[1])
                ref, alts = parts[3], parts[4].split(',')
                
                annotations = []
                for alt in alts:
                    effect, gene = annotator._get_coding_effect(pos, ref, alt, chrom)
                    annotations.append(f"GENE={gene};EFFECT={effect}")
                
                parts[7] = ';'.join([parts[7]] + annotations) if parts[7] != '.' else ';'.join(annotations)
                out.write('\t'.join(parts) + '\n')

if __name__ == "__main__":
    main()
