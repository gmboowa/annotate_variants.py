# Annotate Variants with Gene Information

`annotate_variants.py` is a Python script that annotates a VCF file with gene information from a GenBank reference file. It matches variant positions to known genes & adds an extra column to indicate the corresponding Gene ID or labels it as **Intergenic**.

## Features
- Extracts **gene annotations** from a GenBank reference file.
- Matches **variant positions** in a VCF file to genes.
- Adds a **Gene_ID column** to the annotated VCF.
- Supports **command-line arguments** for flexibility.
- **Preserves** the original VCF structure and headers.

---

## Installation

### **1. Install Required Dependencies**
Make sure you have **Python 3+** and `biopython` installed.

```sh
pip install biopython
```

---

## Usage

Run the script with the following arguments:

```sh
python annotate_variants.py -i variants.vcf --ref reference.gbk -o annotated.vcf
```

| Argument | Description |
|----------|------------|
| `-i, --input_vcf` | Input VCF file with variant data |
| `--ref, --reference_gbk` | Reference GenBank file containing gene annotations |
| `-o, --output_vcf` | Output annotated VCF file |

---

## Example

### **1. Download a Reference GenBank File**
You can download a **GenBank file** from NCBI using:
```sh
esearch -db nucleotide -query MH910496.1 | efetch -format Genbank > reference.gbk
```

### **2. Run the Annotation Script**
```sh
python annotate_variants.py -i variants.vcf --ref reference.gbk -o annotated_variants.vcf
```

### **3. Output (Annotated VCF File)**
```
#CHROM  POS     ID      REF     ALT     Gene_ID     Variant_Effect
chr1    500     .       A       T       ABC1        Coding
chr1    1200    .       G       C       Intergenic  Intergenic
chr1    3000    .       T       A       XYZ2        Coding
```

---

## Troubleshooting
- **Error: `FileNotFoundError`** → Ensure the input files exist.
- **Error: `biopython not found`** → Install it using `pip install biopython`.
- **Incorrect gene matches?** → Check if the GenBank file contains complete gene annotations.

---

## License
This project is open-source and available under the **MIT License**.

---

## Contributions
Contributions are welcome! If you have suggestions, feel free to submit a **pull request**.

---
