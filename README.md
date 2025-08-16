# SV Pipeline

A Snakemake pipeline for alignment of whole-genome sequencing (WGS) data and structural variant (SV) calling.

## Steps
1. **Reference preparation** (BWA + samtools + GATK dict).
2. **Read alignment** with BWA-MEM → sorted, indexed BAM.
3. **Structural variant calling** using Delly.
4. **Filtering & export** → CSV with useful fields (CHROM, START, END, SIZE, QUAL, FILTER, Variant_Type).

## Requirements
- Snakemake
- BWA
- samtools
- GATK
- Delly
- Python ≥3.8 (with pandas)

## Usage
1. Edit `config.yaml` with:
   - Reference genome (`reference`)
   - Paired-end FASTQ files (`fastq`)

2. Run the pipeline:
   ```bash
   snakemake --cores 4
