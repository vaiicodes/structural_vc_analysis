#Genomic Pipeline for Structural Variant Calling with Delly

configfile: "config.yaml"

from pathlib import Path
import gzip, csv

SAMPLE = "HG002"
REF = config["reference"]
R1, R2 = config["fastq"]
DICT = str(Path(REF).with_suffix(".dict"))

# ---------- helper: VCF -> CSV ----------
#parse Delly's VCF output into the same CSV format
def sv_vcf_to_csv(in_vcf, out_csv):
    import pandas as pd
    import gzip

    rows = []
    with gzip.open(in_vcf, 'rt') as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            chrom, pos, vid, ref, alt, qual, filt, info = parts[:8]

            #Parse INFO field for SV details
            info_dict = {kv.split("=")[0]: kv.split("=")[1] for kv in info.split(";") if "=" in kv}
            svtype = info_dict.get("SVTYPE", "NA")
            end = int(info_dict.get("END", pos))
            size = abs(end - int(pos))

            rows.append({
                "CHROM": chrom,
                "START": int(pos),
                "END": end,
                "SIZE": size,
                "QUAL": float(qual) if qual != "." else None,
                "FILTER": filt,
                "Variant_Type": svtype
            })

    pd.DataFrame(rows).to_csv(out_csv, index=False)

# ---------- final outputs ----------
rule all:
    input:
        f"results/{SAMPLE}.sorted.bam",
        f"results/{SAMPLE}.sorted.bam.bai",
        f"results/{SAMPLE}.delly.vcf.gz",
        f"results/{SAMPLE}.variants.csv"

# ---------- index the reference ----------
rule index_ref:
    input:
        REF
    output:
        fai = REF + ".fai",
        dic = DICT
    shell:
        r"""
        samtools faidx {input}
        gatk CreateSequenceDictionary -R {input} -O {output.dic}
        """

# ---------- read alignment, sort, index ----------
rule align:
    input:
        ref = REF,
        idx = rules.index_ref.output,
        r1  = R1,
        r2  = R2
    output:
        bam = f"results/{SAMPLE}.sorted.bam",
        bai = f"results/{SAMPLE}.sorted.bam.bai"
    threads: 4
    shell:
        r"""
        mkdir -p results
        bwa mem -t {threads} -R '@RG\tID:{SAMPLE}\tSM:{SAMPLE}\tPL:ILLUMINA\tLB:lib1' \
            {input.ref} {input.r1} {input.r2} \
        | samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """

# ---------- Call structural variants with Delly ----------
rule call_delly:
    input:
        bam = f"results/{SAMPLE}.sorted.bam",
        ref = REF
    output:
        vcf = f"results/{SAMPLE}.delly.vcf.gz"
    threads: 4
    shell:
        r"""
        delly call -g {input.ref} {input.bam} -o {output.vcf}
        """

# ---------- Filter structural variants ----------
rule filter_variants:
    input:
        vcf=f"results/{SAMPLE}.delly.vcf.gz",
        ref=REF
    output:
        vcf=f"results/{SAMPLE}.filtered.vcf.gz"
    shell:
        r"""
        #simple filtering
        bcftools filter -O z \
            -e 'QUAL<20 || INFO/SVTYPE="BND"' \
            {input.vcf} > {output.vcf}
        tabix -p vcf {output.vcf}
        """

# ----------- Convert filtered VCF to CSV --------------
rule vcf_to_csv:
    input:
        vcf=f"results/{SAMPLE}.filtered.vcf.gz"
    output:
        csv=f"results/{SAMPLE}.variants.csv"
    run:
        sv_vcf_to_csv(input.vcf, output.csv)
