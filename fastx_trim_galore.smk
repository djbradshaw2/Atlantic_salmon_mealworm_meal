#snakemake -s 'filepath to snakefile' --cores #

import sys

prefixes = []

with open("manifest.csv") as in_fh:
    # skip the first line
    in_fh.readline()
    for line in in_fh:
        (sample_name, path, direction) = line.split(",")
        (prefix, sep, filename) = path.rpartition("_")
        prefixes.append(prefix)
        # skip the next line as they are all in pairs
        in_fh.readline()

print("Found {} samples.".format(len(prefixes)))




rule all:
    input:
        expand("{pathname}_R1_trim_fastqc.zip", pathname=prefixes),
        expand("{pathname}_R2_trim_fastqc.zip", pathname=prefixes)


rule fastqc:
    input:
        "{pathname}_R1.fastq",
        "{pathname}_R2.fastq"
    output:
        "{pathname}_R1_fastqc.zip",
        "{pathname}_R2_fastqc.zip"
    threads: 2
    shell:
        "fastqc -t 2 {input[0]} {input[1]}"


rule trim_primers:
    input:
        "{pathname}_R1.fastq",
        "{pathname}_R2.fastq",
        "{pathname}_R1_fastqc.zip",
        "{pathname}_R2_fastqc.zip"
    output:
        "{pathname}_R1_val_1.fq",
        "{pathname}_R2_val_2.fq"
    shell:
        "trim_galore -a X --fastqc -q 20 --length 20 --clip_R1 19 --clip_R2 18 --dont_gzip -o `dirname {input[0]}` --paired {input[0]} {input[1]}"

rule trim_forward_reads:
    input:
        "{pathname}_R1_val_1.fq"
    output:
        "{pathname}_R1_trim.fq.gz"
    shell:
        "fastx_trimmer -z -l 300 -i {input[0]} -o {output[0]}"

rule trim_reverse_reads:
    input:
        "{pathname}_R2_val_2.fq"
    output:
        "{pathname}_R2_trim.fq.gz"
    shell:
        "fastx_trimmer -z -l 300 -i {input[0]} -o {output[0]}"


rule fastqc_final:
    input:
        "{pathname}_R1_trim.fq.gz",
        "{pathname}_R2_trim.fq.gz"
    output:
        "{pathname}_R1_trim_fastqc.zip",
        "{pathname}_R2_trim_fastqc.zip"
    threads: 2
    shell:
        "fastqc -t 2 {input[0]} {input[1]}"
