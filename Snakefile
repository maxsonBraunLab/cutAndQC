# CutAndQC performs initial QC on CutAndTag projects
import glob
import os 
from pathlib import Path,PurePath,PurePosixPath
from collections import defaultdict

configfile: "src/config.yml"
fastq_dir=config["FASTQDIR"][0]
fastq_ext=config["FASTQEXT"][0]

# map samples to fastqs
def detect_samples(dir):
    '''Decect samples from fastq directory'''
    samps=defaultdict(list)
    files = os.listdir(dir)
    for f in files:
        if f.endswith(fastq_ext):
            s=f.split("R")[0].strip("_.")
            samps[s].append(os.path.join(dir,f))
    return samps 

samps=detect_samples(fastq_dir)
print("Found samples:")
for s in samps:
    print(f"{s}: {samps[s]}")

# samples and reads 
sampdict = detect_samples(fastq_dir)
reads=[Path(Path(f).stem).stem for f in os.listdir(fastq_dir) if f.endswith(fastq_ext)]

fastqScreenDict = {
'database': {
   'hg38': {
     'bowtie2': config["BOWTIE2"]["HG38"][0]},
   'mm10': {
     'bowtie2': config["BOWTIE2"]["MM10"][0]}, 
   'ecoli': {
     'bowtie2': config["BOWTIE2"]["ECOLI"][0]}, 
   'myco': {
     'bowtie2': config["BOWTIE2"]["MYCO"][0]}, 
 },
 'aligner_paths': {'bowtie2': 'bowtie2'}
}


rule all:
    input:
        expand("data/fastqc/{read}.html", read=reads),
        expand("data/fastq_screen/{read}.fastq_screen.txt", read=reads),
        expand(["data/aligned/{sample}.bam",
        "data/preseq/lcextrap_{sample}.txt",
        "data/dtools/fingerprint_{sample}.tsv",
        ], sample=sampdict.keys()),
        "data/multiqc/multiqc_report.html"

# fastqc for each read 
rule fastqc:
    input:
        "data/raw/{read}.fastq.gz"
    output:
        html="data/fastqc/{read}.html",
        zip="data/fastqc/{read}_fastqc.zip"
    log:
        "data/logs/fastqc_{read}.log"
    threads: 4
    wrapper:
        "0.49.0/bio/fastqc"

# detect contaminants
rule fastq_screen:
    input:
        "data/raw/{read}.fastq.gz"
    output:
        txt="data/fastq_screen/{read}.fastq_screen.txt",
        png="data/fastq_screen/{read}.fastq_screen.png"
    params:
        fastq_screen_config=fastqScreenDict,
        subset=100000,
        aligner='bowtie2'
    log:
        "data/logs/fastq_screen_{read}.log"
    threads: 8
    wrapper:
        "0.60.0/bio/fastq_screen"

# align samples to genome
rule bowtie2:
    input:
        lambda wildcards: sampdict[wildcards.sample]
    output:
        "data/aligned/{sample}.bam"
    log:
        err="data/logs/bowtie2_{sample}.err"
    conda:
        "envs/align.yml"
    threads: 8
    shell:
        "bowtie2 --local --very-sensitive-local "
        "--no-unal --no-mixed --threads {threads} "
        "--no-discordant --phred33 "
        "-I 10 -X 700 -x {config[BOWTIE2][HG38]} "
        "-1 {input[0]} -2 {input[1]} 2>{log.err} | samtools view -Sbh - > {output}"

rule collate:
    input:
        rules.bowtie2.output
    output:
        temp("data/aligned/{sample}.collate.bam")
    conda:
        "envs/preseq.yml"
    log:
        "data/logs/samtools_collate_{sample}.log"
    shell:
        "samtools collate -o {output} {input} > {log} 2>&1"

rule fixmate:
    input:
        rules.collate.output
    output:
        temp("data/aligned/{sample}.fixmate.bam")
    conda:
        "envs/preseq.yml"
    log:
        "data/logs/samtools_fixmate_{sample}.log"
    shell:
        "samtools fixmate -m {input} {output} > {log} 2>&1"

rule sort:
    input:
        rules.fixmate.output
    output:
        temp("data/aligned/{sample}.sort.bam")
    conda:
        "envs/preseq.yml"
    log:
        "data/logs/samtools_sort_{sample}.log"
    shell:
        "samtools sort -o {output} {input} > {log} 2>&1"

rule markdup:
    input:
        rules.sort.output
    output:
        "data/aligned/{sample}.sort.dedup.bam"
    conda:
        "envs/preseq.yml"
    log:
        "data/logs/samtools_markdup_{sample}.log"
    shell:
        "samtools markdup {input} {output} > {log} 2>&1"

rule index:
    input:
        rules.markdup.output
    output:
        "data/aligned/{sample}.sort.dedup.bam.bai"
    conda:
        "envs/preseq.yml"
    log:
        "data/logs/samtools_index_{sample}.log"
    shell:
        "samtools index {input} > {log} 2>&1"

rule preseq:
    input:
       rules.markdup.output
    output:
        "data/preseq/estimates_{sample}.txt"
    conda:
        "envs/preseq.yml"
    log:
        "data/logs/preseq_{sample}.log"
    shell:
        "preseq c_curve -B -P -o {output} {input} > {log} 2>&1" 

rule preseq_lcextrap:
    input:
        rules.markdup.output
    output:
        "data/preseq/lcextrap_{sample}.txt"
    conda:
        "envs/preseq.yml"
    log:
        "data/logs/preseq_{sample}.log"
    shell:
        "preseq lc_extrap -B -P -e 1000000000 -o {output} {input} > {log} 2>&1"
    

rule plotFinger:
    input:
        expand("data/aligned/{sample}.sort.dedup.bam.bai", sample=sampdict.keys())
    output:
        "data/dtools/fingerprint_{sample}.tsv"
    conda:
        "envs/dtools.yml"
    params:
        bam="data/aligned/{sample}.sort.dedup.bam"
    log:
        "data/logs/fingerprint_{sample}.log"
    shell:
        "plotFingerprint -b {params.bam} --smartLabels --outRawCounts {output}"

rule multiqc:
    input:
        expand("data/dtools/fingerprint_{sample}.tsv", sample=sampdict.keys()), directory("data/")
    output:
        "data/multiqc/multiqc_report.html"
    conda:
        "envs/multiqc.yml"
    log:
        "data/logs/multiqc.log"
    shell:
        "multiqc --force -o data/multiqc {input} > {log} 2>&1"

