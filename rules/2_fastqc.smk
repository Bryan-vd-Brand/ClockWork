#This rule file performs fastqc/multiqc on the original data set.

import pandas as pd
from os import path
import os 

DATA_DIR = config['input_dir']
BARCODES = config['demultiplex_barcode_tsv']

#opens FASTA format barcode file, returns all barcode names for input filename generation
def getBarcodes(barcodeFile):
    print(F"Opening {barcodeFile} and constructing barcodeIDs")
    BarcodeIDs = []
    lines = [line.strip() for line in open(barcodeFile).readlines()]
    for line in lines:
        if '>' in line:
            barcodeID = line.replace('>','')
            BarcodeIDs.append(barcodeID)
    return BarcodeIDs

checkpoint fixed:
    input:
        did1 = "results/1_DemultiplexTrim/DidDemux.touch"
    output:
        passed = "results/2_fastqc/passed_fixed.touch"
    shell:
        """
        touch results/2_fastqc/passed_fixed.touch
        echo "Checkpoint fixed"
        """

#Checkpoint method, causes re-evaluation of the DAG/Jobs and returns all paired files for fastqc
def getDemultiplexed_FQpairs(wildcards):
    #Forces checkpoint
    check = checkpoints.fixed.get(**wildcards).output[0]
    samplename = wildcards.sample
    barcode = getBarcodes(BARCODES)
    #fastqc expects output directory to already exist
    if not os.path.isdir(f"results/2_fastqc/{samplename}/"):
        os.system(f"echo Did not find {samplename} dir ")
        os.system(f"mkdir results/2_fastqc/{samplename}/")


    return expand("results/1_DemultiplexTrim/{Barcode}_{sample}_{pair}.fastq", sample = samplename, Barcode = barcode, pair = [1,2])


rule fastqc_paired_2:
    input:
        fqs = getDemultiplexed_FQpairs
    output:
        outdir = directory("results/2_fastqc/{sample}/"),
        completed = "results/2_fastqc/fastqc_{sample}.touch"
    threads: 20
    log:
        stdout = "logs/2_fastqc/{sample}/{sample}_fastqc.stdout",
        stderr = "logs/2_fastqc/{sample}/{sample}_fastqc.stderr"
    shell:
        """
        mkdir "results/2_fastqc/{wildcards.sample}/"
        touch "results/2_fastqc/fastqc_{wildcards.sample}.touch"
        fastqc -o {output.outdir} -t {threads} {input.fqs} 1> {log.stdout} 2> {log.stderr}
        """

rule multiqc_2:
    input:
        dir = expand("results/2_fastqc/{sample}/", sample = config.get("samples").keys()),
        touch = expand("results/2_fastqc/fastqc_{sample}.touch", sample = config.get("samples").keys())
    output:
        "results/2_fastqc/multiqc_report.html"
    params:
        outdir = "results/2_fastqc"
    shell:
        """
        multiqc -v -f {input.dir} -o {params.outdir}
        """