#this file performs demultiplexing and trimming of adapters and bad quality bases on paired end fastq files.

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

rule demultiplexing_paired:
    input:
        r1 = f"{DATA_DIR}/{{sample}}/{{sample}}_L4_1.fq.gz",
        r2 = f"{DATA_DIR}/{{sample}}/{{sample}}_L4_2.fq.gz"
    output:
        "results/1_DemultiplexTrim/demux_{sample}.touch"
    threads: 2
    log:
        stdout = "logs/1_DemultiplexTrim/{sample}/{sample}_demux.stdout",
        stderr = "logs/1_DemultiplexTrim/{sample}/{sample}_demux.stderr"
    shell:
        """
        cutadapt -j {threads} -e 0.15 -g ^file:{BARCODES} -o {{name}}_{wildcards.sample}_1.fastq -p {{name}}_{wildcards.sample}_2.fastq {input.r1} {input.r2}
        touch results/1_DemultiplexTrim/demux_{wildcards.sample}.touch
        """

checkpoint demux:
    input:
        did = expand("results/1_DemultiplexTrim/demux_{sample}.touch", sample = config.get("samples").keys())
    output:
        passed = "results/1_DemultiplexTrim/passed_demux.touch"
    shell:
        """
        touch results/1_DemultiplexTrim/passed_demux.touch
        echo "Checkpoint Demux"
        """

#Checkpoint method, allows for re-evaluation of the DAG after demultiplexing_paired
def getBarcodedFQ(wildcards):
    #Forces checkpoint
    check = checkpoints.demux.get(**wildcards).output[0]
    samplename = wildcards.sample
    barcodes = getBarcodes(BARCODES)
    return expand("{Barcode}_{sample}_{pair}.fastq", sample = samplename, Barcode = barcodes, pair = [1,2])

#Occasionally the paired read might have a barcode for a different file, unpairing the reads -> bbtools re-pair.
rule repair:
    input:
        list = getBarcodedFQ
    output:
        "results/1_DemultiplexTrim/repair_{sample}.touch"
    params:
        script = srcdir("../bbmap/repair.sh")
    log:
        stdout = "logs/1_DemultiplexTrim/{sample}/{sample}_repair.stdout",
        stderr = "logs/1_DemultiplexTrim/{sample}/{sample}_repair.stderr"
    shell:
        """
        echo "List is"
        echo {input.list}
        ../bbmap/repair.sh in={input.list[0]} in2={input.list[1]} out=fixed_{input.list[0]} out2=fixed_{input.list[1]} outs=unpaired_{wildcards.sample}.fastq
        touch results/1_DemultiplexTrim/repair_{wildcards.sample}.touch
        """

rule DidDemux:
    input:
        samples = expand("results/1_DemultiplexTrim/demux_{sample}.touch", sample = config.get("samples").keys())
#        repaired = expand("results/1_DemultiplexTrim/repair_{sample}.touch", sample = config.get("samples").keys())
    output:
        "results/1_DemultiplexTrim/DidDemux.touch"
    shell:
        """
        mv *.fastq results/1_DemultiplexTrim/
        touch results/1_DemultiplexTrim/DidDemux.touch
        """