#This rule file performs crispresso on paired fastqc files seperated by barcode in the Barcode_SampleName_1.fastq format

BARCODES = config['demultiplex_barcode_tsv']
AMPLICONREF = config['AmpliconReference']
ALLELES = config['Alleles']
ANALYSIS = config['Analysis']
STRELKA = config['Strelka']

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

checkpoint conditional:
    input:
        did2 = "results/2_fastqc/multiqc_report.html"
    output:
        passed = "results/3_crispresso/passed_conditional.touch"
    shell:
        """
        touch results/3_crispresso/passed_conditional.touch
        echo "Checkpoint conditional"
        """

#Input function that uses config to use forward, reverse or paired fastq's for CRISPResso
def conditional_crispresso(wildcards):
    check = checkpoints.conditional.get(**wildcards).output[0]
    samplename = wildcards.sample
    lines = [line.strip() for line in open(ANALYSIS).readlines()]
    samples = config.get("samples").keys()
    bcs = getBarcodes(BARCODES)
    
    for line in lines:
        AmpliconName = line.split("\t")[0]
        AnalysisType = line.split("\t")[1]
        if AmpliconName in samplename:
            if AnalysisType == "Paired":
                #Use both fq in CRISPResso
                return expand("results/3_crispresso/crispressoPaired_{bc}_{sample}.touch", sample = samplename, bc = bcs)
            if AnalysisType == "Forward":
                #Use only forward fq in CRISPResso
                return expand("results/3_crispresso/crispressoForward_{bc}_{sample}.touch", sample = samplename, bc = bcs)
                
ruleorder: convert > sort > index

rule crispressoForward:
    input:
        r1 = "results/1_DemultiplexTrim/{bc}_{sample}_1.fastq",
        r2 = "results/1_DemultiplexTrim/{bc}_{sample}_2.fastq"
    output:
        "results/3_crispresso/crispressoForward_{bc}_{sample}.touch"
    params:
        script = srcdir("../scripts/runCRISPResso.py")
    threads: 2
    log:
        stdout = "logs/3_crispresso/{sample}/{bc}_{sample}_crispresso.stdout",
        stderr = "logs/3_crispresso/{sample}/{bc}_{sample}_crispresso.stderr"
    shell:
        """
        python {params.script} -a {ALLELES} -r1 {input.r1} -r2 {input.r2} -o results/3_crispresso/ -n {wildcards.bc}_{wildcards.sample}_Forward -t Forward
        touch results/3_crispresso/crispressoForward_{wildcards.bc}_{wildcards.sample}.touch
        """

rule crispressoPaired:
    input:
        r1 = "results/1_DemultiplexTrim/{bc}_{sample}_1.fastq",
        r2 = "results/1_DemultiplexTrim/{bc}_{sample}_2.fastq"
    output:
        "results/3_crispresso/crispressoPaired_{bc}_{sample}.touch"
    params:
        script = srcdir("../scripts/runCRISPResso.py")
    threads: 2
    log:
        stdout = "logs/3_crispresso/{sample}/{bc}_{sample}_crispresso.stdout",
        stderr = "logs/3_crispresso/{sample}/{bc}_{sample}_crispresso.stderr"
    shell:
        """
        python {params.script} -a {ALLELES} -r1 {input.r1} -r2 {input.r2} -o results/3_crispresso/ -n {wildcards.bc}_{wildcards.sample}_Paired -t Paired
        touch results/3_crispresso/crispressoPaired_{wildcards.bc}_{wildcards.sample}.touch
        """

rule align:
    input:
        r1 = "results/1_DemultiplexTrim/{bc}_{sample}_1.fastq",
        r2 = "results/1_DemultiplexTrim/{bc}_{sample}_2.fastq"
    output:
        "{bc}_{sample}.sam"
    threads: 2
    log:
        stdout = "logs/3_crispresso/{sample}/{bc}_{sample}_align.stdout",
        stderr = "logs/3_crispresso/{sample}/{bc}_{sample}_align.stderr"
    shell:
        """
        bwa mem -o {wildcards.bc}_{wildcards.sample}.sam {AMPLICONREF} {input.r1} {input.r2}
        """

rule convert:
    input:
        sam = "{bc}_{sample}.sam"
    output:
        "{bc}_{sample}.bam"
    threads: 2
    log:
        stdout = "logs/3_crispresso/{sample}/{bc}_{sample}_align.stdout",
        stderr = "logs/3_crispresso/{sample}/{bc}_{sample}_align.stderr"
    shell:
        """
        samtools view -b -o {wildcards.bc}_{wildcards.sample}.bam {input.sam}
        rm {input.sam}
        """

rule sort:
    input:
        "{bc}_{sample}.bam"
    output:
        "sorted_{bc}_{sample}.bam"
    threads: 2
    log:
        stdout = "logs/3_crispresso/{sample}/{bc}_{sample}_align.stdout",
        stderr = "logs/3_crispresso/{sample}/{bc}_{sample}_align.stderr"
    shell:
        """
        samtools sort -o sorted_{wildcards.bc}_{wildcards.sample}.bam {input}
        rm {input}
        """

rule index:
    input:
        sorted = "sorted_{bc}_{sample}.bam"
    output:
        "results/3_crispresso/align_{bc}_{sample}.touch"
    threads: 2
    log:
        stdout = "logs/3_crispresso/{sample}/{bc}_{sample}_align.stdout",
        stderr = "logs/3_crispresso/{sample}/{bc}_{sample}_align.stderr"
    shell:
        """
        samtools index {input.sorted}
        mv {input.sorted} results/3_crispresso/
        mv {input.sorted}.bai results/3_crispresso/
        touch results/3_crispresso/align_{wildcards.bc}_{wildcards.sample}.touch
        """

rule strelka:
    input:
        aligned = "results/3_crispresso/align_{bc}_{sample}.touch"
    output:
        "results/3_crispresso/strelka_{bc}_{sample}.touch"
    threads: 2
    shell:
        """
        {STRELKA}/bin/configureStrelkaGermlineWorkflow.py --exome --bam "results/3_crispresso/sorted_{wildcards.bc}_{wildcards.sample}.bam" --referenceFasta {AMPLICONREF} --runDir ./results/3_crispresso/Strelka/{wildcards.bc}_{wildcards.sample}/
        ./results/3_crispresso/Strelka/{wildcards.bc}_{wildcards.sample}/runWorkflow.py -m local -j 2
        touch results/3_crispresso/strelka_{wildcards.bc}_{wildcards.sample}.touch
        """

rule did_crispresso:
    input:
        crp = conditional_crispresso,
        align = expand("results/3_crispresso/align_{bc}_{sample}.touch", sample = config.get("samples").keys(), bc = getBarcodes(BARCODES)),
        strelka = expand("results/3_crispresso/strelka_{bc}_{sample}.touch", sample = config.get("samples").keys(), bc = getBarcodes(BARCODES)),
    output:
        "results/3_crispresso/finished_crispresso_{sample}.touch"
    shell:
        """
        CRISPRessoAggregate -p 10 --name "ClockWorkPaired" --prefix results/3_crispresso/ --suffix _Paired
        CRISPRessoAggregate -p 10 --name "ClockWorkForward" --prefix results/3_crispresso/ --suffix _Forward
        touch results/3_crispresso/finished_crispresso_{wildcards.sample}.touch
        """
