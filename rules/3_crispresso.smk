#This rule file performs crispresso on paired fastqc files seperated by barcode in the Barcode_SampleName_1.fastq format

BARCODES = config['demultiplex_barcode_tsv']
AMPLICONREF = config['AmpliconReference']
GUIDESEQ = config['GuideSequences']

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

#Opens FASTA format amplicon file, takes sequences/names and returns strings in CRISPResso's preferred format x,y,z.
def getAmpliconSeqs(fastaFile):

    AmpliconReferenceNames = []
    AmpliconReferenceSequences = []
    GuideSequences = []

    lines = [line.strip() for line in open(fastaFile).readlines()]
    for line in lines:
        if '>' in line:
            ampliconName = line.replace('>','')
            AmpliconReferenceNames.append(ampliconName)
            #Retrieve corresponding guide sequence
            guides = [l.strip() for l in open(GUIDESEQ).readlines()]
            for guideLine in guides:
                if ampliconName in guideLine:
                    guideSeq = guideLine.split('\t')[1]
                    GuideSequences.append(guideSeq)
        if not '>' in line:
            AmpliconReferenceSequences.append(line)
    if len(AmpliconReferenceSequences) != len(GuideSequences):
        print("ERROR: unequal lengths for AmpliconReference and GuideSequences")

    #Compose csv 
    ARN = ""
    ARS = ""
    GS = ""
    for i in range(0, len(AmpliconReferenceSequences)):
        ARN = ARN + AmpliconReferenceNames[i] + ","
        ARS = ARS + AmpliconReferenceSequences[i] + ","
        GS = GS + GuideSequences[i] + ","
    return [ARN[0:-1],ARS[0:-1],GS[0:-1]]

rule crispresso:
    input:
        r1 = "results/1_DemultiplexTrim/{bc}_{sample}_1.fastq",
        r2 = "results/1_DemultiplexTrim/{bc}_{sample}_2.fastq"
    output:
        "results/3_crispresso/crispresso_{bc}_{sample}.touch"
    threads: 2
    params:
        AR = getAmpliconSeqs(AMPLICONREF)
    log:
        stdout = "logs/3_crispresso/{sample}/{bc}_{sample}_crispresso.stdout",
        stderr = "logs/3_crispresso/{sample}/{bc}_{sample}_crispresso.stderr"
    shell:
        """
        CRISPResso --expand_allele_plots_by_quantification -r1 {input.r1} -r2 {input.r2} -a {params.AR[1]} -an {params.AR[0]} -g {params.AR[2]} -o results/3_crispresso/
        touch results/3_crispresso/crispresso_{wildcards.bc}_{wildcards.sample}.touch
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
        touch results/3_crispresso/align_{wildcards.bc}_{wildcards.sample}.touch
        """

rule did_crispresso:
    input:
        crp = expand("results/3_crispresso/crispresso_{bc}_{sample}.touch", sample = config.get("samples").keys(), bc = getBarcodes(BARCODES)),
        algn = expand("results/3_crispresso/align_{bc}_{sample}.touch", sample = config.get("samples").keys(), bc = getBarcodes(BARCODES))
    output:
        "results/3_crispresso/finished_crispresso.touch"
    shell:
        """
        mv sorted_*.bam results/3_crispresso/
        mv sorted_*.bam.bai results/3_crispresso/
        rm *.sam
        rm *.bam
        CRISPRessoAggregate --name "ClockWork" --prefix results/3_crispresso/ --suffix _2
        touch results/3_crispresso/finished_crispresso.touch
        """