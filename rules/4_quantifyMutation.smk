#This rule takes in a directory structure created by CRISPResso in the previous ruleset, aggregates files, determines samples with sufficient modified reads and quantifies the effectiveness of the mutation (Frameshift, non-synonymous mutation)

BARCODES = config['demultiplex_barcode_tsv']
AMPLICONREF = config['AmpliconReference']
ANALYSIS = config['Analysis']

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


rule generateGraphs:
    input:
        crp = expand("results/3_crispresso/finished_crispresso_{sample}.touch", sample = config.get("samples").keys())
    params:
        script = srcdir("../scripts/generateGraphs.R")
    output:
        "results/4_quantifyMutation/finished_Graphs.touch"
    shell:
        """
	    Rscript {params.script}
        touch results/4_quantifyMutation/finished_Graphs.touch
        """


rule finished_quantification:
    input:
        crp = expand("results/3_crispresso/finished_crispresso_{sample}.touch", sample = config.get("samples").keys()),
        graphs = "results/4_quantifyMutation/finished_Graphs.touch"
    params:
        script = srcdir("../scripts/quantifyMutations.py"),
	    dir = srcdir("../results/3_crispresso/")
    output:
        "results/4_quantifyMutation/finished_quantification.touch"
    shell:
        """
	    python {params.script} -dir {params.dir}
        touch results/4_quantifyMutation/finished_quantification.touch
        """
