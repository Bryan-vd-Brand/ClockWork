#This rule takes in a directory structure created by CRISPResso in the previous ruleset, aggregates files, determines samples with sufficient modified reads and quantifies the effectiveness of the mutation (Frameshift, non-synonymous mutation)

BARCODES = config['demultiplex_barcode_tsv']
AMPLICONREF = config['AmpliconReference']
ANALYSIS = config['Analysis']
ALLELES = config['Alleles']

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

rule generateKnockOutReport:
    input:
        crp = expand("results/3_crispresso/finished_crispresso_{sample}.touch", sample = config.get("samples").keys())
    params:
        script = srcdir("../scripts/quantifyMutations.py"),
	    dir = srcdir("../results/3_crispresso/")
    output:
        "results/4_quantifyMutation/finished_knockout.touch"
    shell:
        """
	    python {params.script} -dir {params.dir}
        touch results/4_quantifyMutation/finished_knockout.touch
        """

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

rule generate96Plate:
    input:
        crp = expand("results/3_crispresso/finished_crispresso_{sample}.touch", sample = config.get("samples").keys())
    params:
        script = srcdir("../scripts/generatePlate.R")
    output:
        "results/4_quantifyMutation/finished_Plate.touch"
    shell:
        """
	    Rscript {params.script}
        touch results/4_quantifyMutation/finished_Plate.touch
        """

rule generateIGVscreens:
    input:
        crp = expand("results/3_crispresso/finished_crispresso_{sample}.touch", sample = config.get("samples").keys())
    params:
        script = srcdir("../scripts/generateIGV.py"),
	    dir = srcdir("../results/3_crispresso/"),
        outDir = srcdir("../results/4_quantifyMutation/")
    output:
        "results/4_quantifyMutation/finished_IGV.touch"
    shell:
        """
	    python {params.script} -dir {params.dir} -a {ALLELES} -outDir {params.outDir} -ref {AMPLICONREF}
        touch results/4_quantifyMutation/finished_IGV.touch
        """

rule generatePDFs:
    input:
        crp = expand("results/3_crispresso/finished_crispresso_{sample}.touch", sample = config.get("samples").keys()),
        wellpdf = "results/4_quantifyMutation/finished_Plate.touch"
    params:
        script = srcdir("../scripts/generatePDF.py"),
	    dir = srcdir("../results/3_crispresso/"),
        outDir = srcdir("../results/4_quantifyMutation/"),
        KOReport = srcdir("../KnockoutReport.tsv"),
        wellPDF = srcdir("../wellPlot.pdf")
    output:
        "results/4_quantifyMutation/finished_generatePDFs.touch"
    shell:
        """
	    python {params.script} -dir {params.dir} -outDir {params.outDir} -KO {params.KOReport}
        rm *_name.pdf
        touch results/4_quantifyMutation/finished_generatePDFs.touch
        """

rule finished_quantification:
    input:
        ko = "results/4_quantifyMutation/finished_knockout.touch",
        graphs = "results/4_quantifyMutation/finished_Graphs.touch",
        plate = "results/4_quantifyMutation/finished_Plate.touch",
        IGV = "results/4_quantifyMutation/finished_IGV.touch",
        PDF = "results/4_quantifyMutation/finished_generatePDFs.touch"
    output:
        "results/4_quantifyMutation/finished_quantification.touch"
    shell:
        """
        touch results/4_quantifyMutation/finished_quantification.touch
        """
