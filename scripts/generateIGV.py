import argparse
import os.path
import pandas as pd
import glob

def parse_args():
    parser = argparse.ArgumentParser()

    requiredArgs = parser.add_argument_group("Required arguments")


    requiredArgs.add_argument(
        "-a",
        "--Alleles",
        dest = "Alleles",
        nargs="+",
        required=True,
        help="TSV specifying alleles"
    )

    requiredArgs.add_argument(
        "-ref",
        "--Reference",
        dest = "Reference",
        nargs="+",
        required=True,
        help="Fasta file specifying references"
    )

    requiredArgs.add_argument(
        "-dir",
        "--directory",
        dest = "directory",
        nargs="+",
        required=True,
        help="directory containing bam files"
    )

    requiredArgs.add_argument(
        "-outDir",
        "--outDir",
        dest = "outDir",
        nargs="+",
        required=True,
        help="place to store screenshots"
    )

    return parser.parse_args()

#Opens Alleles.tsv, returns samplenames
def getSampleNames(alleleFile):
    
    if os.path.exists(alleleFile):
        aDF = pd.read_table(alleleFile, sep="\t", header=None)
        return aDF.iloc[:,0]
    else:
        print("ERROR: did not find alleles.tsv")

#Loops a dir for a certain samplename
def loopBamfiles(samplename, directory, ref, outDir):
    batchedCommands = f"{samplename}_IGVBatch.txt"
    outDir = outDir + f"/{samplename}/"

    for bamfile in glob.iglob(directory + f"sorted_*_{samplename}_*.bam"):
        baseName = os.path.basename(os.path.normpath(bamfile))
        saveName = baseName.split('_')[1] + "_" + baseName.split('_')[2]
        vcfFile = glob.glob(directory + f"Strelka/{saveName}_*/results/variants/variants.vcf.gz")

        if len(vcfFile) > 1 or len(vcfFile) == 0:
            print("ERROR: Unexpected # of vcfFiles from Strelka") 
        
        vcfFile = vcfFile[0]

        with open (batchedCommands, 'a') as igv:
            igv.write("new\n")
            igv.write(f"genome {ref}\n")
            igv.write(f"load {bamfile}\n")
            igv.write(f"load {vcfFile}\n")
            igv.write(f"snapshotDirectory {outDir}\n")
            igv.write(f"goto {samplename}\n")
            igv.write("collapse\n")
            igv.write(f"snapshot {saveName}.png\n")

#method generates a text file runnable by IGV: https://software.broadinstitute.org/software/igv/batch https://github.com/igvteam/igv/wiki/Batch-commands
def generateIGV():

    args = parse_args()
    AlleleFile = args.Alleles[0]
    bamDir = args.directory[0]
    refFasta = args.Reference[0]
    outDir = os.path.abspath(args.outDir[0])
    Alleles = getSampleNames(AlleleFile)

    if not os.path.exists(bamDir):
        print("ERROR: directory supplied does not exist.")

    for allele in Alleles:
        loopBamfiles(allele,bamDir,refFasta,outDir)
   
if __name__ == "__main__":
	generateIGV()