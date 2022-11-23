import argparse
import subprocess
import os.path

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
        "-r1",
        "--forward",
        dest = "forward",
        nargs="+",
        required=True,
        help="forward FQ"
    )
    requiredArgs.add_argument(
        "-r2",
        "--reverse",
        dest = "reverse",
        nargs="+",
        required=True,
        help="reverse FQ"
    )
    requiredArgs.add_argument(
        "-o",
        "--output",
        dest = "output",
        nargs="+",
        required=True,
        help="output Directory"
    )
    requiredArgs.add_argument(
        "-n",
        "--name",
        dest = "name",
        nargs="+",
        required=True,
        help="output name of the report"
    )
    requiredArgs.add_argument(
        "-t",
        "--AnalysisType",
        dest = "type",
        nargs="+",
        required=True,
        help="Run analysis on Forward fq or both FQ"
    )
    
    return parser.parse_args()

#Opens Alleles.tsv, takes sequences/names and returns strings in CRISPResso's preferred format x,y,z.
def getAmpliconSeqs(alleleFile,fqFile):

    basename = os.path.basename(fqFile)

    lines = [line.strip() for line in open(alleleFile).readlines()]
    for line in lines:
        split = line.split('\t')
        if split[0] in basename:
            #AlleleNames, AlleleSequence, GuideSequence, CodingSequence
            return [split[1],split[2],split[3],split[4]]
        else:
            continue
    return []
#plot_window_size 12 due to the edge of forwards
#CRISPResso -r1 {input.r1} -r2 {input.r2} --amplicon_seq {params.AR[1]} --amplicon_name {params.AR[0]} --guide_seq {params.AR[2]} --coding_seq {params.AR[3]} -o results/3_crispresso/ -n {wildcards.bc}_{wildcards.sample}_Paired -w 0 --exclude_bp_from_left 0 --exclude_bp_from_right 0
def runCRISPResso():

    args = parse_args()
    AlleleFile = args.Alleles[0]
    r1 = args.forward[0]
    r2 = args.reverse[0]
    outputDir = args.output[0]
    reportName = args.name[0]
    analysisType = args.type[0]
    #AlleleNames, AlleleSequence, GuideSequence, CodingSequence
    params = getAmpliconSeqs(AlleleFile,r1)

    if len(params) == 0:
        print("ERROR: did not find allele information")

    if analysisType == "Paired":
        subprocess.run(["CRISPResso","-r1",r1,"-r2",r2,"--amplicon_seq",params[1],"--amplicon_name",params[0],"--guide_seq",params[2],"--coding_seq",params[3],"-o","results/3_crispresso/","-n",reportName])
    elif analysisType == "Forward":
        subprocess.run(["CRISPResso","-r1",r1,"--amplicon_seq",params[1],"--amplicon_name",params[0],"--guide_seq",params[2],"--coding_seq",params[3],"-o","results/3_crispresso/","-n",reportName,"-w 2","--exclude_bp_from_left","0","--exclude_bp_from_right","0","--plot_window_size","12"])

if __name__ == "__main__":
	runCRISPResso()