import argparse
from importlib.util import module_for_loader
import os
import glob
from random import sample 
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser()

    requiredArgs = parser.add_argument_group("Required arguments")

    requiredArgs.add_argument(
        "-dir",
        "--directory",
        dest = "directory",
        nargs="+",
        required=True,
        help="directory containing folders with CRISPRessoPooled results"
    )
    
    return parser.parse_args()

#Script crawls directory of CRISPRessoPooled results, 
# retrieves % (un-)modified and #reads from SAMPLES_QUANTIFICATION_SUMMARY.txt -> determines sucess of CRISPR-cas9
# If sufficient modified reads, crawls directory for CRISPResso_on_xxx directories, retrieves the Frameshift_analysis.txt and Alleles_frequency_table_around_*.txt
# Uses the table provided to determine the size of indel and call the mutation as a frameshift or not (# of indel bases % 3 != 0)
def getMutations():
    args = parse_args()
    sourceDir = args.directory[0]
    print(F"Opening {sourceDir} and quantifiying mutation")

    if not os.path.exists(sourceDir):
        print("ERROR: directory supplied for quantification does not exist.")

    ResultDataFrame = pd.DataFrame(columns=["SampleName","Reference","Unmodified","Reads_aligned","Knockout","Reason","Insertions","Deletions","Substitutions"])

    print("glob")
    for directory in glob.iglob(sourceDir + "CRISPResso_on_*"):

        #Dirty Hack for bad iglobs
        if ".html" in directory:
            continue

        lastDir = os.path.basename(os.path.normpath(directory))
        sampleName = lastDir.split('_')[2] + "_" + lastDir.split('_')[3]
        sampleSummaryFile = os.path.join(directory,"CRISPResso_quantification_of_editing_frequency.txt")
        sampleSummary = pd.read_table(sampleSummaryFile, sep='\t',header = 0)
        amplicons = sampleSummary['Amplicon']
        unmodified = sampleSummary['Unmodified%']
        readsAligned = sampleSummary['Reads_aligned']

        for i in range(0,len(amplicons)):
            if unmodified[i] < 5.0 and readsAligned[i] > 1000:
                #95% of data is modified reads, and atleast 1K total reads were used
                #Iterate subfolders and determine indels
                
                frameshiftFile = directory + f"/{amplicons[i]}.Frameshift_analysis.txt"
                if not os.path.exists(frameshiftFile):
                    print(f"ERROR: missing frameshiftFile: {frameshiftFile}")
                
                alleleFreqFile = directory + f"/{amplicons[i]}.Alleles_frequency_table_around_*.txt"
                for aff in glob.iglob(alleleFreqFile):
                    alleleFreqFile = aff                               
                if not os.path.exists(alleleFreqFile):
                    print(f"ERROR: missing alleleFreqFile: {alleleFreqFile}")

                print(F"Processing {alleleFreqFile}")    
                
                #Open up allele freq, determine # of indel bases
                alleleFreqTable = pd.read_table(alleleFreqFile, sep='\t', header = 0)
                #Header = Aligned_Sequence	Reference_Sequence	Unedited	n_deleted	n_inserted	n_mutated	#Reads	%Reads
                NoAlleles = 0
                for x in range(0,len(alleleFreqTable)):
                    if alleleFreqTable.iloc[x,7] > 5.0:
                        NoAlleles += 1
                    else:
                        break
                
                #variables across all alleles
                totalModifiedReads = 0
                allEdited = True
                deleted = ""
                inserted = ""
                mutated = ""

                for y in range(0,NoAlleles):
                    Aligned_sequence = alleleFreqTable.iloc[y,0]
                    Unedited = alleleFreqTable.iloc[y,2]
                    deleted = deleted + str(alleleFreqTable.iloc[y,3]) + ","
                    inserted = inserted + str(alleleFreqTable.iloc[y,4]) + ","
                    mutated = mutated + str(alleleFreqTable.iloc[y,5]) + ","
                    n_reads = alleleFreqTable.iloc[y,6]
                    percentage_reads = alleleFreqTable.iloc[y,7]

                    #One of the alleles is functional?
                    if Unedited == False:
                        totalModifiedReads += n_reads
                    else:
                        allEdited = False

                if allEdited:
                    #All alleles are modified, check for frameshift or in-frame
                    file = open(frameshiftFile, 'r')
                    lines = file.readlines()
                    noncoding = lines[1].split(':')[1].split(' ')[0]
                    inFrame = lines[2].split(':')[1].split(' ')[0]
                    frameshift = lines[3].split(':')[1].split(' ')[0]

                    #Inframe - Frameshift ; if atleast 10% of transcripts for each alelles (sub)-variants are inFrame we assume the phenotype to be fit
                    if int(inFrame) > int(frameshift) * 0.1:
                        #report inframe
                        rowDF = pd.DataFrame(data={'SampleName':[sampleName],'Reference':[amplicons[i]],'Unmodified':[unmodified[i]],'Reads_aligned':[readsAligned[i]],'Knockout':["NO"],'Reason':["InFrame"],'Insertions':[inserted],'Deletions':[deleted],'Substitutions':[mutated]})
                        ResultDataFrame = pd.concat([ResultDataFrame,rowDF],sort=False)
                    else:
                    #report knockout
                        rowDF = pd.DataFrame(data={'SampleName':[sampleName],'Reference':[amplicons[i]],'Unmodified':[unmodified[i]],'Reads_aligned':[readsAligned[i]],'Knockout':["YES"],'Reason':["Frameshift"],'Insertions':[inserted],'Deletions':[deleted],'Substitutions':[mutated]})
                        ResultDataFrame = pd.concat([ResultDataFrame,rowDF],sort=False)

                else:

                    #One of the alleles is unEdited
                    rowDF = pd.DataFrame(data={'SampleName':[sampleName],'Reference':[amplicons[i]],'Unmodified':[unmodified[i]],'Reads_aligned':[readsAligned[i]],'Knockout':["NO"],'Reason':["FunctionalAllele"],'Insertions':[inserted],'Deletions':[deleted],'Substitutions':[mutated]})
                    ResultDataFrame = pd.concat([ResultDataFrame,rowDF],sort=False)

            #Failed to pass 1K / 5%
            else:
                if unmodified[i] > 5.0:
                    rowDF = pd.DataFrame(data={'SampleName':[sampleName],'Reference':[amplicons[i]],'Unmodified':[unmodified[i]],'Reads_aligned':[readsAligned[i]],'Knockout':["NO"],'Reason':["Unmodified"],'Insertions':["NA"],'Deletions':["NA"],'Substitutions':["NA"]})
                    ResultDataFrame = pd.concat([ResultDataFrame,rowDF],sort=False)
                    continue
                #less then 1K aligned reads
                rowDF = pd.DataFrame(data={'SampleName':[sampleName],'Reference':[amplicons[i]],'Unmodified':[unmodified[i]],'Reads_aligned':[readsAligned[i]],'Knockout':["NO"],'Reason':["InsufficientData"],'Insertions':["NA"],'Deletions':["NA"],'Substitutions':["NA"]})
                ResultDataFrame = pd.concat([ResultDataFrame,rowDF],sort=False)
                          
        
    ResultDataFrame.to_csv("KnockoutReport.tsv",sep='\t', index = False)
    print(ResultDataFrame)



if __name__ == "__main__":
	getMutations()
