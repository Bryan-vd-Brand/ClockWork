import argparse
from importlib.util import module_for_loader
import os
import glob 
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
# If sufficient modified reads, crawls directory for CRISPResso_on_xxx directories, retrieves the Quantification_window_modification_count_vectors.txt file (tsv)
# Uses the table provided to determine the size of indel and call the mutation as a frameshift or not (# of indel bases % 3 != 0)
def getMutations():
    args = parse_args()
    pooledDirectory = args.directory[0]
    print(F"Opening {pooledDirectory} and quantifiying mutation")

    if not os.path.exists(pooledDirectory):
        print("ERROR: directory supplied for quantification does not exist.")

    ResultDataFrame = pd.DataFrame(columns=["SampleName","Reference","Unmodified","Reads_aligned","Knockout","Insertions","Deletions","Substitutions"])

    print("glob")
    for directory in glob.iglob(pooledDirectory + "CRISPRessoPooled_on_*_Paired"):
        lastDir = os.path.basename(os.path.normpath(directory))
        sampleName = lastDir.split('_')[2] + "_" + lastDir.split('_')[3]
        sampleSummaryFile = os.path.join(directory,"SAMPLES_QUANTIFICATION_SUMMARY.txt")
        sampleSummary = pd.read_table(sampleSummaryFile, sep='\t',header = 0)
        amplicons = sampleSummary['Name']
        unmodified = sampleSummary['Unmodified%']
        readsTotal = sampleSummary['Reads_total']
        for i in range(0,len(amplicons)):
            if unmodified[i] < 5.0 and readsTotal[i] > 1000:
                #95% of data is modified reads, and atleast 1K total reads were used
                #Iterate subfolders and determine indels
                if not os.path.exists(directory + "/CRISPResso_on_" + amplicons[i]):
                    print("ERROR: missing expected CRISPResso folder")
                quantificationFile = directory + "/CRISPResso_on_" + amplicons[i] + "/Quantification_window_modification_count_vectors.txt"
                if not os.path.exists(quantificationFile):
                    print(f"ERROR: missing quantification file: {quantificationFile}")
                #Open up quantification, determine # of indel bases
                quantificationTable = pd.read_table(quantificationFile, sep='\t', header = 0)
                #each column is a base around the window, iterate columns. Header = G/T/A/C, 1st row = Insertions, 2rd row = Deletions, 3th row = Substitutions, 4th row = All_modifications, 5th row = Total
                Insertions = 0
                Deletions = 0
                Substitutions = 0
                for column in quantificationTable.columns[1:]:
                    values = quantificationTable[column].values
                    insertionCount = values[0]
                    deletionCount = values[1]
                    substitutionCount = values[2]
                    allModifications = values[3]
                    totalCount = values[4]
                    #if less then 1000 modified reads support the indel or sub, skip.
                    if allModifications < 1000:
                        continue
                    #Insertion
                    if insertionCount > 0.50 * totalCount:
                        Insertions+=1
                    #Deletion
                    if deletionCount > 0.50 * totalCount:
                        Deletions+=1
                    #Sub
                    if substitutionCount > 0.50 * totalCount:
                        Substitutions+=1
                #only interested in frameshifts, so %3
                moduloIndel = (Insertions+Deletions)%3

                if not moduloIndel == 0:
                    rowDF = pd.DataFrame(data={'SampleName':[sampleName],'Reference':[amplicons[i]],'Unmodified':[unmodified[i]],'Reads_aligned':[readsTotal[i]],'Knockout':["YES"],'Insertions':[Insertions],'Deletions':[Deletions],'Substitutions':[Substitutions]})
                    ResultDataFrame = pd.concat([ResultDataFrame,rowDF],sort=False)
                else:
                    rowDF = pd.DataFrame(data={'SampleName':[sampleName],'Reference':[amplicons[i]],'Unmodified':[unmodified[i]],'Reads_aligned':[readsTotal[i]],'Knockout':["NO"],'Insertions':[Insertions],'Deletions':[Deletions],'Substitutions':[Substitutions]})
                    ResultDataFrame = pd.concat([ResultDataFrame,rowDF],sort=False)
            #Failed to pass 1K / 5%
            else:
                rowDF = pd.DataFrame(data={'SampleName':[sampleName],'Reference':[amplicons[i]],'Unmodified':[unmodified[i]],'Reads_aligned':[readsTotal[i]],'Knockout':["NO"],'Insertions':["NA"],'Deletions':["NA"],'Substitutions':["NA"]})
                ResultDataFrame = pd.concat([ResultDataFrame,rowDF],sort=False)
                          
        
    ResultDataFrame.to_csv("KnockoutReport.tsv",sep='\t', index = False)
    print(ResultDataFrame)



if __name__ == "__main__":
	getMutations()
