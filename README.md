# ClockWork Usage Guide

1. Clone the git at https://github.com/Bryan-vd-Brand/ClockWork
2. Open command line in the cloned directory
2a. Ensure bioconda is available as a channel for anaconda
2b. Create the conda environment with the required libraries ; conda create --name YOURNAME --file ClockWorkEnvSpecFile.txt
3. Load the env; Test for snakemake by running "snakemake -c1" -> expected result is "Building DAG of jobs ..." into WorkflowError:
4a. Add your fastq files to the data directory (Or change the input_dir variable to your data location)
4b. Ensure your fastq files are placed in a directory with the same name as the sample, a.k.a X/data/SampleName/SampleName.fq
5. Change the path of STRELKA to your strelka installation
6. To generate the list of samples in the config file Run: python ./scripts/GenerateConfigFile.py -c config.yaml  
7. Verify that all paths are correct
8. Replace Barcodes.tsv with your barcodes for this sample, use the format described in config.yaml
9. Replace FASTA format reference sequence
10. Edit Alleles.tsv and describe your "Samplename \t AlleleNames \t AlleleSequence \t GuideSequence \t CodingSequence" ; Seperate by tab, specify multiple alleles by seperating them by comma ; Example for CodingSequence: ATAT,ATGT
11. Edit Analysis.tsv (NAME \t PAIRED) and describe the usage for each dataset, where the name is a subsequence of the sampleset's name
12. run the pipeline by snakemake -c{NUM OF CORES} {SNAKEFILE'.}
13. If you run into issues with wildcard matching (happens if the fq files contain extra _'s) Comment (#) the last 3 lines of rule all and run the pipeline for the first 3 steps, then uncomment and run the last 3 steps