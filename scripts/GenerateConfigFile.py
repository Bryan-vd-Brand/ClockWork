#This script opens the config.yaml, reads out the containing lines until #START
#Then takes the data directory, gathers all of the different fastq and adds them to the config file. #TODO: Support merge in fq (2 forward, 2 reverse fq, same sequencing experiment)
#Requires the following structure DataFile/SampleName/SampleName.fastq.gz (or paired ; SampleName_1.fastq.gz)
#  HsREV7: data/HsREV7/HsREV7_FKDN220301829-1A_HHFGYDSX3_L4.fq.gz
#  HsREV7: data/HsREV7/HsREV7_FKDN220301829-1A_HHFGYDSX3_L4_1.fq.gz, data/HsREV7/HsREV7_FKDN220301829-1A_HHFGYDSX3_L4_2.fq.gz

import glob
import os
import argparse

ConfigFile = "config.yaml"
DataFolder = ""

def parse_args():
    parser = argparse.ArgumentParser()

    requiredArgs = parser.add_argument_group("Required arguments")

    requiredArgs.add_argument(
        "-c",
        "--config_yaml",
        dest = "config_file",
        nargs="+",
        required=True,
        help="Config file for snakemake"
    )
    
    return parser.parse_args()


def main():
    args = parse_args()
    ConfigFile = args.config_file[0]
    with open(F'{ConfigFile}','rt') as configFile:
        line = configFile.readline()
        while("input_dir:" not in line):
            line = configFile.readline()
        DataFolder = configFile.readline().strip()
    with open("./temp_config.yaml",'w') as newConfigFile:
        with open(F'{ConfigFile}','rt') as configFile:
            line = configFile.readline()
            while("#START" not in line):
                newConfigFile.write(line)
                line = configFile.readline()
            newConfigFile.write(line) #START line, add it again
            #generate new entry for all samples in DataFolder
            newConfigFile.write("samples:\n")
            for sampleFolder in os.listdir(DataFolder):
                forward = glob.glob(F"{DataFolder}/{sampleFolder}/{sampleFolder}_*_1.fq.gz")
                reverse = glob.glob(F"{DataFolder}/{sampleFolder}/{sampleFolder}_*_2.fq.gz")
                #save with the entire file name for upstream snakemake input DAG
                if not "_L4_1.fq.gz" in str(os.path.basename(forward[0])):
                    print("ERROR: Failed to strip fq extension, L2/L4?")
                    exit
                strippedForward = str(os.path.basename(forward[0])).replace("_L4_1.fq.gz","")
                #Check if fq is accessable/permissions good
                if os.path.exists(forward[0]): #paired
                    newConfigFile.write(F"  {strippedForward}: [{forward[0]},{reverse[0]}]\n")
                elif os.path.exists(glob.glob(F"{DataFolder}/{sampleFolder}/*.fq.gz")): #unpaired
                    single = glob.glob(F"{DataFolder}/{sampleFolder}/*.fq.gz")
                    newConfigFile.write(F"  {sampleFolder}: {single}\n")
                else:
                    print(F"ERROR failed to create config entry for {sampleFolder}")
    os.system(F"rm -f {ConfigFile}")
    os.rename("temp_config.yaml","config.yaml")
    print(DataFolder)
    print(F"Finished generating new config file out of ^")

if __name__ == "__main__":
	main()
