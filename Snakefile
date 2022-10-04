configfile: "config.yaml"
report: "report/workflow.rst"
include: "rules/0_fastqc.smk"
include: "rules/1_demultiplex.smk"

# both methods return the same information. Notice how for the paired end samples,
# there is only one string containing the two files. Because of that we have to
# split by "," later to get the two files
print(config["samples"])
print(config.get("samples"))

def get_samples(wildcards):
    print(config["samples"][wildcards.sample])
    return config["samples"][wildcards.sample].split(",") # "," split


print(config.get('input_dir'))
DATA_DIR = config.get('input_dir')


rule all:
    input:
        "results/0_fastqc/multiqc_report.html"