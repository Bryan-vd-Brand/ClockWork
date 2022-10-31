configfile: "config.yaml"
report: "report/workflow.rst"
include: "rules/0_fastqc.smk"
include: "rules/1_DemultiplexTrim.smk"
include: "rules/2_fastqc.smk"
include: "rules/3_crispresso.smk"

print(config.get("samples"))

def get_samples(wildcards):
    print(config["samples"][wildcards.sample])
    return config["samples"][wildcards.sample].split(",") # "," split


print(config.get('input_dir'))
DATA_DIR = config.get('input_dir')

rule all:
    input:
        "results/0_fastqc/multiqc_report.html",
        "results/1_DemultiplexTrim/DidDemux.touch",
        "results/2_fastqc/multiqc_report.html",
        expand("results/3_crispresso/finished_crispresso_{sample}.touch", sample = config.get("samples").keys())