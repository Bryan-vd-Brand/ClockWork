#This rule file performs fastqc/multiqc on the original data set.

DATA_DIR = config['input_dir']

rule fastqc_paired:
    input:
        r1 = f"{DATA_DIR}/{{sample}}/{{sample}}_L4_1.fq.gz",
        r2 = f"{DATA_DIR}/{{sample}}/{{sample}}_L4_2.fq.gz"
    output:
        expand("results/0_fastqc/{{sample}}/{{sample}}_L4_{read}_fastqc.zip", read=["1", "2"]),
        expand("results/0_fastqc/{{sample}}/{{sample}}_L4_{read}_fastqc.html", read=["1", "2"]),
        outdir = directory("results/0_fastqc/{sample}/")
    threads: 2
    log:
        stdout = "logs/0_fastqc/{sample}/{sample}_fastqc.stdout",
        stderr = "logs/0_fastqc/{sample}/{sample}_fastqc.stderr"
    shell:
        """
        fastqc -o {output.outdir} -t {threads} {input} 1> {log.stdout} 2> {log.stderr}
        """

rule multiqc:
    input:
        expand("results/0_fastqc/{sample}/", sample = config.get("samples").keys())
    output:
        "results/0_fastqc/multiqc_report.html"
    params:
        outdir = "results/0_fastqc"
    shell:
        """
        multiqc -v -f {input} -o {params.outdir}
        """