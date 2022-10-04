#This rule file performs demultiplexing on paired fastq files, xxx_L4_1_fq.gz

DATA_DIR = config['input_dir']

rule demultiplexing_paired:
    input:
        r1 = f"{DATA_DIR}/{{sample}}/{{sample}}_L4_1.fq.gz",
        r2 = f"{DATA_DIR}/{{sample}}/{{sample}}_L4_2.fq.gz"
    output:
    
    threads: 2
    log:
        stdout = "logs/1_demultiplex/{sample}/{sample}_fastqc.stdout",
        stderr = "logs/1_demultiplex/{sample}/{sample}_fastqc.stderr"
    shell:
        """
        fastqc -o {output.outdir} -t {threads} {input} 1> {log.stdout} 2> {log.stderr}
        """
