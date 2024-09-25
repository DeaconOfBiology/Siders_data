shell.prefix("""
# http://linuxcommand.org/wss0150.php
PROGNAME=$(basename $0)

function error_exit
{{
#   ----------------------------------------------------------------
#   Function for exit due to fatal program error
#       Accepts 1 argument:
#           string containing descriptive error message
#   ----------------------------------------------------------------
    echo "${{PROGNAME}}: ${{1:-"Unknown Error"}}" 1>&2
    exit 1
}}
""")

rule bowtie2_index:
    input:
        scaffold=config["assemblies"] + "{sample}-cell-{group}/virsorter_first/final-viral-combined.fa"
    output:
        scaffold1=config["assemblies"] + "{sample}-cell-{group}/cell_viral_sequences_scaffold.1.bt2"
    params:
        scaffold=config["assemblies"] + "{sample}-cell-{group}/cell_viral_sequences_scaffold"
    threads: 8
    log:
        scaffold=config["logs"] + "bowtie2/vcontig_{sample}-cell-{group}_sequences_bowtie2_index.log"
    shell:
        """
            (bowtie2-build --threads {threads} {input.scaffold} {params.scaffold}) 2> {log.scaffold}
        """


rule bowtie2_sam:
    input:
        scaffold1=config["assemblies"] + "{sample}-cell-{group}/cell_viral_sequences_scaffold.1.bt2",
    output:
        sam=temp(config["temp_files"] + "vContigs_{sample}-cell-{group}_{fraction}.sam")
    threads:16
    log:
        log=config["logs"] + "bowtie2/vContigs_{sample}-cell-{group}_{fraction}_samfile.log"
    params:
        scaffold=config["assemblies"] + "{sample}-cell-{group}/cell_viral_sequences_scaffold",
        r1=config["clean_reads"] + "clean_reads/clean_{sample}-cell-{group}_{fraction}_r1.fq.gz",
        r2=config["clean_reads"] + "clean_reads/clean_{sample}-cell-{group}_{fraction}_r2.fq.gz"
    shell:
        """
            (bowtie2 -p {threads} -x {params.scaffold} -1 {params.r1} -2 {params.r2} -S {output.sam})2>{log.log}
        """

rule samtools_view:
    input:
        sam=config["temp_files"] + "vContigs_{sample}-cell-{group}_{fraction}.sam"
    output:
        raw=temp(config["temp_files"] + "RAW_vContigs_{sample}-cell-{group}_{fraction}.bam")
    threads: 8
    shell:
        """
            samtools view -@ {threads} -F 2316 -bS {input.sam} -o {output.raw}
        """

rule samtools_sort:
    input:
        raw=config["temp_files"] + "RAW_vContigs_{sample}-cell-{group}_{fraction}.bam"
    output:
        bam=config["assemblies"] + "{sample}-cell-{group}/bam/vContigs_{sample}-cell-{group}_{fraction}.bam"
    threads: 8
    shell:
        """
            samtools sort -@ {threads} {input.raw} -o {output.bam}
        """

rule samtools_index: 
    input:
        bam=config["assemblies"] + "{sample}-cell-{group}/bam/vContigs_{sample}-cell-{group}_{fraction}.bam",
    output:
        bai=config["assemblies"] + "{sample}-cell-{group}/bam/vContigs_{sample}-cell-{group}_{fraction}.bam.bai",
    threads: 8
    shell:
        """
            samtools index -@ {threads} {input.bam}
        """

rule vcotig_reads:
    input:
        bam=config["assemblies"] + "{sample}-cell-{group}/bam/vContigs_{sample}-cell-{group}_{fraction}.bam.bai",
    output:
        r1=config["clean_reads"] + "vcontig_cell_enriched_reads/vContig_{sample}-cell-{group}_{fraction}_R1.fq.gz",
        r2=config["clean_reads"] + "vcontig_cell_enriched_reads/vContig_{sample}-cell-{group}_{fraction}_R2.fq.gz",
    threads: 8
    params:
        bam=config["assemblies"] + "{sample}-cell-{group}/bam/vContigs_{sample}-cell-{group}_{fraction}.bam",
    log:
        log=config["logs"] + "bowtie2/vContigs_{sample}-cell-{group}_{fraction}_reads.log",
    shell:
        """
            (samtools fastq -@ {threads}  {params.bam} -1 {output.r1} -2 {output.r2} -0 /dev/null -s /dev/null -n)2> {log.log}
        """