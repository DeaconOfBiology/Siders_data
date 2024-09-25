#################################################
# This is the first snake file that needs to be #
# generated before the analysis can move on.    #
# In it, we preproses the samples by:           #
# 1. evaluating the quality of the raw reads    #
# 2. trimming the raw reads                     #
# 3, evaluting the trimmed (clean) reads        #
#################################################

configfile:"/projects/luo_lab/Siders_data/smk_workflow/config/configure.yml"

rule all_run:
    input:
        multiqc=config["quality_checks"] + "multiqc/multiqc_report.html",
        

rule multiqc:
    input:
        vr1 = [expand([config["quality_checks"] + "fastqc/merged_{sample}-{organism}-{group}_r1_fastqc.html"],group=g,organism=o,sample=s)
                for g, org in config["Samples"].items()
                for o, samp in org.items()
                for s in samp.keys()],        
        vr2=  [expand([config["quality_checks"] + "fastqc/merged_{sample}-{organism}-{group}_r2_fastqc.html"],group=g,organism=o,sample=s)
                for g, org in config["Samples"].items()
                for o, samp in org.items()
                for s in samp.keys()],
        scaffolds=config["quality_checks"] + "quast/report.tsv"
    output:
        config["quality_checks"] + "multiqc/multiqc_report.html"
    params:
        in_put=config["quality_checks"],
        out_put=config["quality_checks"] + "multiqc/"
    shell:
        """
           multiqc -d -dd 1 {params.in_put} -o {params.out_put} --export
        """

rule fastqc:
    output:
        merged=[config["quality_checks"] + "fastqc/merged_{sample}-{organism}-{group}_r1_fastqc.html", config["quality_checks"] + "fastqc/merged_{sample}-{organism}-{group}_r2_fastqc.html"]
    input:  
        merged=[config["clean_reads"] + "merged_reads/merged_{sample}-{organism}-{group}_r1.fq.gz",config["clean_reads"] + "merged_reads/merged_{sample}-{organism}-{group}_r2.fq.gz"],
    params:
        config["quality_checks"] + "fastqc/"
    log:
        merged=config["logs"] + "fastqc/merged_{sample}-{organism}-{group}_fastqc_clean.log",
    shell:
        """
           (fastqc -o {params} {input.merged}) 2> {log.merged}
        """

rule quast:
    input:
        cell_control=[expand([config["assemblies"] + "{sample}-{organism}-{group}/anvio/clean_fasta_file/contigs_clean_headers.fa"],group=g,organism=o,sample=s)
                for g, org in config["Samples"].items()
                for o, samp in org.items()
                for s in samp.keys()]
    output:
        config["quality_checks"] + "quast/report.tsv"
    threads: 8
    params:
        out="quality_checks/quast/"
    shell:
        """
            quast.py {input} --threads {threads} -o {params.out} --space-efficient -L 
        """

#Trim addapters from raw reads
rule bbduk_adp:
    input: 
        r1=config["raw_reads"] + "{sample}-{organism}-{group}_{fraction}_r1.fq.gz",
        r2=config["raw_reads"] + "{sample}-{organism}-{group}_{fraction}_r2.fq.gz"

    output:
        r1=temp(config["temp_files"] + "first_pass_{sample}-{organism}-{group}_{fraction}_r1.fq.gz"),
        r2=temp(config["temp_files"] + "first_pass_{sample}-{organism}-{group}_{fraction}_r2.fq.gz")
    resources:
        mem_mb = 25000
    log:
        config["logs"] + "bbduk/{sample}-{organism}-{group}_{fraction}_bbduk_adp.log"
    params:
        bbduk=config["bbduk"] + "bbduk.sh",
        adapters=config["adapters"] + "adapters.fa"
    shell:
        """
            ({params.bbduk} -Xmx1g in1={input.r1} in2={input.r2} \
                out1={output.r1} \
                out2={output.r2} \
                ref={params.adapters} ktrim=r k=21 qin=auto mink=11 hdist=2 tbo tpe) 2> {log}
        """

# Further trimming step to increase the quality
rule bbduk_phix:
    input: 
        r1=config["temp_files"] + "first_pass_{sample}-{organism}-{group}_{fraction}_r1.fq.gz",
        r2=config["temp_files"] + "first_pass_{sample}-{organism}-{group}_{fraction}_r2.fq.gz"

    output:
        r1=config["clean_reads"] + "clean_reads/clean_{sample}-{organism}-{group}_{fraction}_r1.fq.gz",
        r2=config["clean_reads"] + "clean_reads/clean_{sample}-{organism}-{group}_{fraction}_r2.fq.gz"
    params:
        bbduk=config["bbduk"] + "bbduk.sh",
        ref=config["reference"] + "phix174_ill.ref.fa.gz"
    resources:
        mem_mb = 25000
    log:
        config["logs"] + "bbduk/{sample}-{organism}-{group}_{fraction}_bbduk_qal.log"
    shell:
        """
            ({params.bbduk} -Xmx1g in1={input.r1} in2={input.r2} \
                out1={output.r1} \
                out2={output.r2} \
                ref={params.ref} k=27 hdist=1 qtrim=rl trimq=17 cardinality=t mingc=0.05 maxgc=0.95) 2> {log}
        """

rule merge_clean_reads:
    input:
        r1=lambda wildcards: [config["clean_reads"] + f"clean_reads/clean_{wildcards.sample}-{wildcards.organism}-{wildcards.group}_{fraction}_r1.fq.gz" for fraction in config["Samples"][wildcards.group][wildcards.organism][wildcards.sample]],
        r2=lambda wildcards: [config["clean_reads"] + f"clean_reads/clean_{wildcards.sample}-{wildcards.organism}-{wildcards.group}_{fraction}_r2.fq.gz" for fraction in config["Samples"][wildcards.group][wildcards.organism][wildcards.sample]]
    output:
        r1=config["clean_reads"] + "merged_reads/merged_{sample}-{organism}-{group}_r1.fq.gz",
        r2=config["clean_reads"] + "merged_reads/merged_{sample}-{organism}-{group}_r2.fq.gz"
    shell:
        """
            cat  {input.r1} > {output.r1}
            cat  {input.r2} > {output.r2}
        """