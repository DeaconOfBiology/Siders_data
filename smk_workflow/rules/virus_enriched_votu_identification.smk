rule merge_viral_reads:
    input:
        r1 = lambda wildcards: [config["clean_reads"] + f"vcontig_cell_enriched_reads/vContig_{wildcards.sample}-cell-{wildcards.group}_{fraction}_R1.fq.gz" for fraction in config["Samples"][wildcards.group]["cell"][wildcards.sample]],
        r2 = lambda wildcards: [config["clean_reads"] + f"vcontig_cell_enriched_reads/vContig_{wildcards.sample}-cell-{wildcards.group}_{fraction}_R2.fq.gz" for fraction in config["Samples"][wildcards.group]["cell"][wildcards.sample]],
        viral_retentate_r1=config["clean_reads"] + "merged_reads/merged_{sample}-virus-{group}_r1.fq.gz",
        viral_retentate_r2=config["clean_reads"] + "merged_reads/merged_{sample}-virus-{group}_r2.fq.gz"
    output:
        r1=config["clean_reads"] + "merged_reads/merged_{sample}-virus-{group}_vContig_r1.fq.gz",
        r2=config["clean_reads"] + "merged_reads/merged_{sample}-virus-{group}_vContig_r2.fq.gz"
    shell:
        """
            cat  {input.r1} {input.viral_retentate_r1}> {output.r1}
            cat  {input.r2} {input.viral_retentate_r2}> {output.r2}
        """
#metaspades_vContigs
rule megahit_vContigs:
    input:
        r1=config["clean_reads"] + "merged_reads/merged_{sample}-virus-{group}_vContig_r1.fq.gz",
        r2=config["clean_reads"] + "merged_reads/merged_{sample}-virus-{group}_vContig_r2.fq.gz"
    output:
        meta=config["assemblies"] + "{sample}-virus-{group}/megahit/final.contigs.fa",
    threads: 12
    # params:
    #     out=directory(config["assemblies"] + "{sample}-virus-{group}/megahit"),
    log:
        log=config["logs"] + "megahit/{sample}-vContig-{group}_megahit.log",
    shell:
        """ 
            (rm -rf data/processed/Assemblies/{wildcards.sample}-virus-{wildcards.group}/megahit/temp
            megahit -t {threads} -1 {input.r1} -2 {input.r2} -o  data/processed/Assemblies/{wildcards.sample}-virus-{wildcards.group}/megahit/temp
            mv data/processed/Assemblies/{wildcards.sample}-virus-{wildcards.group}/megahit/temp/* data/processed/Assemblies/{wildcards.sample}-virus-{wildcards.group}/megahit/
            ) 2> {log.log}
        """
            #(metaspades.py --threads {threads} -1 {input.r1} -2 {input.r2} -o {params.out}) 2> {log.log}



rule vContig_virsorter2_id:
    input:
        metasp=config["assemblies"] + "{sample}-virus-{group}/megahit/final.contigs.fa",
    output:
        folder = directory(config["assemblies"] + "{sample}-virus-{group}/virsorter_first/"),
        metasp = config["assemblies"] + "{sample}-virus-{group}/virsorter_first/final-viral-combined.fa"
    threads: 8
    log:
        meta=config["logs"] + "virsorter/{sample}-{group}-vContig_vs_run1.log",
    shell:
        """
            virsorter config --set HMMSEARCH_THREADS={threads}
            (virsorter run --keep-original-seq -i {input.metasp} -w {output.folder}  --min-length 5000 --min-score 0.5 --include-groups "dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae" -j {threads} all --scheduler greedy) 2> {log.meta}
        """

rule vContig_checkV:
    input:
        virsorter_meta = config["assemblies"] + "{sample}-virus-{group}/virsorter_first/final-viral-combined.fa"
    output:
        folder = directory(config["assemblies"] + "{sample}-virus-{group}/checkv/"),
        outfiles = [config["assemblies"] + "{sample}-virus-{group}/checkv/proviruses.fna",config["assemblies"] + "{sample}-virus-{group}/checkv/viruses.fna"]
    threads: 16
    log: 
        first=config["logs"] + "checkv/{sample}-{group}-vContig_fist_metaspades_checkv.log",
    shell:
        """
            (checkv end_to_end {input.virsorter_meta} {output.folder} -t {threads} -d  /users/troger50/checkv-db-v1.5) 2> {log.first}
        """

rule vContig_virsorter2_annotation_prep:
    input:
        v1=config["assemblies"] + "{sample}-virus-{group}/checkv/proviruses.fna",
        v2=config["assemblies"] + "{sample}-virus-{group}/checkv/viruses.fna"
    output:
        folder = directory(config["assemblies"] + "{sample}-virus-{group}/virsorter_final/"),
        fasta = config["assemblies"] + "{sample}-virus-{group}/virsorter_final/for-dramv/final-viral-combined-for-dramv.fa",
        tab = config["assemblies"] + "{sample}-virus-{group}/virsorter_final/for-dramv/viral-affi-contigs-for-dramv.tab"
    params:
        v=config["assemblies"] + "{sample}-virus-{group}/checkv/vContig-clean_viral_combined.fasta"
    log: 
        log=config["logs"] + "virsorter/{sample}-{group}-vContig_vs_run2.log",
    threads: 16
    shell:
        """
            cat {input.v1} {input.v2} > {params.v}

            virsorter config --set HMMSEARCH_THREADS={threads}
            (virsorter run --seqname-suffix-off --viral-gene-enrich-off --provirus-off --prep-for-dramv -i {params.v} -w {output.folder} --min-length 5000 --min-score 0.5 --include-groups "dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae" -j {threads} all --scheduler greedy) 2> {log.log}
        """

# rule vContig_quast:
#     input:
#         samps="data/processed/Viral_Assemblies/vContig_VirSorter2_pass2/virsorter/for-dramv/final-viral-combined-for-dramv.fa",
#     output:
#         "quality_checks/quast/vContig/report.tsv"
#     threads: 8
#     params:
#         out="quality_checks/quast/vContig/"
#     shell:
#         """
#             quast.py {input} --threads {threads} -o {params.out} --space-efficient -L 
#         """

rule vContig_clean_headers:
    input:
        meta=config["assemblies"] + "{sample}-virus-{group}/virsorter_final/for-dramv/final-viral-combined-for-dramv.fa",
    output:
        fa=config["assemblies"] + "{sample}-virus-{group}/anvio/clean_fasta_file/contigs_clean_headers.fa",
        report=config["assemblies"] + "{sample}-virus-{group}/anvio/clean_fasta_file/reformat-report.txt"
    params:
        v=lambda wildcards: sanitize_cell_name(f"{wildcards.sample}_virus_{wildcards.group}"),
    threads: 16
    log:
        log=config["logs"] + "anvio/vContigs_{sample}-{group}-fixed_fasta_headers.log",
    shell:
        """
            module load anvio/7
            (anvi-script-reformat-fasta {input.meta} -o {output.fa} -l 1000 --simplify-names --prefix {params.v} --report-file {output.report})2>{log.log}
        """