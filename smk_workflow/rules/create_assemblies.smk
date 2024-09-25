# ##################################################
# # This is the second snake file that needs to be #
# # generated before the analysis can move on. In  #
# # it, we create and assess the assemblies by:    #
# # 1. running metaspades                          #
# # 2. running metaviralspades                     #
# # 3. Runing QUANT for assembly stats             #
# # 4. Running multiqc                             #
# ##################################################
#metaspades_1st_run
rule megahit_1st_run:
    input:
        r1=config["clean_reads"] + "merged_reads/merged_{sample}-cell-{group}_r1.fq.gz",
        r2=config["clean_reads"] + "merged_reads/merged_{sample}-cell-{group}_r2.fq.gz"
    output:
        meta=config["assemblies"] + "{sample}-cell-{group}/megahit/final.contigs.fa"
    threads: 12
    # params:
    #     out=directory(config["assemblies"] + "{sample}-cell-{group}/megahit/")
    log:
        cell=config["logs"] + "megahit/{sample}-cell-{group}_megahit.log"
    shell:
        """
            (rm -rf data/processed/Assemblies/{wildcards.sample}-cell-{wildcards.group}/megahit/temp
            megahit -t {threads} -1 {input.r1} -2 {input.r2} -o data/processed/Assemblies/{wildcards.sample}-cell-{wildcards.group}/megahit/temp
            mv data/processed/Assemblies/{wildcards.sample}-cell-{wildcards.group}/megahit/temp/* data/processed/Assemblies/{wildcards.sample}-cell-{wildcards.group}/megahit/) 2> {log.cell}
        """
            #(metaspades.py --threads {threads} -1 {input.r1} -2 {input.r2} -o {params.out}) 2> {log.cell}

# Function to sanitize the cell name
def sanitize_cell_name(cell_name):
    return re.sub(r'[^a-zA-Z0-9]', '_', cell_name)

rule clean_headers:
    input:
        meta=config["assemblies"] + "{sample}-cell-{group}/megahit/final.contigs.fa"
    output:
        fa=config["assemblies"] + "{sample}-cell-{group}/anvio/clean_fasta_file/contigs_clean_headers.fa",
        report=config["assemblies"] + "{sample}-cell-{group}/anvio/clean_fasta_file/scaffolds-clean_headers.txt"
    params:
        cell=lambda wildcards: sanitize_cell_name(f"{wildcards.sample}_cell_{wildcards.group}")
    threads: 16
    log:
        names=config["logs"] + "anvio/{sample}-cell-{group}_fixed_fasta_headers.log"
    shell:
        """
            module load anvio/7
            (anvi-script-reformat-fasta {input.meta} -o {output.fa} -l 1000 --simplify-names --prefix {params.cell} --report-file {output.report})2>{log.names}
        """
