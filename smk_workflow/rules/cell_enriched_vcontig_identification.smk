# Verify viral contigs created in metaviralspades
rule virsorter2_id:
    input:
        metasp=config["assemblies"] + "{sample}-cell-{group}/anvio/clean_fasta_file/contigs_clean_headers.fa"
    output:
        folder=directory(config["assemblies"] + "{sample}-cell-{group}/virsorter_first/"),
        emetasp=config["assemblies"] + "{sample}-cell-{group}/virsorter_first/final-viral-combined.fa"
    threads: 8
    log:
        meta=config["logs"] + "virsorter/{sample}-cell-{group}_vcontigs.log"
    shell:
        """
            virsorter config --set HMMSEARCH_THREADS={threads}
            (virsorter run --keep-original-seq -i {input.metasp} -w {output.folder}  --min-length 1500 --min-score 0.5 --include-groups "dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae" -j {threads} all --scheduler greedy) 2> {log.meta} #Set min length to 1500 to recover as many vContigs as possible
        """