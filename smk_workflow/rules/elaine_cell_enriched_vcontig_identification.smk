# For this, lets consider steps in https://www.protocols.io/view/viral-sequence-identification-sop-with-virsorter2-5qpvoyqebg4o/v3?step=3

##################################################
# This is the third snake file that needs to be  #
# generated before the analysis can move on. In  #
# it, we look for viral contigs by:              #
# 1. using VIBRANT to look for viral signatures  #
# 2. using VirSorter2 to look for viral          #
# signatures                                     #
# 3. using ViralVerify to look for viral         #
# signatures                                     #
# 4. using genomad to .....                      #
# 5. using vCheck as an additional quality       #
# check from the above outputs                   #
##################################################

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


# #CheckV for quality control
# rule checkV:
#     input:
#         virsorter_meta="data/processed/Viral_Assemblies/virsorter2/cell-clean_headers_pass1/final-viral-combined.fa",
#     output:
#         folder=directory("data/processed/Viral_Assemblies/checkv/cell-clean_headers"),
#         outfiles=["data/processed/Viral_Assemblies/checkv/cell-clean_headers/proviruses.fna","data/processed/Viral_Assemblies/checkv/cell-clean_headers/viruses.fna"]
#     threads: 16
#     log: 
#         first="logs/cell-{group}_fist_metaspades_checkv.log",
#     shell:
#         """
#             (checkv end_to_end {input.virsorter_meta} {output.folder} -t {threads} -d  /users/troger50/checkv-db-v1.5) 2> {log.first}
#         """
