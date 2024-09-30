cell_eaf_12C3 <- fread("/projects/luo_lab/Siders_data/results/tables/cell_eaf_taxa.csv") %>%
  subset(!is.na(`12C`)) %>%
  subset(!is.na(`13C`)) %>%
  rename("contig"="organism")

original_contig_names <- fread("/projects/luo_lab/Siders_data/data/processed/Assemblies/clean_hearders.txt",
                               header = FALSE, col.names = c("contig", "original_name"))
annotations <- read_tsv("/projects/luo_lab/Siders_data/data/processed/metabolic_predictions/MAGs/combined_final_annotation_summary.tsv") %>%
  filter(!is.na(best_hit))%>%
  subset(grepl("PFAM",HMM)|
           grepl("phrog",HMM)|
           grepl("TIGR",HMM)|
           grepl("KOFam",HMM)) %>%
  as.data.frame %>%
  rename("contig"="target")

pfam <- annotations%>%
  select(contig, best_hit, score)%>%
  pivot_wider(names_from = best_hit, values_from = score)%>%
  mutate(contig = sub("_[^_]*$", "", contig)) %>%
  as.data.table()


pfam_aggregated <- pfam[, lapply(.SD, max, na.rm = TRUE), by = contig] %>%
  inner_join(original_contig_names,.) %>% #Used to get the length of contig.
  extract(original_name, 
          into = "length", 
          regex = "len=(\\d+)", 
          convert = TRUE) %>%
  subset(length>=5000)
write_csv(pfam_eaf, "/projects/luo_lab/Siders_data/results/tables/best_gene_hit_eaf_cell_contigs_5K.csv")

pfam_aggregated <- pfam[, lapply(.SD, max, na.rm = TRUE), by = contig] %>%
  inner_join(original_contig_names,.) %>% #Used to get the length of contig.
  extract(original_name, 
          into = "length", 
          regex = "len=(\\d+)", 
          convert = TRUE) %>%
  subset(length>=10000)

pfam_eaf <- pfam_aggregated %>%
  inner_join(cell_eaf_12C3[,c("contig","EAF")],.)
head(pfam_eaf[,1:10])
write_csv(pfam_eaf, "/projects/luo_lab/Siders_data/results/tables/best_gene_hit_eaf_cell_contigs_10K.csv")
