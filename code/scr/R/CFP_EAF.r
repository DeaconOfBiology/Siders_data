#################################
# Metabolic analysis pathviewer #
#################################
set.seed(03281986)
source("/projects/luo_lab/Siders_data/code/scr/R/packages.r")

bin.KEGG <- fread("/projects/luo_lab/Siders_data/data/processed/metabolic_predictions/MAGs/step_10-visualizeData/combined/counts_KOFam_all_KEGG.tsv") %>%
  pivot_longer(-ID)%>%
  select(name,ID,value)%>%
  rename("KO_ID"="ID")%>%
  subset(value>0)%>%
  select(-value) 

  
bin_ids <-colnames(bin.KEGG)


#Carbon fixation pathways
WL.mod=keggLink("ko","M00377") %>%
  gsub("ko:", "", .) %>%
  as.data.frame() %>%
  rename("KO_ID"='.')%>%
  mutate(pathway_percent_wo_key=ifelse(KO_ID=="K00198"|KO_ID=="K14138",0,1/(nrow(.)-2)),
         acetyl_coa_synthase=ifelse(KO_ID=="K14138",1,0),
         co_dehydrogenase=ifelse(KO_ID=="K00198",1,0))


rTCA.mod=keggLink("ko","M00173")%>%
  gsub("ko:", "", .) %>%
  as.data.frame() %>%
  rename("KO_ID"='.')%>%
  mutate(pathway_percent_wo_key=ifelse(KO_ID=="K15230"|KO_ID=="K15231"
                                       |KO_ID=="K15232"|KO_ID=="K15233"|KO_ID=="K15234",
                                       0,1/(nrow(.)-5)),
         op1_atp_citrate_lyase=ifelse(KO_ID=="K15230"|KO_ID=="K15231",0.5,0),
         op2_citryl_coa_synthetase=ifelse(KO_ID=="K15232"|KO_ID=="K15233"|KO_ID=="K15234",0.333,0))


CBB.mod=keggLink("ko","M00165")%>%
  unique%>%
  gsub("ko:", "", .) %>%
  as.data.frame() %>%
  rename("KO_ID"='.')%>%
  mutate(pathway_percent_wo_key=ifelse(KO_ID=="K01601"|KO_ID=="K01602",0,1/(nrow(.)-2)),
         large=ifelse(KO_ID=="K01601",1,0),
         small=ifelse(KO_ID=="K01602",1,0))


#Find bins with the three carbon fixation pathways
bin_CBB<- inner_join(bin.KEGG, CBB.mod) %>%
  select(-KO_ID)%>%
  group_by(name)%>% 
  summarise_if(is.numeric, sum, na.rm = TRUE)%>%
  mutate(cfp_present=ifelse(pathway_percent_wo_key>=0.39&large==1,"CBB",NA_character_))%>%
  subset(cfp_present=="CBB")%>%
  select(name,cfp_present)
  
bin_WL<- inner_join(bin.KEGG, WL.mod)%>%
  select(-KO_ID)%>%
  group_by(name)%>% 
  summarise_if(is.numeric, sum, na.rm = TRUE)%>%
  mutate(cfp_present=ifelse(pathway_percent_wo_key>=0.39&acetyl_coa_synthase==1&co_dehydrogenase==1,"WL",NA_character_))%>%
  subset(cfp_present=="WL")%>%
  select(name,cfp_present)


bin_rtca<- inner_join(bin.KEGG, rTCA.mod)%>%
  select(-KO_ID)%>%
  group_by(name)%>% 
  summarise_if(is.numeric, sum, na.rm = TRUE)%>%
  mutate(cfp_present=ifelse(pathway_percent_wo_key>0.39&(op1_atp_citrate_lyase==1|op2_citryl_coa_synthetase==1),"rTCA",NA_character_))%>%
  subset(cfp_present=="rTCA")%>%
  select(name,cfp_present)


auto_trophs<- rbind(bin_CBB,bin_rtca,bin_WL)
heterotrophs <- bin.KEGG%>%
  filter(!(name %in% auto_trophs$name))%>%
  mutate(cfp_present="heterotroph")%>%
  select(name,cfp_present)%>%
  unique

trophic_level <- rbind(auto_trophs,heterotrophs)%>%
  mutate(name=sub("bin.","bin_",name)) %>%
  group_by(name) %>%
  summarise(cfp_present = str_c(cfp_present, collapse = ", ")) %>%
  ungroup()

#Now, add EAF values to the trophic_level df
cell_eaf <- fread("/projects/luo_lab/Siders_data/results/tables/cell_eaf_taxa.csv") %>%
  subset(!is.na(`12C`)) %>%
  subset(!is.na(`13C`)) %>%
  dplyr::rename("contig"="organism",
         "name"="MAG")%>%
  select(name,contig,EAF,taxa)%>%
  unique%>%
  mutate(taxa2 = if_else(
    str_detect(taxa, "\\(NCBI: [^)]+\\)"),
    str_extract(taxa, "(?<=\\(NCBI: )[^)]+"),
    taxa)) %>%
  mutate(taxa2 = if_else(
    str_starts(taxa2, "S_") & str_detect(taxa2, "uncultured"),
    str_replace(taxa2, "^S_uncultured (.+)$", "G_\\1"),
    if_else(
      str_starts(taxa2, "S_"),
      str_replace(str_extract(taxa2, "^S_[^ ]+"), "^S_", "G_"),
      taxa2))) %>%
  mutate(taxa2 = if_else(
    str_starts(taxa2, "ST_"),
    str_replace(str_extract(taxa2, "^ST_[^ ]+"), "^ST_", "G_"),
    taxa2)) %>%
  mutate(taxa2 = if_else(
    str_starts(taxa2, "SK_"),
    str_replace(str_extract(taxa2, "^SK_[^ ]+"), "^SK_", "D_"),
    taxa2))%>%
  filter(!str_detect(taxa2, '^NR_')) %>%
  filter(!str_detect(taxa2, 'Homo')) %>%
  filter(!str_detect(taxa2, 'Human')) %>%
  mutate(taxa2 = str_replace(taxa2, "^(.)_", "(\\1)"))


trophic_EAF <- inner_join(trophic_level,cell_eaf)%>%
  select(name,contig,taxa2,cfp_present,EAF)
  
counts <- trophic_EAF[,c("name","contig")] %>%
  group_by(name) %>%
  dplyr::summarize(num_contigs = n()) %>%
  filter(num_contigs >= 10) %>%
  pull(name)

#Remove MAGs with less than 10 contigs
trophic_EAF_filter <- trophic_EAF%>%
  filter(name %in% counts) 

taxon_pool <- trophic_EAF_filter %>%
  group_by(taxa2)%>%
  summarise(pool=mean(EAF)<0.019,
            mean=mean(EAF),
            .groups = "drop")

trophic_plot_df <- inner_join(trophic_EAF_filter, taxon_pool)%>%
  mutate(taxa2=if_else(pool,"Other",taxa2)) 
  
  


trophic_plot_df$cfp_present <- factor(trophic_plot_df$cfp_present, levels=c("CBB","rTCA","WL","CBB, rTCA","heterotroph"))