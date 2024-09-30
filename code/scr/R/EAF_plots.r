############################
# Import relevant packages #
############################
source("/projects/luo_lab/Siders_data/code/scr/R/packages.r")

##########################
# Import EAF data frames #
##########################
cell_EAF <- fread("/projects/luo_lab/Siders_data/results/tables/cell_eaf_taxa.csv") %>%
  subset(!is.na(`12C`)) %>%
  subset(!is.na(`13C`))

viral_EAF_ce <- fread("/projects/luo_lab/Siders_data/results/tables/viral_eaf_ce_taxa.csv") %>%
  subset(!is.na(`12C`)) %>%
  subset(!is.na(`13C`)) %>%
  rename("taxa"="lowest_taxonomy")

viral_EAF_ve <- fread("/projects/luo_lab/Siders_data/results/tables/viral_eaf_ve_taxa.csv") %>%
  subset(!is.na(`12C`)) %>%
  subset(!is.na(`13C`))%>%
  rename("taxa"="lowest_taxonomy")

###################################
# Modify data frames for plotting #
###################################

plot_cell_df <- cell_EAF %>%
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
    taxa2
  )) %>%
  mutate(taxa2 = if_else(
    str_starts(taxa2, "SK_"),
    str_replace(str_extract(taxa2, "^SK_[^ ]+"), "^SK_", "D_"),
    taxa2))%>%
  filter(!str_detect(taxa2, '^NR_')) %>%
  filter(!str_detect(taxa2, 'Homo')) %>%
  filter(!str_detect(taxa2, 'Human'))  %>%
  subset(phylum!="") %>%
  group_by(phylum) %>%
  filter(n() > 10) %>%
  ungroup()

length(unique(plot_cell_df$taxa2)) #390
length(unique(plot_cell_df$phylum)) #41
length(unique(plot_cell_df$class)) #110
sort(unique(plot_cell_df$taxa2))
sort(table(plot_cell_df$taxa2))