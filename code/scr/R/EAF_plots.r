############################
# Import relevant packages #
############################
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggtext)
library(RColorBrewer)
library(gridExtra)

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
  group_by(taxa2) %>%
  filter(n() >= 10) %>%
  ungroup()

##########################
# Set up plot parameters #
##########################
#Parameters for cell contigs
threshold<-0.04907736
rect <- data.frame(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=0.06576442,
                   gxmin=-Inf, gxmax=Inf, gymin=0.06576442, gymax=Inf)
colors3 <-c("royalblue2","cyan3","steelblue1","lightskyblue1","cadetblue","blue4",
            "gold4","darkorange3","darkmagenta","hotpink1","mediumorchid3","deeppink",           
            "lightskyblue4","honeydew4","navajowhite2","firebrick4","red2","deeppink4",
            "plum4","yellow1","gold1","darkseagreen3","grey","brown","thistle",          
            "goldenrod4","purple1","mediumpurple3","lightslateblue","olivedrab",
            "seagreen","chartreuse1","darkolivegreen1","yellow4","olivedrab3",
            "darkgoldenrod2","darkorange","rosybrown1","burlywood2","orange4",
            "lightsalmon3","coral1","darkseagreen4","lemonchiffon4")

names(colors3) <- unique(plot_cell_df$taxa2)

#Parameters for viral contigs in the cell enrichment
viral_colors_ce <- c("darkmagenta","blue4","cyan3",
  "darkorange3","chartreuse1","yellow4",
  "darkseagreen4","plum4","lightskyblue1")
names(viral_colors_ce) <- unique(viral_EAF_ce$taxa)

#Parameters for viral contigs in the viral enrichment
viral_colors_ve <-c("royalblue2","cyan3","steelblue1","lightskyblue1","cadetblue","blue4",
            "gold4","darkorange3","darkmagenta","hotpink1","mediumorchid3","deeppink",           
            "lightskyblue4","honeydew4","navajowhite2","firebrick4","red2","deeppink4",
            "plum4","yellow1","gold1","darkseagreen3","grey","brown","thistle",          
            "goldenrod4","purple1","mediumpurple3","lightslateblue","olivedrab",
            "seagreen","chartreuse1","darkolivegreen1","yellow4","olivedrab3")
names(viral_colors_ve) <- unique(viral_EAF_ve$taxa)

##############################
# Plot just the cell contigs #
##############################
cell_taxa_p <- plot_cell_df%>%
  ggplot(.,aes(x=`12C`,y=EAF, color=taxa2))+
  geom_point(aes(shape = ifelse(MAG == "unbinned", "triangle", "circle")))+
  geom_hline(yintercept = threshold, color="darkred",linetype='dashed',linewidth=1)+
  scale_shape_manual(values = c("circle" = 19, "triangle" = 23)) +  # Use 21 for filled circles, 24 for filled triangles
  scale_color_manual(values=colors3) +
  ggtitle("Prokaryotic contigs")+
  xlab("DNA density (12C)")+
  ylab("Excess Atom Fraction")+
  theme_classic() +
  theme(legend.position = "none") 
#cell_taxa_p
ggsave(file="/projects/luo_lab/Siders_data/submission/test_figures/eaf_12C_taxa.pdf", cell_taxa_p, width = 8,height = 8)

#####################################################################
# Plot a facet scatter plot at the phylum level of the cell contigs #
#####################################################################
cell_taxa_p_facet <- plot_cell_df%>%
  ggplot(.,aes(x=`12C`,y=EAF, color=taxa2))+
  geom_point()+ 
  geom_hline(yintercept = threshold, color="darkred",linetype='dashed',linewidth=1)+
  facet_wrap(~ taxa2)+ 
  scale_color_manual(values=colors3) +
  ggtitle("Prokaryotic contigs faceted by taxa")+
  xlab("DNA density (12C)")+
  ylab("Excess Atom Fraction")+
  theme_classic() +
  theme(legend.position = "none") 
#cell_taxa_p_facet
ggsave(file="/projects/luo_lab/Siders_data/submission/test_figures/eaf_12C_taxa_facet.pdf", cell_taxa_p_facet, width = 15,height = 15)

###############################
# Plot just the viral contigs #
###############################
#Viral contigs in cell enrichment
viral_taxa_ce_p <- viral_EAF_ce%>%
  ggplot(.,aes(x=`12C`,y=EAF, color=taxa))+
  geom_point()+ 
  geom_hline(yintercept = threshold, color="darkred",linetype='dashed',linewidth=1)+
  scale_color_manual(values=viral_colors_ce) +
  ggtitle("vOTU contigs in cell-enrichment")+
  xlab("DNA density (12C)")+
  ylab("Excess Atom Fraction")+
  theme_classic() +
  theme(legend.position = "none") 
#viral_taxa_ce_p
ggsave(file="/projects/luo_lab/Siders_data/submission/test_figures/eaf_12C_viral_ce.pdf", viral_taxa_ce_p, width = 8,height = 8)

#Viral contigs in viral enrichment
viral_taxa_ve_p <- viral_EAF_ve%>%
  ggplot(.,aes(x=`12C`,y=EAF, color=taxa))+
  geom_point()+ 
  scale_color_manual(values=viral_colors_ve) +
  ggtitle("vOTU contigs in viral-enrichment")+
  xlab("DNA density (12C)")+
  ylab("Excess Atom Fraction")+
  theme_classic() +
  theme(legend.position = "none") 
#viral_taxa_ve_p
ggsave(file="/projects/luo_lab/Siders_data/submission/test_figures/eaf_12C_viral_ve.pdf", viral_taxa_ve_p, width = 8,height = 8)

################################################################
# Plot viral contigs in cell enrichement onto the cell contigs #
################################################################
cell_viral_p <-
  ggplot()+
  geom_point(data=plot_cell_df,aes(x=`12C`,y=EAF, color=taxa2,
    shape = ifelse(MAG == "unbinned", "triangle", "circle")))+
  scale_shape_manual(values = c("circle" = 19, 
    "triangle" = 23)) +  # Use 21 for filled circles, 24 for filled triangles
  scale_color_manual(values=colors3) +
  geom_point(data=viral_EAF_ce, aes(x=`12C`,y=EAF),shape=8,
    size=4)+
  geom_hline(yintercept = threshold, 
    color="darkred",linetype='dashed',linewidth=1)+
  ggtitle("Prokaryotic and vOTU contigs in cell-enrichment")+
  xlab("DNA density (12C)")+
  ylab("Excess Atom Fraction")+
  theme_classic() +
  theme(legend.position = "none") 
#cell_viral_p
ggsave(file="/projects/luo_lab/Siders_data/submission/test_figures/cell_virus_eaf_12C_taxa.pdf", cell_viral_p, width = 8,height = 8)
