source("/projects/luo_lab/Siders_data/code/scr/R/packages.r")
source("/projects/luo_lab/Siders_data/code/scr/R/EAF.r")
source("/projects/luo_lab/Siders_data/code/scr/R/CFP_EAF.r")
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

colors <- c(
  "#1F78B4", "#33A02C", "#E31A1C", "#FF7F00", "#6A3D9A",
  "#B15928", "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F",
  "#CAB2D6", "#FFFF99", "#E41A1C", "#377EB8", "#4DAF4A",
  "#FF7F00", "#6A3D9A", "#B15928", "#FBB4AE", "#B2DF8A",
  "#99C794", "#FFD92F", "#E41A1C", "#FFBB78", "#D9D9D9",
  "#BC80BD", "#FFFF99", "#FF7F00", "#66C2A5", "#FC8D62",
  "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#D9D9D9",
  "#F7F7F7", "#B2DF8A", "#FC8D62", "#8DA0CB", "#E78AC3",
  "#A6D854", "#FFD92F", "#FF7F00", "#A6D854", "#D9D9D9",
  "#BFD3C1", "#F6EB61", "#E31A1C", "#33A02C", "#1F78B4",
  "#6A3D9A", "#B15928", "#D95F02", "#7570B3", "#E41A1C",
  "#F4A582", "#D9D9D9"
)

names(colors3) <- unique(plot_cell_df$taxa2)
names(colors) <- unique(plot_cell_df$phylum)

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
###########
# Figures #
###########

#Phylum level
##all
cell_taxa_p_all <- plot_cell_df %>%
  ggplot(.,aes(x=`12C`,y=EAF, color=phylum))+
  geom_point(aes(shape = ifelse(MAG == "unbinned", "Unbinned", "Binned")))+
  geom_hline(yintercept = threshold, color="darkred",linetype='dashed',linewidth=1)+
  scale_shape_manual(values = c("Binned" = 19, "Unbinned" = 23)) +
  theme(legend.position = "none")+
  scale_color_manual(values=colors) +
  ggtitle("Prokaryotic contigs phylum level")+
  xlab("DNA density (12C)")+
  ylab("Excess Atom Fraction")+
  theme_classic()+ 
  theme(legend.title = element_blank())

cell_taxa_p_all
ggsave(file="/projects/luo_lab/Siders_data/results/figures/eaf_12C_phylum_all.pdf", cell_taxa_p_all, width = 15,height = 8)

##binned
cell_taxa_p_binned <- plot_cell_df%>%
  subset(phylum!=""&MAG!="unbinned") %>%
  ggplot(.,aes(x=`12C`,y=EAF, color=phylum))+
  geom_point(aes(shape = ifelse(MAG == "unbinned", "Unbinned", "Binned")))+
  geom_hline(yintercept = threshold, color="darkred",linetype='dashed',linewidth=1)+
  scale_shape_manual(values = c("Binned" = 19, "Unbinned" = 23)) +
  theme(legend.position = "none")+
  scale_color_manual(values=colors) +
  ggtitle("Prokaryotic contigs phylum level (binned)")+
  xlab("DNA density (12C)")+
  ylab("Excess Atom Fraction")+
  theme_classic()+ 
  theme(legend.title = element_blank())

cell_taxa_p_binned
ggsave(file="/projects/luo_lab/Siders_data/results/figures/eaf_12C_phylum_binned.pdf", cell_taxa_p_binned, width = 15,height = 8)

##unbinned
cell_taxa_p_unbinned<- plot_cell_df%>%
  subset(phylum!=""&MAG=="unbinned") %>%
  ggplot(.,aes(x=`12C`,y=EAF, color=phylum))+
  geom_point(aes(shape = ifelse(MAG == "unbinned", "Unbinned", "Binned")))+
  geom_hline(yintercept = threshold, color="darkred",linetype='dashed',linewidth=1)+
  scale_shape_manual(values = c("Binned" = 19, "Unbinned" = 23)) +
  theme(legend.position = "none")+
  scale_color_manual(values=colors) +
  ggtitle("Prokaryotic contigs phylum level (unbinned)")+
  xlab("DNA density (12C)")+
  ylab("Excess Atom Fraction")+
  theme_classic()+ 
  theme(legend.title = element_blank())

cell_taxa_p_unbinned
ggsave(file="/projects/luo_lab/Siders_data/results/figures/eaf_12C_phylum_unbinned.pdf", cell_taxa_p_unbinned, width = 15,height = 8)

#####################################################################
# Plot a facet scatter plot at the phylum level of the cell contigs #
#####################################################################
cell_phylum_p_facet <- plot_cell_df%>%
  ggplot(.,aes(x=`12C`,y=EAF, color=phylum))+
  geom_point()+ 
  geom_hline(yintercept = threshold, color="darkred",linetype='dashed',linewidth=1)+
  facet_wrap(~ phylum)+ 
  scale_color_manual(values=colors) +
  ggtitle("Prokaryotic contigs faceted by taxa")+
  xlab("DNA density (12C)")+
  ylab("Excess Atom Fraction")+
  theme_classic() +
  theme(legend.position = "none",
  text = element_text(size = 15)) 
cell_phylum_p_facet
ggsave(file="/projects/luo_lab/Siders_data/results/figures/eaf_12C_phylum_facet.pdf", cell_phylum_p_facet, width = 30,height = 30)

###############################
# Plot just the viral contigs #
###############################
#Viral contigs in cell enrichment
viral_taxa_ce_p <- viral_EAF_ce%>%
  mutate(taxa=ifelse(taxa=="","Unclassified",taxa))%>%
  ggplot(.,aes(x=`12C`,y=EAF, color=taxa))+
  geom_point()+ 
  geom_hline(yintercept = threshold, color="darkred",linetype='dashed',linewidth=1)+
  scale_color_manual(values=viral_colors_ce) +
  ggtitle("vOTU contigs in cell-enrichment")+
  xlab("DNA density (12C)")+
  ylab("Excess Atom Fraction")+
  theme_classic()  
viral_taxa_ce_p
ggsave(file="/projects/luo_lab/Siders_data/results/figures/eaf_12C_viral_ce.pdf", viral_taxa_ce_p, width = 8,height = 8)

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
viral_taxa_ve_p
ggsave(file="/projects/luo_lab/Siders_data/results/figures/eaf_12C_viral_ve.pdf", viral_taxa_ve_p, width = 8,height = 8)

################################################################
# Plot viral contigs in cell enrichement onto the cell contigs #
################################################################
cell_viral_p <-
  ggplot()+
  geom_point(data=plot_cell_df,aes(x=`12C`,y=EAF, color=phylum,
    shape = ifelse(MAG == "unbinned", "triangle", "circle")))+
  scale_shape_manual(values = c("circle" = 19, 
    "triangle" = 23)) +  # Use 21 for filled circles, 24 for filled triangles
  scale_color_manual(values=colors) +
  geom_point(data=viral_EAF_ce, aes(x=`12C`,y=EAF),shape=8,
    size=4)+
  geom_hline(yintercept = threshold, 
    color="darkred",linetype='dashed',linewidth=1)+
  ggtitle("Prokaryotic and vOTU contigs in cell-enrichment")+
  xlab("DNA density (12C)")+
  ylab("Excess Atom Fraction")+
  theme_classic() +
  theme(legend.position = "none") 
cell_viral_p
ggsave(file="/projects/luo_lab/Siders_data/results/figures/cell_virus_eaf_12C_taxa.pdf", cell_viral_p, width = 8,height = 8)


################
# CFP EAF Plot #
################
#Find uniqueish colors
P33 = createPalette(33,  c("#ff0000", "#00ff00", "#0000ff"))%>%
  unname()

eaf_limit_rect <- data.frame(xmin=-0.05, xmax=0.04907736, ymin=-Inf, ymax=12)
phototroph_astrick <- data.frame(
  x = c(-0.02),
  y = c(3)
)

cfp_eaf_p <- trophic_plot_df%>%
  ggplot(aes(y=reorder(taxa2,EAF),x=EAF,color = cfp_present))+
  # geom_bar(stat="identity")+
  geom_rect(data=eaf_limit_rect, aes(xmin=-Inf, xmax=xmax, ymin=ymin, ymax=ymax), fill="darkgrey",
            alpha=0.1,
            inherit.aes = FALSE)+
  stat_summary(fun.data = median_hilow, geom="pointrange",
               fun.args = list(conf.int=0.5),
               position = position_dodge(width = 0.5),
               size=0.5,
               linewidth=0.75
               #show.legend = FALSE
               )+
  geom_point(data=phototroph_astrick, aes( x=x, y=y), shape = 8,
             size =3,
             color="black")+
  scale_color_manual(values=c("#DC267F","#FE6100",
                              "#785EF0", "#FFB000",
                              "black")) +
  theme_classic()+
  guides(fill=guide_legend(title=""))+
  xlab("Total MAG 13C assimilation")+
  theme(axis.title.y = element_blank(),
        axis.text = element_text(size=8),
        axis.text.y = element_text(face = "italic"),
        legend.text = element_text(size = 8),
        legend.position = c(0.85, 0.2))
cfp_eaf_p
ggsave("/projects/luo_lab/Siders_data/results/figures/CFP_EAF_pointrage.pdf",plot = cfp_eaf_p, width = 5,height = 5)
