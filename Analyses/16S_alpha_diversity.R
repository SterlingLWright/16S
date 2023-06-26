#################################################
#
#
#     Alpha diversity analyses for 16S paper
#     Create by Sterling Wright 
#     June 26, 2023 
#
#################################################


# set working directory 
setwd("$PATH")

# load necessary packages 
library(vegan)
library(phyloseq)
library(ggplot2)
library(MicrobiotaProcess)

#### feature analysis ####

set.seed(1024)
features<-import_biom("Feature_Analysis/Exported-Feature_Analysis/absolute_all_merged_FEATURES-json.biom")
metadata<-import_qiime_sample_data("16S-enrichment-MappingFile_15JUN23.txt")
features_phylo<-merge_phyloseq(features,metadata)
features_phylo_filt<-rarefy_even_depth(features_phylo)
rareres <- get_rarecurve(obj=features_phylo_filt, chunks=400)
p_rare <- ggrarecurve(obj=rareres,
                      indexNames=c("Observe","Chao1","Shannon"),
) +
  theme(legend.spacing.y=unit(0.01,"cm"),
        legend.text=element_text(size=4))

prare_features <- ggrarecurve(obj=rareres, factorNames="LibraryType",
                      indexNames=c("Observe", "Chao1", "Shannon")
) +
  scale_fill_manual(values=c("#00AED7", "#FD9347", "blue2"))+
  scale_color_manual(values=c("#00AED7", "#FD9347", "blue2"))+
  theme_bw()+
  theme(axis.text=element_text(size=14), panel.grid=element_blank(),
        strip.background = element_rect(colour=NA,fill="grey"),
        strip.text.x = element_text(face="bold"))          

prare2 <- ggrarecurve(obj=rareres,
                      factorNames="Group",
                      shadow=FALSE,
                      indexNames=c("Observe", "Chao1", "ACE")
) +
  scale_color_manual(values=c("#00AED7", "#FD9347"))+
  theme_bw()+
  theme(axis.text=element_text(size=8), panel.grid=element_blank(),
        strip.background = element_rect(colour=NA,fill="grey"),
        strip.text.x = element_text(face="bold"))
p_rare / prare1 / prare2


#### genus analysis ####
genus<-import_biom("Genus_Analysis/Exported-Genus_Analysis/absolute_all_merged_GENERA-json.biom")
genus_phylo<-merge_phyloseq(genus,metadata)
genus_phylo_filt<-rarefy_even_depth(genus_phylo)
rareres2 <- get_rarecurve(obj=genus_phylo_filt, chunks=400)
p_rare2 <- ggrarecurve(obj=rareres2,
                      indexNames=c("Observe","Shannon","ACE"),
) +
  theme(legend.spacing.y=unit(0.01,"cm"),
        legend.text=element_text(size=4))

prare1_genus <- ggrarecurve(obj=rareres2, factorNames="LibraryType",
                      indexNames=c("Observe", "Chao1", "Shannon")
) +
  scale_fill_manual(values=c("#00AED7", "#FD9347", "blue2"))+
  scale_color_manual(values=c("#00AED7", "#FD9347", "blue2"))+
  theme_bw()+
  theme(axis.text=element_text(size=14), panel.grid=element_blank(),
        strip.background = element_rect(colour=NA,fill="grey"),
        strip.text.x = element_text(face="bold"))          

prare_features / prare1_genus
