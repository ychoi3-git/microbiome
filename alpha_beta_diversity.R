## plotting alpha diversity and beta diversity using ileum samples.
# I used pipelines I learned in workshops (Youn-Jeong Choi).

library("tidyverse")

# alpha diversity plotting and statistics
metadata1 <- read.table("WT_KO_I_data/Miseq_08062020_metadata_KO_I_rarefy_nohash.txt", header=TRUE, check.names=FALSE, sep="\t", comment="")
metadata2 <- select(metadata1, -BarcodeSequence, -LinkerPrimerSequence, -Description)
adiv <- read.table("WT_KO_I_data/adiv_chao1_pd_5000.txt", header=T,comment="",sep="\t", check.names=F)
adiv_meta <- left_join(adiv, metadata2, by="SampleID") 

ggplot(adiv_meta, aes (x=Treatment, y=shannon))+
  geom_boxplot()+
  geom_point()+  
  ylab("Shannon diversity") + 
  xlab("Genotype")

ggplot(adiv_meta, aes (x=Treatment, y=PD_whole_tree))+
  geom_boxplot()+
  geom_point()+  
  ylab("Faith's phylogenetic Diversity") + 
  xlab("Genotype")

ggplot(adiv_meta, aes (x=Treatment, y=chao1))+
  geom_boxplot()+
  geom_point()+  
  ylab("chao1") + 
  xlab("Genotype")
ggsave("WT_KO_I_data/figure/shannon.pdf", height=3, width=4)
ggsave("WT_KO_I_data/figure/PD_whole_tree.pdf", height=3, width=4)
ggsave("WT_KO_I_data/figure/chao1.pdf", height=3, width=4)

wilcox.test(shannon~Treatment, data=adiv_meta)
wilcox.test(PD_whole_tree~Treatment, data=adiv_meta)
wilcox.test(chao1~Treatment, data=adiv_meta)

#load packages for beta diversity plotting and statistics
library(tidyverse)
library(vegan)
library(ape)
library(lmerTest)
library(broom)
library(dada2)
library(ggtree)
library(devtools)

##beta diversity

MYMetaf <- read.table("WT_KO_I_data/Miseq_08062020_metadata_KO_I_rarefy_nohash.txt", header=TRUE, check.names=FALSE, sep="\t", comment="")
MYMetaf2 <- select(MYMetaf, -BarcodeSequence, -LinkerPrimerSequence, -Description)
dfile <- read.table("WT_KO_I_data/bdiv_even_5000/unweighted_unifrac_dm.txt", header=TRUE, check.names=FALSE, sep="\t", comment.char="", row.names=1)
dfile2 <- read.table("WT_KO_I_data/bdiv_even_5000/weighted_unifrac_dm.txt", header=TRUE, check.names=FALSE, sep="\t", comment.char="", row.names=1)
dfile3 <- read.table("WT_KO_I_data/beta_div_b5000/bray_curtis_rarefied_otutable_5000.txt", header=TRUE, check.names=FALSE, sep="\t", comment.char="", row.names=1)

as.matrix(dfile)[1:10, 1:10]
pcbraycurtis<-pcoa(dfile)
pcbraycurtis<-pcoa(dfile2)
pcbraycurtis<-pcoa(dfile3)
names(pcbraycurtis)
head(pcbraycurtis)

pcbraycurtis$values %>%
  as.data.frame() %>%
  rownames_to_column("PC") %>%
  mutate(PC=as.numeric(PC)) %>%
  select(PC, Relative_eig) %>%
  mutate(PercentVariation=Relative_eig*100) %>%
  ggplot(aes(x=PC, y=PercentVariation)) +
  geom_point()

plotbraycurtis<-
  pcbraycurtis$vectors %>%
  as.data.frame() %>%
  rownames_to_column("SampleID") %>%
  select(SampleID, Axis.1, Axis.2) %>% 
  left_join(MYMetaf2)

head(plotbraycurtis)

plotbraycurtis %>%
  ggplot(aes(x=Axis.1, y=Axis.2, color=Treatment)) +
  geom_point()

plotbraycurtis %>%
  ggplot(aes(x=Axis.1, y=Axis.2, color=Treatment)) +
  geom_point() +
  xlab(paste0("PC1:", round(pcbraycurtis$values$Relative_eig[1]*100, 2),"%")) +
  ylab(paste0("PC2:", round(pcbraycurtis$values$Relative_eig[2]*100, 2),"%")) +
  theme_bw() +
  scale_color_manual(values=c("indianred","cornflowerblue")) +
  ggtitle("Unweighted_unifrac")

plotbraycurtis %>%
  ggplot(aes(x=Axis.1, y=Axis.2, color=Treatment)) +
  geom_point() +
  xlab(paste0("PC1:", round(pcbraycurtis$values$Relative_eig[1]*100, 2),"%")) +
  ylab(paste0("PC2:", round(pcbraycurtis$values$Relative_eig[2]*100, 2),"%")) +
  theme_bw() +
  scale_color_manual(values=c("indianred","cornflowerblue")) +
  ggtitle("Weighted_unifrac")

plotbraycurtis %>%
  ggplot(aes(x=Axis.1, y=Axis.2, color=Treatment)) +
  geom_point() +
  xlab(paste0("PC1:", round(pcbraycurtis$values$Relative_eig[1]*100, 2),"%")) +
  ylab(paste0("PC2:", round(pcbraycurtis$values$Relative_eig[2]*100, 2),"%")) +
  theme_bw() +
  scale_color_manual(values=c("indianred","cornflowerblue")) +
  ggtitle("bray curtis")

ggsave("WT_KO_I_data/figure/unweighted_unifrac.pdf", height=3, width=4)
ggsave("WT_KO_I_data/figure/wunifrac.pdf", height=3, width=4)
ggsave("WT_KO_I_data/figure/bray-curtis.pdf", height=3, width=4)

#statistics for beta diversity
adonis(dfile~Treatment, data = MYMetaf2)
adonis(dfile2~Treatment, data = MYMetaf2)
adonis(dfile3~Treatment, data = MYMetaf2)
