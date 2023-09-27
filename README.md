# Microbiome data analysis
  Experimental design and data analysis was conducted by Youn-Jeong Choi.
  
## Strategic questions for microbiome data analysis 
Q1. Does SD (Sphingadiene, chemopreventing agent) treatment increase the diversity of microbial communities and selectively enrich for species that outcompete pathogenic species in the gut? 
     
Q2. Does SPL (Sphingosine-1-phophate Lyase) Knockout reduce the diversity of microbial communities and increase the number of pathogens and pathogen growth in the gut? 
 
Q3. How does SD treatment or SPL Knockout affect to the microbial communities in the feces vs ileum samples? 

## Experimental schemes
Experiment 1: SD (or vehicle only) was administered by gavage to wild type mice daily for consecutive 7 days.  Feces were collected on one day before SD treatment and one day after treatment and mice were sacrificed immediately after the last fecal collection and ileum were collected. 
 
Experiment 2: gut SPL disruption was induced by daily NF injection intra-peritoneally for consecutive 5 days.  Feces were collected on one day before the start of NF injection and collected on day 14 after completion of NF injection.  Mice were sacrificed after fecal collection and Ileum was collected. 

## Key outputs of Data analysis and statistics 
1. Alpha diversity metrics: Shannon index, chao1, phylogenetic diversity

2. Beta diversity metrics: bray curtis, weighted unifrac, unweighted unifrac, canberra 

3. Relative abundance by group

4. Differential abundance of taxa by group

## The order of data processing and analysis
1. Upload a data folder on HPC, run fastQCs, and demultiplex sequences: described in split_libraries.rtf

2. Process data, annotate sequences, assign ASV, generate ASV tables: 16S_dada2.R

3. Filter non-bacterial and low variant species, determine alpha rarefaction, remove negative control, rarefy the data, measure alpha diversity and beta diversity: described in qiime.rtf

4. Plots and statistics of alpha diversity and beta diversity: alpha_beta_diversity.R

5. Differential abundance of taxa was identified using poisson, negative binomial, zero-inflated negative binomial models.
