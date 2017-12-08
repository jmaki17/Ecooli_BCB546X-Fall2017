#  Toward a Metagenomic Understanding on the Bacterial Composition and Resistome in Hong Kong Banknotes


Yoshitaro Heshiki 1, Thrimendra Dissanayake 1, Tingting Zheng 1, Kang Kang 1,Ni Yueqiong 1, Zeling Xu 2, Chinmoy Sarkar 3, Patrick C. Y. Woo 4, 5, 6, 7, 8, Billy K. C. Chow 2,David Baker 2, Aixin Yan 2, Christopher J. Webster 3, Gianni Panagiotou1, 9* and Jun Li 1*



##### Click [here]([https://www.frontiersin.org/articles/10.3389/fmicb.2017.00632/full](https://www.frontiersin.org/articles/10.3389/fmicb.2017.00632/full)) to be directed to paper
---

#  **Contents**

### 1. Introduction of Paper
### 2. Technical Details of Replication Analysis
### 3.  Summary of Replication
#Introduction to Paper
###  &nbsp; &nbsp; &nbsp; &nbsp; Heshiki et al. 2017 wanted to assess to what degree Hong Kong (HK) banknotes serve as a vector for transmitting pathogens and drug resistance throughout the metropolis. They characterized the bacterial community _in vitro_ and by metagenomic analysis that included data on microbial diversity, profiling of potential pathogens and antibiotic-resistance markers on banknotes collected from both hospitals and metro station cashiers. 

### &nbsp; &nbsp; &nbsp; &nbsp;  They used metagenomic shotgun sequencing to profile the bacterial communities of HK banknotes with high resolution. In total 15 samples were used generating 98.5 Gbp worth of data. 22% of these reads (HiSeq) contained DNA of bacterial origin. To highlight the presumed unique microbiome of HK banknotes they did a comparative study of publicly available HK environmental samples, and compared the HK microbiome to that of publicly available Indian banknote profiles.

### Bioinformatic results are as followed:
1. The most abundant bacteria is the genus _Acinetobacter_ (16.27%), which harbours several potential pathogenic species.
2. At the species level _Propionibacterium acnes_ was at the highest abundance
3. The relative abundance of potential pathogenic species was on average 35.6%of the total bacteria of the 15 samples.
4. When compared to environmental samples (water, palms of hands etc..) they found that the HK banknote microbiomes are more diverse and that the banknotes contained around half of the total genera in all the environment samples. They also found that the banknotes have a higher abundance of potential pathogenic species than the environmental samples (p = 9.6e-10). Results also showed that the banknote micobiomes had a larger abundance of antibiotic-resistance genes compared to the environmental samples (p = 6.1e-5)
5. When comparing their samples to that of Indian banknotes they found that HK notes had significantly less abundance of potential pathogens (p = 0.0049)

### Their conclusions
1. Their comparative analysis of bacterial composition, pathogenic content, and resistome profile in banknotes obtained from different regions, facilities and city
network indices revealed no statistically significant differences between sites.
2. Their comparison of HK banknotes with other local environmental samples indicates higher bacterial diversity of banknotes.

#### They claim that these findings together suggest "some evidence that the high circulation of banknotes and their continuous exposure to the microbes that surrounds us may result to one of the most useful surveillance platforms for monitoring the “city’s microbiome status.”" They do mention that much more in-depth sampling is required to confirm their concluding hypothesis. 

# Technical Details of Replication Analysis

### File Acquisition, Inspection, Cleaning

The first step in the technical processing of this data is the acquisition of files by accession number using SRA toolkit, allowing us to acquire them as fastq files.

---

sratoolkit.2.8.2-1-mac64/bin/fastq-dump SRR5312479	#Downloading the hospital reads as a fastq file using the SRA toolkit provided through NCBI

mv SRR5312479.fastq Hospital_11_raw.fastq	# renaming the hospital file to something a little less awful

	# Do same with Metro1:
sratoolkit.2.8.2-1-mac64/bin/fastq-dump SRR5312476
mv SRR5312476.fastq Metrostation_01.fastq

	# Analyzing sequences:
	# To run this you need to have the JDK environment available. To do this, you must go to Java.com and select the “see all java downloads”. On the next page, select the “looking for JDK?” option on the left side of the page. Then select 
curl -O http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.10.1.zip
unzip fastqc_v0.10.1.zip
cd FastQC/
chmod +x fastqc
cd ..

	# installing FASTQC for visualizing sequences and unzipping
	# if you already have this loaded then you can disregard this step

	# Taking a look at the sequence read quality
./FastQC/fastqc Hospital_11_raw.fastq
./FastQC/fastqc Metrostation_01.fastq
	# The hospital_11 reads look to be pretty good quality but tail off slightly at the end, the Metro reads are also fairly good quality but seem to tail off sharply at the end, will attempt to trim off the adaptors to increase the sequence quality

	# installing BBMerge from source forge
wget http://downloads.sourceforge.net/project/bbmap/BBMap_37.68.tar.gz
tar zvxf BBMap_37.68.tar.gz
	# used wget to download bbmap from sourcefourge, then unzipped the file

bbmap/bbduk.sh in=Metrostation_01.fastq out=metro_clean.fq ref=adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo
	# trimming the Metrostation interleaved file, outputting a metro_clean.fq file

./FastQC/fastqc metro_clean.fq
	# compare the fastqc outputs and check for improvement (see a slight improvement between the 2 sets)

	# repeat with the hospital files
bbmap/bbduk.sh in=Hospital_11_raw.fastq out=hosp_clean.fq ref=adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo
./FastQC/fastqc hosp_clean.fq


---


### MetaPhlAn and Relative Abundance

---

wget https://bitbucket.org/nsegata/metaphlan/get/default.tar.bz2
tar xjvf default.tar.bz2
mv *-metaphlan-* metaphlan
	# Downloading MetaPhlan package and decompressing

cd metaphlan
mkdir hkbanknotes
	# Making a new directory in the metaphlan folder to store the fastq files

brew install bowtie2
	# install bowtie2 because it is needed for metaphlan

mv hosp_clean.fq metaphlan/hkbanknotes/
mv metro_clean.fq metaphlan/hkbanknotes/
	# move the files into the hkbanknotes folder in metaphlan

---

---

brew install python
conda search python
conda install python=2.7.8
	# revert to python 2.7.8 (required to run metaphlan)

mkdir metaphlan_results
	# make a directory for the metaphlan output

metaphlan joelmaki$ $ python metaphlan.py hkbanknotes/hosp_clean.fq --bowtie2db bowtie2db/mpa --nproc 5 --bowtie2out hosp_metagenome.bt2out.txt > metaphlan_results/HK_hosp_real.txt
	# This command aligns the reads to the metahlan’s database, which is then able to generate an abundance table giving the relative abundances of the various OTUs within the datatable. This is nice because this abundance table can be used to generate figures (rarefaction curves, etc.) using the R package Vegan. Metaphlan uses the program bowtie2 to align the sequences with the reference database. Metaphlan requires you to be working in a python >2.6 environment but not python 3 (issues with tabs).

	# Repeat script with the metrostation fastq this time

wc -l hosp_metagenome.bt2out.txt
	# determining the number of taxon “hits” in the dataset for creation of an abundance table in R (2070)


metaphlan joelmaki $ python metaphlan.py hkbanknotes/metro_clean.fq --bowtie2db bowtie2db/mpa --nproc 5 --bowtie2out metro_metagenome.bt2out.txt > metaphlan_results/HK_metro_real.txt

wc -l metro_metagenome.bt2out.txt
	# determining the number of taxon “hits” in the dataset for creation of an abundance table in R (157815)


mkdir figures
	# make a directory for the figures about to be generated 
utils/merge_metaphlan_tables.py metaphlan_results/*.txt > figures/merged_abundance_table.txt
	# produce an output table that contains the relative abundances with the otus as rows and the sample as the columns

plotting_scripts/metaphlan_hclust_heatmap.py -c bbcry --top 25 -s log --in figures/merged_abundance_table.txt --out figures/abundance_heatmap.png
	# generates (somewhat) of a heatmap. The dendrogram of the microbes present is facing the wrong way. This seems to have been an issue that hasn’t been resolved since the package hasn’t been updated since 2014. We still get somewhat of a heatmap though so that’s pretty neat.

plotting_scripts/metaphlan_hclust_heatmap.py -c bbcry -f braycurtis --top 25 -s log --in figures/merged_abundance_table.txt --out figures
abundance_heatmap_2.png
-got the dendrogram lines to show up, but they’re still backwards (might have to go into photoshop and mirror the dendrogram) <- this is exactly what I did

	# Hospital 6 scripts:
	# After going through the analysis, we decided to arbitrarily run another (small) sample to generate a better species accumulation curve and to determine whether our initial (hospital 11) sample was corrupted (it gave a very small output)

sratoolkit.2.8.2-1-mac64/bin/fastq-dump SRR5312484
	# downloading the hosp_6 reads using the sra toolkit
mv SRR5312484.fastq hosp6.fastq
	# Renaming the SRR for hospital 6

./FastQC/fastqc hosp6.fastq
	# checking the quality of the hosp6 reads

bbmap/bbduk.sh in=hosp6.fastq out=hosp6_clean.fq ref=adapters.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo
	# trimming the fastq seqs using the bbduk function in bbmap (this should remove the adaptors from the Illumina sequencing and trim any low-quality reads)

mv hosp6_clean.fq metaphlan/hkbanknotes/
	# moving hosp6_clean.fq into the hkbanknotes folder within metaphlan

cd metaphlan
	# change working directory to metaphlan

python metaphlan.py hkbanknotes/hosp6_clean.fq --bowtie2db bowtie2db/mpa --nproc 5 --bowtie2out hosp6_metagenome.bt2out.txt > metaphlan_results/HK_hosp6_real.txt
	# running metaphlan on the hosp6 data

utils/merge_metaphlan_tables.py metaphlan_results/*.txt > figures/merged_abundance_table.txt
	# merging the abundance table from hosp6 with hosp11 and the metrostation

plotting_scripts/metaphlan_hclust_heatmap.py -c bbcry -f braycurtis --top 25 -s log --in figures/merged_abundance_table.txt --out figures/abundance_heatmap_3.png
	# plotting the heatmap dendrogram, the dendrogram is still backwards but it looks right so that’s neat

wc -l hosp6_metagenome.bt2out.txt
	# determining the number of taxon “hits” in the dataset for creation of an abundance table in R (9048)

---

### Preparation of Files for VEGAN

---

library(dplyr)
library(tidyr)
library(stringr)
install.packages("splitstackshape")
library(splitstackshape)
	## importing dplyr and tidyr to get the files into the right format
hosp <- as.data.frame(HK_hosp_real)
colnames(hosp) <- c("taxa", "per_abund")
	## adding column headings
hosp_columns <- cSplit(hosp, "taxa", "|")
	## splitting the string "taxa" into individual columns for kingdom, phylum....species
na.omit(hosp_columns) -> hosp_columns
hosp_columns$abundance <- hosp_columns$per_abund*.01*2070
	## extrapolating the abundance counts from the abundance percentages and the number
	## of taxa "hits"
round(hosp_columns$abundance) ->hosp_columns$abundance
	## rounding the count to the nearest whole number
hosp_columns <- hosp_columns[ ,-c(1:7)]
	## removing all the columns that are not abundance and species

trans_hosp = data.frame(t(hosp_columns[,-1]))
colnames(trans_hosp) <- hosp_columns[,hosp_columns$taxa_7] 
row.names(trans_hosp) <- "Hospital"
##

	## Repeat with the Metro Data
metro <- as.data.frame(HK_metro_real)
colnames(metro) <- c("taxa", "per_abund")
	## adding column headings
metro_columns <- cSplit(metro, "taxa", "|")
	## splitting the string "taxa" into individual columns for kingdom, phylum....species
na.omit(metro_columns) -> metro_columns
metro_columns$abundance <- metro_columns$per_abund*.01*157815
	## extrapolating the abundance counts from the abundance percentages and the number
	## of taxa "hits"
round(metro_columns$abundance) ->metro_columns$abundance
	## rounding the count to the nearest whole number
metro_columns <- metro_columns[ ,-c(1:7)]

trans_metro = data.frame(t(metro_columns[,-1]))
colnames(trans_metro) <- metro_columns[,metro_columns$taxa_7]
row.names(trans_metro) <- "Metro"

	## merging the 2 dataframes into 1 super dataframe for species accumulation goodness
fastmerge <- function(d1, d2) {
  d1.names <- names(d1)
  d2.names <- names(d2)
  
  # columns in d1 but not in d2
  d2.add <- setdiff(d1.names, d2.names)
  
  # columns in d2 but not in d1
  d1.add <- setdiff(d2.names, d1.names)
  
  # add blank columns to d2
  if(length(d2.add) > 0) {
    for(i in 1:length(d2.add)) {
      d2[d2.add[i]] <- NA
    }
  }
  
  	# add blank columns to d1
  if(length(d1.add) > 0) {
    for(i in 1:length(d1.add)) {
      d1[d1.add[i]] <- NA
    }
  }
  
  return(rbind(d1, d2))
}
	## Using a fastmerge function I found on stack overflow
fastmerge(trans_hosp, trans_metro) ->merged_df
	## generating a merged abundance table for importing into vegan

write.csv(merged_df,'hkbn_abundances.csv')

adding metro environmental samples to abundance data
Env <- as.data.frame(Metro_environmental_abundance)
Env_columns <- cSplit(Env, "X1", ",")
	## splitting the column that contains both the abundance data and the kingdom level

Env_columns[Env_columns$X7 == "s__"] <- NA
	## Changing rows not at the species level to NA's
Env_columns$species<- with(Env_columns, paste0(X6, X7))
Env_columns[Env_columns$species == "NANA"] <- NA
na.omit(Env_columns) -> Env_columns
Env_columns$species<-sub("s__","_", Env_columns$species)
Env_columns$species<-sub("g","s", Env_columns$species)
Env_columns$species<-sub("[*]","*", Env_columns$species)
	## removing non-species level data and getting the naming conventions the same
	## between the 2 files
Env_columns_new <- Env_columns[ ,-c(1:6,8)]
Env_columns_new[,c("species", "X1_1")] -> Env_columns_new
	## removing all the columns that are not species or abundance and rearranging
trans_env = data.frame(t(Env_columns_new[,-1]))
	## transposing the environmental samples
colnames(trans_env) <- Env_columns[,Env_columns$species]
	## setting column names to the species

row.names(trans_env) <- "Environment"
	## changing the row names


	## removing names from the column headings of the trans_env_cut df
library(plyr)

rbind.fill(trans_env, merged_df) ->env_merged_df
	## Using rbind.fill to append the environmental samples to the merged dataframe
	## containing the merged metro and hospital abundance data
rownames(env_merged_df) <- c("Environment", "Hospital", "Metro")
	## changing the row names to reflect the datasets they are associated with
write.csv(env_merged_df,'hk_metro_air_environment_abundances.csv')


	## adding hospital6 data in
library(readr)
HK_hosp6_real <- read_delim("~/metaphlan/metaphlan_results/HK_hosp6_real.txt", 
                            "\t", escape_double = FALSE, trim_ws = TRUE)
View(HK_hosp6_real)
	## importing hospital6 dataframe
hosp6 <- as.data.frame(HK_hosp6_real)
colnames(hosp6) <- c("taxa", "per_abund")
	## adding column headings
hosp6_columns <- cSplit(hosp6, "taxa", "|")
	## splitting the string "taxa" into individual columns for kingdom, phylum....species
na.omit(hosp6_columns) -> hosp6_columns
hosp6_columns$abundance <- hosp6_columns$per_abund*.01*9048
	## extrapolating the abundance counts from the abundance percentages and the number
	## of taxa "hits" from the unix command
round(hosp6_columns$abundance) ->hosp6_columns$abundance
	## rounding the count to the nearest whole number
hosp6_columns <- hosp6_columns[ ,-c(1:7)]

trans_hosp6 = data.frame(t(hosp6_columns[,-1]))
colnames(trans_hosp6) <- hosp6_columns[,hosp6_columns$taxa_7]
row.names(trans_hosp6) <- "Hospital_6"

	## merging the 2 dataframes into 1 super dataframe for species accumulation goodness
fastmerge <- function(d1, d2) {
  d1.names <- names(d1)
  d2.names <- names(d2)
  
  # columns in d1 but not in d2
  d2.add <- setdiff(d1.names, d2.names)
  
  # columns in d2 but not in d1
  d1.add <- setdiff(d2.names, d1.names)
  
  # add blank columns to d2
  if(length(d2.add) > 0) {
    for(i in 1:length(d2.add)) {
      d2[d2.add[i]] <- NA
    }
  }
  
  # add blank columns to d1
  if(length(d1.add) > 0) {
    for(i in 1:length(d1.add)) {
      d1[d1.add[i]] <- NA
    }
  }
  
  return(rbind(d1, d2))
}
	## Using a fastmerge function I found on stack overflow
fastmerge(trans_hosp6, env_merged_df) ->hosp6_merged_df
	## generating a merged abundance table for importing into vegan

write.csv(hosp6_merged_df,'hkbn_hosps_env_metro_abundances.csv')
	## writing a CSV


---

### Generation of Rarefaction/Species Accumulation Curve in VEGAN R

---

install.packages("permute", "lattice", "lattice", "MASS", "mgcv", "tcltk", "vegan", "readr") #installing packages

library("permute")  #libraries required as prerequisites to vegan
library("lattice")
library("MASS")
library("mgcv")
library("tcltk")
library("vegan")
library("readr")


hkbn_hosps_env_metro_abundances <- read_csv("C:/Users/jmanast.IASTATE/Downloads/hkbn_hosps_env_metro_abundances.csv")



hkbn_hosps_env_metro_abundances[is.na(hkbn_hosps_env_metro_abundances)]<-0 #changes all the empty cells labed as NA to the number 0, Vegan must have a non-negative number placed in each cell to run.
hkbn_hosps_env_metro_abundances <- hkbn_hosps_env_metro_abundances[,-1] #gitting rid of the first column, Vegan did not like our rownames

new_merged_df<-hosp6_merged_df[c("Hospital", "Hospital_6", "Metro", "Environment"), ]

spa <- specaccum(new_merged_df, "collector") #measure of species richness (y axis) and sample (y axis), collector tells vegan to read and plot in order data is in dataframe
plot(spa, main="Species Accumilation Curve/Rarefaction Curve", ylab ="Richness") #ploting spa

---



# Summary of Replication


Due to time and space limitations and lack of documentation, it was only possible for us to  make use of a subset of the data to generate relevant figures.  Two figures of note were the result of our work.  The first, a heat map representing abundance of microbes present in samples taken from particular locations, was generated using MetaPhlAn.






