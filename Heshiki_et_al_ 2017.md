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

### Generation of Rarefaction Curve

---

install.packages() #installing packages

library("permute")  #librarys required as prerequisites (spelt wrong) to vegan
library("lattice")
library("MASS")
library("mgcv")
library("tcltk")
library("vegan")
library("readr")

hkbn_abundances <- read_csv("C:/Users/lanas/Desktop/r for final project/R_project_final/abundance table/hkbn_abundances.csv") #
hkbn_abundances[is.na(hkbn_abundances)]<-0 #changes all the empty cells labed as NA to the number 0, Vegan must have a non-negative number placed in each cell to run.
hkbn_abundances<- hkbn_abundances[,-1] #gitting rid of the first column, Vegan did not like our



spa <- specaccum(hkbn_abundances) #messure of species richness (y axis) and sample (y axis)
plot(spa)

---



# Summary of Replication









