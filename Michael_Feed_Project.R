###PREPPING YOUR R ENVIRONMENT###

##Set your working directory##
#Use getwd() to find out where you currently are
getwd()
#Use setwd() to set your working directory 
setwd("C:/Users/dbrad/Documents/Bioinformatic_Analysis/Michael_Secret_Feed/R")
#Rerun getwd() to check it worked
getwd()

##Load up SummarySE function##

## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summarized
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

##Load packages##
library(phyloseq)#
library(vegan)#
library(ggplot2)#
library(plyr)#
library(tidyverse)#
library(FSA)#
library(dplyr)#
library(reshape)#
library(DESeq2)
library(metagenomeSeq)
library(pairwiseAdonis)#
library(rstatix)
library(ggvenn)
library(lefser)

###GETTING YOUR DATA INTO R###

##Import files##

BIOM <- import_biom(file.choose()) #Section_1_QIIME2/Test/phyloseq.biom
TREE =  read_tree(file.choose()) #Section_1_QIIME2/Test/tree/tree.nwk
META <- import_qiime_sample_data(file.choose()) #Section_1_QIIME2/test_map
META2  = read.delim(file.choose(), row.names=1)

#check that sample names are the same between metadata and abundance 
sample_names(META)
sample_names(BIOM)

#Merge three items into one phyloseq object called data#
data <- merge_phyloseq (BIOM,TREE,META)
data
###DATA MANIPULATION AND SUMMARIZING SEQUENCES### 

#See what present names are#
colnames(tax_table(data))

#Switch them to classic format#
colnames(tax_table(data))=c("Kingdom","Phylum","Class","Order","Family","Genus", "Species")

#Check that it worked#
colnames(tax_table(data))

##Subsetting and merging phyloseqs##

#Check number of samples
nsamples(data) #44#

##Filtering ASVs##

#Check number of taxa#
ntaxa(data) #2000#

#Filter low abundance ASVs

Fdata = filter_taxa(data, function(x) sum(x) >10, TRUE)
ntaxa(Fdata) #940#

#Copy out asv abundance (phyloseq uses otu, probably because it is just an aesthetics thing) table as a separate matrix from phyloseq object into R environment#
otutable = as(otu_table(Fdata), "matrix")

#Bind to it the taxonomy table which you copy out from phyloseq in the same step via the same rownames#
taxotutable = cbind(as(otutable, "matrix"), as(tax_table(data)[rownames(otutable), ], "matrix"))

#Export ASV & Taxonomy table out of R using the write function. Can write it as many things, but csv is my preferred method, cause easier to open in Excel#
write.csv(taxotutable , "test.otu.tax.csv")

##Summarize sequences by ASVs 

#Unfiltered version first#

#export table of sequence sums#
data_seqs_per_ESV <- as.data.frame(taxa_sums(data))

#Change colnames to Sequences#
colnames(data_seqs_per_ESV) <- c("Sequences")
sum(data_seqs_per_ESV$Sequences) #305935#

#Use cbind to get taxonomy added in#
data_seqs_per_ESV = cbind(as(data_seqs_per_ESV, "data.frame"), as(tax_table(data)[rownames(data_seqs_per_ESV), ], "matrix"))

#export it as a csv, and tell it that the the rownames should be named OTUID#
write.csv(data.frame("OTUID" =rownames(data_seqs_per_ESV), data_seqs_per_ESV) , "seqs_per_ESV.csv", row.names=FALSE)

#Filtered version second#

Fdata_seqs_per_ESV <- as.data.frame(taxa_sums(Fdata))
colnames(Fdata_seqs_per_ESV) <- c("Sequences")
sum(Fdata_seqs_per_ESV$Sequences) #302661#
Fdata_seqs_per_ESV = cbind(as(Fdata_seqs_per_ESV, "data.frame"), as(tax_table(Fdata)[rownames(Fdata_seqs_per_ESV), ], "matrix"))
write.csv(data.frame("OTUID" =rownames(Fdata_seqs_per_ESV), Fdata_seqs_per_ESV) , "filtered_seqs_per_ESV.csv", row.names=FALSE)

##Summarize sequences by samples (filtered and unfiltered)##

#Unfiltered first
data_seqs_per_sample <- as.data.frame(sample_sums(data))
colnames(data_seqs_per_sample) <- c("Full_Sequences")
sum(data_seqs_per_sample$Full_Sequences) #305935#
Fdata_seqs_per_sample <- as.data.frame(sample_sums(Fdata))
colnames(Fdata_seqs_per_sample) <- c("Trimmed_Sequences")
sum(Fdata_seqs_per_sample$Trimmed_Sequences) #302661#

#Sum of sequences should be same as above. Instead of exporting two different tables you can instead combine them with cbind and export that#
Comdata_seqs_per_sample = cbind(as(data_seqs_per_sample, "data.frame"), as(Fdata_seqs_per_sample, "data.frame"))
write.csv(data.frame("OTUID" =rownames(Comdata_seqs_per_sample), Comdata_seqs_per_sample) , "sequences_per_sample.csv", row.names=FALSE)


#All feed samples before and after filtering have <1000 sequences, thus they will be removed from further analysis

Fdata_Gut <- subset_samples(Fdata, Type=="Gut")
Fdata_Gut #32 samples, 940 taxa
Fdata_Gut = filter_taxa(Fdata_Gut, function(x) sum(x) >0, TRUE)
Fdata_Gut # 32 samples, 846 taxa

#Copy out asv abundance (phyloseq uses otu, probably because it is just an aesthetics thing) table as a separate matrix from phyloseq object into R environment#
Fdata_Gut_otutable = as(otu_table(Fdata_Gut), "matrix")

#Bind to it the taxonomy table which you copy out from phyloseq in the same step via the same rownames#
Fdata_Gut_tax_otu_table = as.data.frame(cbind(as(Fdata_Gut_otutable, "matrix"), as(tax_table(Fdata_Gut)[rownames(Fdata_Gut_otutable), ], "matrix")))


##Determine numbers of various taxonomic levels: total and defined

##INFOMATION USED FOR SUPPLEMENTARY TABLE 2

#Total numbers are all the unique ids given to organisms at that level while defined means that the id is actually biologically relevant instead of being uncultured, or unknown, or a repeat of the previous taxonomic level

#Extract out the taxonomy table from filtered dataset
Fdata_Guttaxa <- as.data.frame(tax_table(Fdata_Gut))

#Remove the labels at each taxonomic level to allow comparisons between taxonomic levels without it being a compounding variable
#For examples f__Unknown Family and g__Unknown Family would not be comparable otherwise
Fdata_Guttaxa <- data.frame(lapply(Fdata_Guttaxa, function(x) {gsub("d__", "", x)}))
Fdata_Guttaxa <- data.frame(lapply(Fdata_Guttaxa, function(x) {gsub("p__", "", x)}))
Fdata_Guttaxa <- data.frame(lapply(Fdata_Guttaxa, function(x) {gsub("c__", "", x)}))
Fdata_Guttaxa <- data.frame(lapply(Fdata_Guttaxa, function(x) {gsub("o__", "", x)}))
Fdata_Guttaxa <- data.frame(lapply(Fdata_Guttaxa, function(x) {gsub("f__", "", x)}))
Fdata_Guttaxa <- data.frame(lapply(Fdata_Guttaxa, function(x) {gsub("g__", "", x)}))
Fdata_Guttaxa <- data.frame(lapply(Fdata_Guttaxa, function(x) {gsub("s__", "", x)}))

#Kingdom
#Total
Fdata_Guttaxa%>%distinct(Kingdom, .keep_all=TRUE)%>%nrow #1
#Defined
Fdata_Guttaxa%>%distinct(Kingdom, .keep_all=TRUE)%>% filter(!is.na(Kingdom), Kingdom!="uncultured") %>%nrow #1

#Phylum
#Total
Fdata_Guttaxa%>%distinct(Kingdom, Phylum, .keep_all=TRUE)%>%nrow #12
#Defined
Fdata_Guttaxa%>%distinct(Kingdom, Phylum, .keep_all=TRUE)%>% filter(!is.na(Phylum), Phylum!="uncultured") %>%nrow #12

#Class
#Total
Fdata_Guttaxa%>%distinct(Kingdom, Phylum, Class, .keep_all=TRUE)%>%nrow #18
#Defined
Fdata_Guttaxa%>%distinct(Kingdom, Phylum, Class, .keep_all=TRUE)%>% filter(!is.na(Class), !grepl(pattern = "uncultured", x = Class)) %>% filter(!grepl(pattern = "unidentified", x = Class))%>% filter(!grepl(pattern = "metagenome", x = Class))%>%filter(as.character(Phylum) != as.character(Class))%>%nrow #18

#Order
#Total
Fdata_Guttaxa%>%distinct(Kingdom, Phylum, Class, Order, .keep_all=TRUE)%>%nrow #51
#Defined
Fdata_Guttaxa%>%distinct(Kingdom, Phylum, Class, Order, .keep_all=TRUE)%>% filter(!is.na(Order), !grepl(pattern = "uncultured", x = Order)) %>% filter(!grepl(pattern = "unidentified", x = Order))%>% filter(!grepl(pattern = "metagenome", x = Order))%>%filter(!grepl(pattern = "Unknown", x = Order))%>%filter(as.character(Class) != as.character(Order))%>%nrow #49

#Family
#Total
Fdata_Guttaxa%>%distinct(Kingdom, Phylum, Class, Order, Family, .keep_all=TRUE)%>%nrow #81
#Defined
Fdata_Guttaxa%>%distinct(Kingdom, Phylum, Class, Order, Family, .keep_all=TRUE)%>% filter(!is.na(Family), !grepl(pattern = "uncultured", x = Family)) %>% filter(!grepl(pattern = "unidentified", x = Family))%>% filter(!grepl(pattern = "metagenome", x = Family))%>%filter(!grepl(pattern = "Unknown", x = Family))%>%filter(as.character(Order) != as.character(Family))%>%nrow #70

#Genus
#Total
Fdata_Guttaxa%>%distinct(Kingdom, Phylum, Class, Order, Family, Genus, .keep_all=TRUE)%>%nrow #148
#Defined
Fdata_Guttaxa%>%distinct(Kingdom, Phylum, Class, Order, Family, Genus, .keep_all=TRUE)%>% filter(!is.na(Genus), !grepl(pattern = "uncultured", x = Genus)) %>% filter(!grepl(pattern = "unidentified", x = Genus))%>% filter(!grepl(pattern = "metagenome", x = Genus))%>%filter(!grepl(pattern = "Unknown", x = Genus))%>%filter(as.character(Family) != as.character(Genus))%>%nrow #114

#Species
#Total
Fdata_Guttaxa%>%distinct(Kingdom, Phylum, Class, Order, Family, Genus, Species, .keep_all=TRUE)%>%nrow #195
#Defined
Fdata_Guttaxa%>%distinct(Kingdom, Phylum, Class, Order, Family, Genus, Species, .keep_all=TRUE)%>% filter(!is.na(Species), !grepl(pattern = "uncultured", x = Species)) %>% filter(!grepl(pattern = "unidentified", x = Species))%>% filter(!grepl(pattern = "uncultured", x = Genus))%>% filter(!grepl(pattern = "metagenome", x = Species))%>%filter(!grepl(pattern = "Unknown", x = Species))%>%filter(as.character(Genus) != as.character(Species))%>% filter(!grepl(pattern = "bacterium", x = Species))%>% filter(!grepl(pattern = "enrichment", x = Species))%>% filter(!grepl(pattern = "sp.", x = Species))%>% filter(!grepl(pattern = "endosymbiont", x = Species))%>%nrow #20



###BETA DIVERSITY ANALYSIS OTHER TRANSFORMATIONS VERSION###

##Transform data using square root
SRFdata_Gut <- transform_sample_counts(Fdata_Gut, function(x){x^(1/2)})

##Transform data using DESeq2
#Alternative:
#https://github.com/joey711/phyloseq/issues/299
#https://github.com/joey711/phyloseq/issues/229

#Load alternative gm_mean function to handle zeros based upon github issue shown above#
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# Use `~1` as the experimental design so that the actual design doesn't influence your tranformation.
dds = phyloseq_to_deseq2(Fdata_Gut, ~1)

#Calculate geometric means prior to estimate size factors#
geoMeans = apply(counts(dds), 1, gm_mean)
dds = estimateSizeFactors(dds, geoMeans = geoMeans)

#Conduct DESEQ2 test#
dds = DESeq(dds, fitType="local")

#Make a copy of Fdata_Gut so you can have a designated DESeq transformed copy
Fdata_Gut
DSFdata_Gut = Fdata_Gut
DSFdata_Gut

#Switch the asv table with the DESeq2 transformed data
otu_table(DSFdata_Gut) <- otu_table(getVarianceStabilizedData(dds), taxa_are_rows = TRUE)

#Check to see if your basic phyloseq information was not changed
DSFdata_Gut

#https://bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis
#https://github.com/joey711/phyloseq/issues/445

#Make any negative values equal to 0 since the more negative they are the more likely they were to be zero over very small and unlikely to affect your results

#Make a copy of DSFdata_Gut to manipulate
ZDSFdata_Gut <- DSFdata_Gut
ZDSFdata_Gut

#extract out asv table
DESeq2_otu_table <- as.data.frame(otu_table(ZDSFdata_Gut))

#Change all negatives to zero
DESeq2_otu_table[DESeq2_otu_table < 0.0] <- 0.0

#Switch out the asv table in phyloseq object
otu_table(ZDSFdata_Gut) <-otu_table(DESeq2_otu_table, taxa_are_rows = TRUE)

#Check to make sure basic phyloseq info did not change
ZDSFdata_Gut

#Show how amount of positive numbers changed throughout the transformation 
z <- otu_table(Fdata_Gut)
table(as.vector(z) > 0) / prod(dim(z))

# FALSE       TRUE 
# 0.92723109 0.07276891  

z <- otu_table(DSFdata_Gut)
table(as.vector(z) > 0) / prod(dim(z))

#FALSE      TRUE 
#0.92723109 0.07276891 

z <- otu_table(ZDSFdata_Gut)
table(as.vector(z) > 0) / prod(dim(z))

#FALSE      TRUE 
#0.92723109 0.07276891 

##Standardize with CSS:
#https://forum.qiime2.org/t/css-normalization-of-sequence-counts-with-metagenomeseq-in-r/4288

#Convert file types to use in metagenomeSeq
MGS <- phyloseq_to_metagenomeSeq(Fdata_Gut) 

#Perform normalization following: https://bioconductor.org/packages/release/bioc/vignettes/metagenomeSeq/inst/doc/metagenomeSeq.pdf
CSSp <- cumNormStatFast(MGS)
CSSp
MGS <- cumNorm(MGS, p = CSSp)

#Store post-transformation abundance table
normmybiom <- MRcounts(MGS, norm = T)

#Create copy of filtered data to be modified
MGSFdata_Gut = Fdata_Gut
MGSFdata_Gut

#Switch out old ASV table for transformed one
otu_table(MGSFdata_Gut) <-otu_table(normmybiom, taxa_are_rows = TRUE)
MGSFdata_Gut

#Plot each of the different methods to see the differences#
plot_ordination(SRFdata_Gut, ordinate(SRFdata_Gut, "PCoA", "bray"), color = "Treatment..", shape = "Treatment..") 
plot_ordination(ZDSFdata_Gut, ordinate(ZDSFdata_Gut, "PCoA", "bray"), color = "Treatment..", shape = "Treatment..")
plot_ordination(MGSFdata_Gut, ordinate(MGSFdata_Gut, "PCoA", "bray"), color = "Treatment..", shape = "Treatment..")

#Check your library size differences by dividing the sequences sum of the sample with the most sequences by the one with the least sequences typically want it to be <10x difference


max(sample_sums(Fdata_Gut))/min(sample_sums(Fdata_Gut))#16.24108
max(sample_sums(SRFdata_Gut))/min(sample_sums(SRFdata_Gut)) #3.219232
max(sample_sums(ZDSFdata_Gut))/min(sample_sums(ZDSFdata_Gut)) #3.15008
max(sample_sums(MGSFdata_Gut))/min(sample_sums(MGSFdata_Gut)) #32.17254

#Check to see if they are statistically normal, if >0.05 then it is likely normal
shapiro.test(sample_sums(Fdata_Gut)) #2.757e-06
shapiro.test(sample_sums(SRFdata_Gut)) #0.175
shapiro.test(sample_sums(ZDSFdata_Gut)) #0.4932
shapiro.test(sample_sums(MGSFdata_Gut)) #3.267e-10

#Visual representation of the sample sums
hist(sample_sums(Fdata_Gut))
hist(sample_sums(SRFdata_Gut))
hist(sample_sums(ZDSFdata_Gut))
hist(sample_sums(MGSFdata_Gut))

#ZDSFdata_Gut performed the best, having lowest ratio, and more obvious normally distributed histogram

###BETA DIVERSITY ANALYSIS SQUARE ROOT VERSION###

##Export ESV table for PRIMER7##

# Extract abundance matrix from the phyloseq object#
QIIME2_WS_Bio = as(otu_table(Fdata_Gut), "matrix")

# Coerce to data.frame#
QIIME2_WS_Bio = as.data.frame(QIIME2_WS_Bio)

#Make a txt file#
write.table(QIIME2_WS_Bio, file = 'QIIME2_WS_Bio.txt', sep = "\t")

#Extract Distance Matrix for PRIMER7
Bray<-phyloseq::distance(ZDSFdata_Gut, "bray")
Bray_ZDS_df<-as.data.frame(as.matrix(Bray))
write.csv(Bray_ZDS_df, "morista.csv")


##PERMANOVA Analysis##

#make a distance matrix#
Bray<-phyloseq::distance(ZDSFdata_Gut, "bray")

#convert metadata of the phyloseq to a dataframe#
sampledf<-data.frame(sample_data(ZDSFdata_Gut))

#run a PERMANVOA using vegan's adonis function#
adonis2(Bray~Treatment.., data=sampledf, permutations=9999)

# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 9999
# 
# adonis2(formula = Bray ~ Treatment.., data = sampledf, permutations = 9999)
# Df SumOfSqs      R2      F Pr(>F)    
# Treatment..  3   1.4237 0.15001 1.6472  1e-04 ***
#   Residual    28   8.0668 0.84999                  
# Total       31   9.4905 1.00000                  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

pairwise.adonis(Bray, sampledf$Treatment.., perm = 9999)


##PCoA Analysis##

#Make a dissimilarity/similarity matrix, by default Bray Curtis is used#
ZDSFdata_GutBrayPCoAO <- ordinate(ZDSFdata_Gut,"PCoA")

PCoA_Bray_Treatment <- plot_ordination(ZDSFdata_Gut, ZDSFdata_GutBrayPCoAO, color="Treatment..", shape="Treatment..")+
  geom_point(size=7, alpha=0.75)+
  scale_colour_manual(name="Experimental Diets", values=c("#FFFFFF", "#92D050", "#ED7D31", "#0070C0"), labels=c("Control (100% FM)", "50% DMM", "100% DMM", "50% WMM")) +
  scale_shape_manual(name="Experimental Diets", values =c(15, 16, 17, 18), labels=c("Control (100% FM)", "50% DMM", "100% DMM", "50% WMM")) +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=11), legend.key.size = unit(0.7,"cm"), legend.title=element_text(size=13))
PCoA_Bray_Treatment

##Save Graph for Publication##
tiff('Figure 3 - PCoA by Treatment.tiff', units="in", width=10, height=6, res=300)
PCoA_Bray_Treatment
dev.off()

###TAXONOMIC BAR PLOTS###

# get abundance in relative percentage#
phy <- transform_sample_counts(Fdata_Gut, function(x) 100*x/sum(x))

##Create table ready for making stacked bar graph for Phyla >1% by Site-Survey combo##

# agglomerate taxa such that the ASVs with the same upper taxonomic ranks are merged at whatever rank you choose#
glom <- tax_glom(phy, taxrank = 'Phylum')

# Create dataframe from phyloseq object, just melts the phyloseq object to a full table#
dat <- psmelt(glom)

# convert Phylum to a character vector from a factor because R#
dat$Phylum <- as.character(dat$Phylum)

# group dataframe by Phylum, calculate mean rel. abundance#
means <- ddply(dat, ~Phylum, function(x) c(mean=mean(x$Abundance)))

# find Phyla whose rel. abund. is less than 1%#
Other <- means[means$mean <= 1,]$Phylum

# change their name to "Other Prokaryotes"#
dat[dat$Phylum %in% Other,]$Phylum <- 'Other Prokaryotes'

#remove all Phylums labeled Other Prokaryotes#
dat <-dat[!dat$Phylum=="Other Prokaryotes",]

#remove unnessary columns#
dat <- subset(dat, select=c(Treatment.., Abundance, Phylum))

#Summarize based upon target parameter#
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("Treatment..","Phylum"), na.rm = TRUE)

#Arrange by Abundance#
dat <- arrange(dat, Abundance)

#Create a table that is the leftover Prokaryotes#
Abundance <- ddply(dat, ~Treatment.., function(x) c(Abundance=100-sum(x$Abundance)))

#Add a column labeling the leftover Prokaryotes#
Abundance$Phylum<- "Other Prokaryotes"

#remove unnessary columns#
dat <- subset(dat, select=c(Treatment.., Abundance, Phylum))

#combine with original table#
Phylum_Fdata_Gut_Treatment <- rbind(dat, Abundance)

##Create table ready for making stacked bar graph for Genuses <0.5%##
# agglomerate taxa
glom <- tax_glom(phy, taxrank = 'Genus')
# create dataframe from phyloseq object
dat <- psmelt(glom)
#Combine Family and Genus columns#
dat<- unite(dat, Genus, Family:Genus, sep=';')
# convert Genus to a character vector from a factor because R
dat$Genus <- as.character(dat$Genus)
# group dataframe by Genus, calculate mean rel. abundance
means <- ddply(dat, ~Genus, function(x) c(mean=mean(x$Abundance)))
# find Phyla whose rel. abund. is less than 1%
Other <- means[means$mean <= 1,]$Genus
# change their name to "Other Prokaryotes"
dat[dat$Genus %in% Other,]$Genus <- 'Other Prokaryotes'
#remove all Genera labeled Other Prokaryotes
dat <-dat[!dat$Genus=="Other Prokaryotes",]
#remove unnessary columns
dat <- subset(dat, select=c(Treatment.., Abundance, Genus))
#Summarize based upon target parameter
dat <- summarySE(data=dat, measurevar="Abundance", groupvars=c("Treatment..","Genus"), na.rm = TRUE)
#Arrange by Abundance
dat <- arrange(dat, Abundance)
#Create a table that is the leftover Prokaryotes
Abundance <- ddply(dat, ~Treatment.., function(x) c(Abundance=100-sum(x$Abundance)))
#Add a column labeling the leftover Prokaryotes
Abundance$Genus<- "Other Prokaryotes"
#remove unnessary columns
dat <- subset(dat, select=c(Treatment.., Abundance, Genus))
#combine with original table
Genus_Fdata_Gut_Treatment <- rbind(dat, Abundance)


##Make a series of graphs using ggplot based on above tables##

#Make Treatment Phylum graph#

spatial_plot_Phylum_Fdata_Gut_Treatment <- ggplot(data=Phylum_Fdata_Gut_Treatment, aes(x=Treatment.., y=Abundance, fill=Phylum)) + 
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("Phylum", values=c("orange", "lightblue", "purple", "firebrick3"))+
  scale_x_discrete(name="Experimental Diets", labels=c("Control (100% FM)", "50% DMM", "100% DMM", "50% WMM")) +
  xlab("Treatment") +
  ylab("Percentage")+
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=11), legend.key.size = unit(0.7,"cm"), legend.title=element_text(size=13))
  #ggtitle("Phyla by Treatment")+
  #theme(plot.title = element_text(hjust = 0.5)) +
  #scale_x_discrete(limits=c("FP-Pre-Hurricane", "FP-Post-Hurricane", "HT-Pre-Hurricane", "HT-Post-Hurricane"), labels=c("Fort Pierce before Irma", "Fort Pierce after Irma", "Harbortown Marina before Irma", "Harbortown Marina after Irma"))
spatial_plot_Phylum_Fdata_Gut_Treatment

##Save Graph for Publication##
tiff('Supplementary Figure 1 - Phylum by Treatment.tiff', units="in", width=10, height=6, res=300)
spatial_plot_Phylum_Fdata_Gut_Treatment
dev.off()

##Make Genus Graph##

spatial_plot_Genus_Fdata_Gut_Treatment <- ggplot(data=Genus_Fdata_Gut_Treatment, aes(x=Treatment.., y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("Family;Genus", values=c("deeppink", "midnightblue", "blue", "red", "yellow", "orange", "green", "maroon", "dodgerblue", "purple", "lightblue", "orangered", "firebrick3", "darkseagreen3", "navy", "black", "seagreen", "turquoise", "white", "rosybrown"))+
  #ggtitle("Genera by Treatment") +
  scale_x_discrete(name="Experimental Diets", labels=c("Control (100% FM)", "50% DMM", "100% DMM", "50% WMM")) +
  xlab("Treatment") +
  ylab("Percentage") +
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=11), legend.key.size = unit(0.7,"cm"), legend.title=element_text(size=13))
  #theme(plot.title = element_text(hjust = 0.5)) +
  #scale_x_discrete(limits=c("FP-Pre-Hurricane", "FP-Post-Hurricane", "HT-Pre-Hurricane", "HT-Post-Hurricane"), labels=c("Fort Pierce before Irma", "Fort Pierce after Irma", "Harbortown Marina before Irma", "Harbortown Marina after Irma"))
spatial_plot_Genus_Fdata_Gut_Treatment

##Save Graph for Publication##
tiff('Figure 4 - Family;Genera by Treatment.tiff', units="in", width=10, height=10, res=300)
spatial_plot_Genus_Fdata_Gut_Treatment
dev.off()

#You will have to probably adjust the plot area 

#Determine Genus Fold Changes

#Switch the structure of dataframe from long to wide
wide_Genus_Fdata_Gut_Treatment <- spread(Genus_Fdata_Gut_Treatment, Treatment.., Abundance) 
wide_Genus_Fdata_Gut_Treatment

#Add a new column that is the 10 ppt abundance of a particular genera divided by the 20 ppt abundance of that same genera
wide_Genus_Fdata_Gut_Treatment$TM1TM4 <- wide_Genus_Fdata_Gut_Treatment$TM1 / wide_Genus_Fdata_Gut_Treatment$TM4 

#Same thing except that its TM4 divided by TM1
wide_Genus_Fdata_Gut_Treatment$TM4TM1 <- wide_Genus_Fdata_Gut_Treatment$TM4 / wide_Genus_Fdata_Gut_Treatment$TM1

#Same thing except that its TM4 divided by TM1
wide_Genus_Fdata_Gut_Treatment$TM4TM1 <- wide_Genus_Fdata_Gut_Treatment$TM4 / wide_Genus_Fdata_Gut_Treatment$TM1

#TM1 divided by TM2
wide_Genus_Fdata_Gut_Treatment$TM1TM2 <- wide_Genus_Fdata_Gut_Treatment$TM1 / wide_Genus_Fdata_Gut_Treatment$TM2 

#TM2 divided by TM1
wide_Genus_Fdata_Gut_Treatment$TM2TM1 <- wide_Genus_Fdata_Gut_Treatment$TM2 / wide_Genus_Fdata_Gut_Treatment$TM1

#TM1 divided by TM3
wide_Genus_Fdata_Gut_Treatment$TM1TM3 <- wide_Genus_Fdata_Gut_Treatment$TM1 / wide_Genus_Fdata_Gut_Treatment$TM3 

#TM3 divided by TM1
wide_Genus_Fdata_Gut_Treatment$TM3TM1 <- wide_Genus_Fdata_Gut_Treatment$TM3 / wide_Genus_Fdata_Gut_Treatment$TM1


###ALPHA DIVERSITY###


##Calculate alpha diversity##

Falphadiv = estimate_richness(Fdata_Gut, split = TRUE)
rownames(Falphadiv) <- c("1-1FE",  "1-2FE",  "10-1FE", "10-2FE", "11-1FE", "11-2FE", "12-1FE", "12-2FE", "13-1FE", "13-2FE", "14-1FE", "14-2FE", "15-1FE", "15-2FE", "16-1FE", "16-2FE", "2-1FE",  "2-2FE",  "3-1FE",  "3-2FE", "4-1FE",  "4-2FE",  "5-1FE",  "5-2FE",  "6-1FE",  "6-2FE",  "7-1FE",  "7-2FE",  "8-1FE",  "8-2FE", "9-1FE",  "9-2FE")
Falphadiv
write.csv(Falphadiv, file='Falphadiv.csv')

#cbind in the metadata from the filtered phyloseq#
Falphadiv_metadata = cbind(as(Falphadiv, "data.frame"), as(sample_data(Fdata_Gut)[rownames(Falphadiv), ], "data.frame"))

##Test for normality##
shapiro.test(Falphadiv_metadata$Shannon) #2.049e-07
shapiro.test(Falphadiv_metadata$Observed) #0.4953
shapiro.test(Falphadiv_metadata$Fisher) #0.1027
shapiro.test(Falphadiv_metadata$Simpson) #1.804e-09

#Shannon and Simpson are not normally distributed thus must use Spearman for correlation tests

##Correlation Tests##
cor.test(Falphadiv_metadata$Shannon, Falphadiv_metadata$Observed, method = "spearman", exact = FALSE) #1.582e-06
cor.test(Falphadiv_metadata$Shannon, Falphadiv_metadata$Fisher, method = "spearman", exact = FALSE) #5.722e-08
cor.test(Falphadiv_metadata$Shannon, Falphadiv_metadata$Simpson, method = "spearman", exact = FALSE) #<2.2e-16

##Multiple test adjustment##
#Make a list of all the p vales from the correlation tests, make it a dataframe, add column of adjusted p values#
alpha_div_p.value=c(1.582e-06,5.722e-08,2.2e-16)
alpha_div_p.value=data.frame(alpha_div_p.value)
alpha_div_p.value$padj <- p.adjust(alpha_div_p.value$alpha_div_p.value, method = "BH")
alpha_div_p.value

# alpha_div_p.value      padj
# 1         1.582e-06 1.582e-06
# 2         5.722e-08 8.583e-08
# 3         2.200e-16 6.600e-16

#Shannon is significantly correlated to other alpha diversity metrics

#Besides normality, another good thing to test for statistical analysis is heteroscedacity
bartlett.test(Shannon ~ Treatment.., Falphadiv_metadata) #0.0006905

#Shannon is non-normal and heteroscedastic, will use Welch's ANOVA and Games-Howell to test

##Basic visualization of Shannon statistics##
ddply(Falphadiv_metadata, .(Treatment..), summarise, mean=mean(Shannon), sd=sd(Shannon), median=median(Shannon), IQR=IQR(Shannon))


ddply(Falphadiv_metadata, .(Treatment..), summarise, mean=mean(Shannon), sd=sd(Shannon), median=median(Shannon), IQR=IQR(Shannon))

# Treatment..     mean        sd   median       IQR
# 1         TM1 3.416990 0.7986628 3.677094 0.3904523
# 2         TM2 2.984462 1.0495770 3.331221 0.8505816
# 3         TM3 2.677683 1.4891391 3.430333 1.0766987
# 4         TM4 3.512199 0.2211072 3.539258 0.1943155

#Make a simple boxplot to get a visual summary of these statistics#
ggplot(Falphadiv_metadata, aes(x=Treatment.., y=Shannon, color=Treatment..)) + 
  geom_boxplot() +
  theme(axis.text.x  = element_text(angle=45, hjust=1))

###TESTS FOR TOTAL AND PAIRWISE SIGNIFICANCE###
#Welch's ANOVA
oneway.test(Shannon ~ Treatment.., data = Falphadiv_metadata, var.equal = FALSE)

# One-way analysis of means (not assuming equal variances)
# 
# data:  Shannon and Treatment..
# F = 1.2974, num df = 3.000, denom df = 12.628, p-value = 0.3182

#Alpha diversity is not significant, so will not test pair-wise Treatments
games_howell_test(Falphadiv_metadata, Shannon ~ Treatment..)

###COMBINING STATISTICAL AND ALPHA DIVERSITY TESTING###

##Making a more informative boxplot##

#Most of this you've seen before with bar plot#
Shannon_Treatment_Boxplots <- ggplot(Falphadiv_metadata, aes(x=Treatment.., y=Shannon)) + 
  geom_boxplot() +
  theme(axis.text.x  = element_text(angle=45, hjust=1)) + 
  xlab("Treatment") +
  scale_x_discrete(name="Experimental Diets", labels=c("Control (100% FM)", "50% DMM", "100% DMM", "50% WMM")) +
  annotate("text", x = c(1:4) , y = c(4.1, 3.9, 3.8, 4.0), label = c("a", "a", "a", "a"), size=5)+
    theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=11), legend.key.size = unit(0.7,"cm"), legend.title=element_text(size=13))
Shannon_Treatment_Boxplots

tiff('Figure 2 - Shannon by Treatment.tiff', units="in", width=10, height=6, res=300)
Shannon_Treatment_Boxplots
dev.off()


mean_se(Falphadiv_metadata$Shannon) #3.147833
sd(Falphadiv_metadata$Shannon) #1.010666







###Venn Diagrams###
##Generate Sample Type based Venn Diagrams
Fdata_Feed <- subset_samples(Fdata, Type=="Feed")
Fdata_Feed #960
Fdata_Feed = filter_taxa(Fdata_Feed, function(x) sum(x) >0, TRUE)
Fdata_Feed #164

#Extract out ASVs as a list from the asv table stored in phyloseq objects
feed_asvs <- as.list(rownames(otu_table(Fdata_Feed)))
gut_asvs <- as.list(rownames(otu_table(Fdata_Gut)))

#Create a table agglomerated to genus level
Fdatataxa_Feed_genus <- Fdata_Feedtaxa%>%distinct(Kingdom, Phylum, Class, Order, Family, Genus, .keep_all=TRUE)

#Fuse the entire taxonomy together into one column, cannot just compare the genus level because multiple occurrences of "same" genus name, like g__uncultured, with different preceding taxonomy strings
Fdatataxa_Feed_genus$Taxonomy <- paste(Fdatataxa_Feed_genus$Kingdom,Fdatataxa_Feed_genus$Phylum,Fdatataxa_Feed_genus$Class,Fdatataxa_Feed_genus$Order,Fdatataxa_Feed_genus$Family,Fdatataxa_Feed_genus$Genus,sep="-")

#Extract out the full taxonomy list column as a list
feed_genera_list <- as.list(Fdatataxa_Feed_genus$Taxonomy)

#Do it for the tissue
Fdatataxa_Gut_genus <- Fdata_Guttaxa%>%distinct(Kingdom, Phylum, Class, Order, Family, Genus, .keep_all=TRUE)
Fdatataxa_Gut_genus$Taxonomy <- paste(Fdatataxa_Gut_genus$Kingdom,Fdatataxa_Gut_genus$Phylum,Fdatataxa_Gut_genus$Class,Fdatataxa_Gut_genus$Order,Fdatataxa_Gut_genus$Family,Fdatataxa_Gut_genus$Genus,sep="-")
gut_genera_list <- as.list(Fdatataxa_Gut_genus$Taxonomy)

#Combine the list into a list of lists
ASV_Genera_Lists <- list('Feed (164)' = feed_asvs,
                         'Fecal (846)' = gut_asvs,
                         'Feed (60)' = feed_genera_list,
                         'Fecal (148)' = gut_genera_list)

#Create a Venn Diagram comparing ASVs and save it
ggvenn(ASV_Genera_Lists, c("Feed (164)", "Fecal (846)"), fill_color = c("lightblue", "gray"), show_percentage = FALSE)+
  ggtitle("ASVs by Sample Type") +
  theme(plot.title = element_text(hjust = 0.5))

tiff('ASV level venn diagram.tiff', units="in", width=10, height=6, res=300)
#plot
dev.off()

#Create a Venn Diagram comparing Genera and save it
ggvenn(ASV_Genera_Lists, c("Feed (60)", "Fecal (148)"), fill_color = c("lightblue", "gray"), show_percentage = FALSE)+
  ggtitle("Genera by Sample Type") +
  theme(plot.title = element_text(hjust = 0.5))

tiff('Genus level venn diagram.tiff', units="in", width=10, height=6, res=300)
#plot
dev.off()

##Generate Treatment based Venn Diagrams for Gut Samples
###Venn Diagrams###

##TM1
#Extract out ASVs as a list from the asv table stored in phyloseq objects
TM1_Gut_Fdata <- subset_samples(Fdata_Gut, Treatment.. == "TM1")
TM1_Gut_Fdata #846
TM1_Gut_Fdata <- filter_taxa(TM1_Gut_Fdata, function(x) sum(x) >0, TRUE)
TM1_Gut_Fdata #333

TM1_Fdata_Gut_taxa <- as.data.frame(tax_table(TM1_Gut_Fdata))

#Create a table agglomerated to genus level
TM1_Fdata_Gut_taxa_genus <- TM1_Fdata_Gut_taxa%>%distinct(Kingdom, Phylum, Class, Order, Family, Genus, .keep_all=TRUE)
TM1_Fdata_Gut_taxa_genus <- subset(TM1_Fdata_Gut_taxa_genus, select=-c(Species))


#Fuse the entire taxonomy together into one column, cannot just compare the genus level because multiple occurrences of "same" genus name, like g__uncultured, with different preceding taxonomy strings
TM1_Fdata_Gut_taxa_genus$Taxonomy <- paste(TM1_Fdata_Gut_taxa_genus$Kingdom,TM1_Fdata_Gut_taxa_genus$Phylum,TM1_Fdata_Gut_taxa_genus$Class,TM1_Fdata_Gut_taxa_genus$Order,TM1_Fdata_Gut_taxa_genus$Family,TM1_Fdata_Gut_taxa_genus$Genus,sep="-")

#Extract out the full taxonomy list column as a list
TM1_genera_list <- as.list(TM1_Fdata_Gut_taxa_genus$Taxonomy) #92

##TM2
#Extract out ASVs as a list from the asv table stored in phyloseq objects
TM2_Gut_Fdata <- subset_samples(Fdata_Gut, Treatment.. == "TM2")
TM2_Gut_Fdata #846
TM2_Gut_Fdata <- filter_taxa(TM2_Gut_Fdata, function(x) sum(x) >0, TRUE)
TM2_Gut_Fdata #316

TM2_Fdata_Gut_taxa <- as.data.frame(tax_table(TM2_Gut_Fdata))

#Create a table agglomerated to genus level
TM2_Fdata_Gut_taxa_genus <- TM2_Fdata_Gut_taxa%>%distinct(Kingdom, Phylum, Class, Order, Family, Genus, .keep_all=TRUE)
TM2_Fdata_Gut_taxa_genus <- subset(TM2_Fdata_Gut_taxa_genus, select=-c(Species))

#Fuse the entire taxonomy together into one column, cannot just compare the genus level because multiple occurrences of "same" genus name, like g__uncultured, with different preceding taxonomy strings
TM2_Fdata_Gut_taxa_genus$Taxonomy <- paste(TM2_Fdata_Gut_taxa_genus$Kingdom,TM2_Fdata_Gut_taxa_genus$Phylum,TM2_Fdata_Gut_taxa_genus$Class,TM2_Fdata_Gut_taxa_genus$Order,TM2_Fdata_Gut_taxa_genus$Family,TM2_Fdata_Gut_taxa_genus$Genus,sep="-")

#Extract out the full taxonomy list column as a list
TM2_genera_list <- as.list(TM2_Fdata_Gut_taxa_genus$Taxonomy) #97

##TM3
#Extract out ASVs as a list from the asv table stored in phyloseq objects
TM3_Gut_Fdata <- subset_samples(Fdata_Gut, Treatment.. == "TM3")
TM3_Gut_Fdata #846
TM3_Gut_Fdata <- filter_taxa(TM3_Gut_Fdata, function(x) sum(x) >0, TRUE)
TM3_Gut_Fdata #262

TM3_Fdata_Gut_taxa <- as.data.frame(tax_table(TM3_Gut_Fdata))

#Create a table agglomerated to genus level
TM3_Fdata_Gut_taxa_genus <- TM3_Fdata_Gut_taxa%>%distinct(Kingdom, Phylum, Class, Order, Family, Genus, .keep_all=TRUE)
TM3_Fdata_Gut_taxa_genus <- subset(TM3_Fdata_Gut_taxa_genus, select=-c(Species))

#Fuse the entire taxonomy together into one column, cannot just compare the genus level because multiple occurrences of "same" genus name, like g__uncultured, with different preceding taxonomy strings
TM3_Fdata_Gut_taxa_genus$Taxonomy <- paste(TM3_Fdata_Gut_taxa_genus$Kingdom,TM3_Fdata_Gut_taxa_genus$Phylum,TM3_Fdata_Gut_taxa_genus$Class,TM3_Fdata_Gut_taxa_genus$Order,TM3_Fdata_Gut_taxa_genus$Family,TM3_Fdata_Gut_taxa_genus$Genus,sep="-")

#Extract out the full taxonomy list column as a list
TM3_genera_list <- as.list(TM3_Fdata_Gut_taxa_genus$Taxonomy) #83

##TM4
#Extract out ASVs as a list from the asv table stored in phyloseq objects
TM4_Gut_Fdata <- subset_samples(Fdata_Gut, Treatment.. == "TM4")
TM4_Gut_Fdata #846
TM4_Gut_Fdata <- filter_taxa(TM4_Gut_Fdata, function(x) sum(x) >0, TRUE)
TM4_Gut_Fdata #302

TM4_Fdata_Gut_taxa <- as.data.frame(tax_table(TM4_Gut_Fdata))

#Create a table agglomerated to genus level
TM4_Fdata_Gut_taxa_genus <- TM4_Fdata_Gut_taxa%>%distinct(Kingdom, Phylum, Class, Order, Family, Genus, .keep_all=TRUE)
TM4_Fdata_Gut_taxa_genus <- subset(TM4_Fdata_Gut_taxa_genus, select=-c(Species))

#Fuse the entire taxonomy together into one column, cannot just compare the genus level because multiple occurrences of "same" genus name, like g__uncultured, with different preceding taxonomy strings
TM4_Fdata_Gut_taxa_genus$Taxonomy <- paste(TM4_Fdata_Gut_taxa_genus$Kingdom,TM4_Fdata_Gut_taxa_genus$Phylum,TM4_Fdata_Gut_taxa_genus$Class,TM4_Fdata_Gut_taxa_genus$Order,TM4_Fdata_Gut_taxa_genus$Family,TM4_Fdata_Gut_taxa_genus$Genus,sep="-")

#Extract out the full taxonomy list column as a list
TM4_genera_list <- as.list(TM4_Fdata_Gut_taxa_genus$Taxonomy) #88




##Sample Type by Location
#Combine the list into a list of lists
Genera_Lists <- list('50% DMM (97)' = TM2_genera_list,
                     '100% DMM (83)' = TM3_genera_list,
                     '50% WMM (88)' = TM4_genera_list,
                     '100% FM (92)' = TM1_genera_list)

trmt_venn <- venn(Genera_Lists, ilab=TRUE, zcolor = "style")

ggvenn(Genera_Lists, c("100% FM (92)", "50% DMM (97)", "100% DMM (83)", "50% WMM (88)"),  show_percentage = FALSE)+
  ggtitle("Genera by Treatment") +
  theme(plot.title = element_text(hjust = 0.5))


tiff('Genera by Treatment Venn.tiff', units="in", width=10, height=6, res=300)
#graph
dev.off()

TM1_exclusives <- anti_join(TM1_Fdata_Gut_taxa_genus, TM2_Fdata_Gut_taxa_genus, by="Taxonomy")
TM1_exclusives <- anti_join(TM1_exclusives, TM3_Fdata_Gut_taxa_genus, by="Taxonomy")
TM1_exclusives <- anti_join(TM1_exclusives, TM4_Fdata_Gut_taxa_genus, by="Taxonomy")


TM2_exclusives <- anti_join(TM2_Fdata_Gut_taxa_genus, TM1_Fdata_Gut_taxa_genus, by="Taxonomy")
TM2_exclusives <- anti_join(TM2_exclusives, TM3_Fdata_Gut_taxa_genus, by="Taxonomy")
TM2_exclusives <- anti_join(TM2_exclusives, TM4_Fdata_Gut_taxa_genus, by="Taxonomy")

TM3_exclusives <- anti_join(TM3_Fdata_Gut_taxa_genus, TM1_Fdata_Gut_taxa_genus, by="Taxonomy")
TM3_exclusives <- anti_join(TM3_exclusives, TM2_Fdata_Gut_taxa_genus, by="Taxonomy")
TM3_exclusives <- anti_join(TM3_exclusives, TM4_Fdata_Gut_taxa_genus, by="Taxonomy")

TM4_exclusives <- anti_join(TM4_Fdata_Gut_taxa_genus, TM1_Fdata_Gut_taxa_genus, by="Taxonomy")
TM4_exclusives <- anti_join(TM4_exclusives, TM2_Fdata_Gut_taxa_genus, by="Taxonomy")
TM4_exclusives <- anti_join(TM4_exclusives, TM3_Fdata_Gut_taxa_genus, by="Taxonomy")


