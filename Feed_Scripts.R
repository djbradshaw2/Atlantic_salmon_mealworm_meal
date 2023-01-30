
##Determine numbers of various taxonomic levels: total and defined

##INFOMATION USED FOR SUPPLEMENTARY TABLE 2

#Total numbers are all the unique ids given to organisms at that level while defined means that the id is actually biologically relevant instead of being uncultured, or unknown, or a repeat of the previous taxonomic level

#Extract out the taxonomy table from filtered dataset
Fdata_Feedtaxa <- as.data.frame(tax_table(Fdata_Feed))

#Remove the labels at each taxonomic level to allow comparisons between taxonomic levels without it being a compounding variable
#For examples f__Unknown Family and g__Unknown Family would not be comparable otherwise
Fdata_Feedtaxa <- data.frame(lapply(Fdata_Feedtaxa, function(x) {gsub("d__", "", x)}))
Fdata_Feedtaxa <- data.frame(lapply(Fdata_Feedtaxa, function(x) {gsub("p__", "", x)}))
Fdata_Feedtaxa <- data.frame(lapply(Fdata_Feedtaxa, function(x) {gsub("c__", "", x)}))
Fdata_Feedtaxa <- data.frame(lapply(Fdata_Feedtaxa, function(x) {gsub("o__", "", x)}))
Fdata_Feedtaxa <- data.frame(lapply(Fdata_Feedtaxa, function(x) {gsub("f__", "", x)}))
Fdata_Feedtaxa <- data.frame(lapply(Fdata_Feedtaxa, function(x) {gsub("g__", "", x)}))
Fdata_Feedtaxa <- data.frame(lapply(Fdata_Feedtaxa, function(x) {gsub("s__", "", x)}))

#Kingdom
#Total
Fdata_Feedtaxa%>%distinct(Kingdom, .keep_all=TRUE)%>%nrow #1
#Defined
Fdata_Feedtaxa%>%distinct(Kingdom, .keep_all=TRUE)%>% filter(!is.na(Kingdom), Kingdom!="uncultured") %>%nrow #1

#Phylum
#Total
Fdata_Feedtaxa%>%distinct(Kingdom, Phylum, .keep_all=TRUE)%>%nrow #7
#Defined
Fdata_Feedtaxa%>%distinct(Kingdom, Phylum, .keep_all=TRUE)%>% filter(!is.na(Phylum), Phylum!="uncultured") %>%nrow #7

#Class
#Total
Fdata_Feedtaxa%>%distinct(Kingdom, Phylum, Class, .keep_all=TRUE)%>%nrow #9
#Defined
Fdata_Feedtaxa%>%distinct(Kingdom, Phylum, Class, .keep_all=TRUE)%>% filter(!is.na(Class), !grepl(pattern = "uncultured", x = Class)) %>% filter(!grepl(pattern = "unidentified", x = Class))%>% filter(!grepl(pattern = "metagenome", x = Class))%>%filter(as.character(Phylum) != as.character(Class))%>%nrow #9

#Order
#Total
Fdata_Feedtaxa%>%distinct(Kingdom, Phylum, Class, Order, .keep_all=TRUE)%>%nrow #23
#Defined
Fdata_Feedtaxa%>%distinct(Kingdom, Phylum, Class, Order, .keep_all=TRUE)%>% filter(!is.na(Order), !grepl(pattern = "uncultured", x = Order)) %>% filter(!grepl(pattern = "unidentified", x = Order))%>% filter(!grepl(pattern = "metagenome", x = Order))%>%filter(!grepl(pattern = "Unknown", x = Order))%>%filter(as.character(Class) != as.character(Order))%>%nrow #23

#Family
#Total
Fdata_Feedtaxa%>%distinct(Kingdom, Phylum, Class, Order, Family, .keep_all=TRUE)%>%nrow #36
#Defined
Fdata_Feedtaxa%>%distinct(Kingdom, Phylum, Class, Order, Family, .keep_all=TRUE)%>% filter(!is.na(Family), !grepl(pattern = "uncultured", x = Family)) %>% filter(!grepl(pattern = "unidentified", x = Family))%>% filter(!grepl(pattern = "metagenome", x = Family))%>%filter(!grepl(pattern = "Unknown", x = Family))%>%filter(as.character(Order) != as.character(Family))%>%nrow #35

#Genus
#Total
Fdata_Feedtaxa%>%distinct(Kingdom, Phylum, Class, Order, Family, Genus, .keep_all=TRUE)%>%nrow #60
#Defined
Fdata_Feedtaxa%>%distinct(Kingdom, Phylum, Class, Order, Family, Genus, .keep_all=TRUE)%>% filter(!is.na(Genus), !grepl(pattern = "uncultured", x = Genus)) %>% filter(!grepl(pattern = "unidentified", x = Genus))%>% filter(!grepl(pattern = "metagenome", x = Genus))%>%filter(!grepl(pattern = "Unknown", x = Genus))%>%filter(as.character(Family) != as.character(Genus))%>%nrow #52

#Species
#Total
Fdata_Feedtaxa%>%distinct(Kingdom, Phylum, Class, Order, Family, Genus, Species, .keep_all=TRUE)%>%nrow #76
#Defined
Fdata_Feedtaxa%>%distinct(Kingdom, Phylum, Class, Order, Family, Genus, Species, .keep_all=TRUE)%>% filter(!is.na(Species), !grepl(pattern = "uncultured", x = Species)) %>% filter(!grepl(pattern = "unidentified", x = Species))%>% filter(!grepl(pattern = "uncultured", x = Genus))%>% filter(!grepl(pattern = "metagenome", x = Species))%>%filter(!grepl(pattern = "Unknown", x = Species))%>%filter(as.character(Genus) != as.character(Species))%>% filter(!grepl(pattern = "bacterium", x = Species))%>% filter(!grepl(pattern = "enrichment", x = Species))%>% filter(!grepl(pattern = "sp.", x = Species))%>% filter(!grepl(pattern = "endosymbiont", x = Species))%>%nrow #14





###BETA DIVERSITY ANALYSIS OTHER TRANSFORMATIONS VERSION###

##Transform data using square root
SRFdata_Feed <- transform_sample_counts(Fdata_Feed, function(x){x^(1/2)})

##Transform data using DESeq2
#Alternative:
#https://github.com/joey711/phyloseq/issues/299
#https://github.com/joey711/phyloseq/issues/229

#Load alternative gm_mean function to handle zeros based upon github issue shown above#
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# Use `~1` as the experimental design so that the actual design doesn't influence your tranformation.
dds = phyloseq_to_deseq2(Fdata_Feed, ~1)

#Calculate geometric means prior to estimate size factors#
geoMeans = apply(counts(dds), 1, gm_mean)
dds = estimateSizeFactors(dds, geoMeans = geoMeans)

#Conduct DESEQ2 test#
dds = DESeq(dds, fitType="local")

#Make a copy of Fdata_Feed so you can have a designated DESeq transformed copy
Fdata_Feed
DSFdata_Feed = Fdata_Feed
DSFdata_Feed

#Switch the asv table with the DESeq2 transformed data
otu_table(DSFdata_Feed) <- otu_table(getVarianceStabilizedData(dds), taxa_are_rows = TRUE)

#Check to see if your basic phyloseq information was not changed
DSFdata_Feed

#https://bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#should-i-normalize-my-data-before-alpha-diversity-analysis
#https://github.com/joey711/phyloseq/issues/445

#Make any negative values equal to 0 since the more negative they are the more likely they were to be zero over very small and unlikely to affect your results

#Make a copy of DSFdata_Feed to manipulate
ZDSFdata_Feed <- DSFdata_Feed
ZDSFdata_Feed

#extract out asv table
DESeq2_otu_table <- as.data.frame(otu_table(ZDSFdata_Feed))

#Change all negatives to zero
DESeq2_otu_table[DESeq2_otu_table < 0.0] <- 0.0

#Switch out the asv table in phyloseq object
otu_table(ZDSFdata_Feed) <-otu_table(DESeq2_otu_table, taxa_are_rows = TRUE)

#Check to make sure basic phyloseq info did not change
ZDSFdata_Feed

#Show how amount of positive numbers changed throughout the transformation 
z <- otu_table(Fdata_Feed)
table(as.vector(z) > 0) / prod(dim(z))

# FALSE       TRUE 
# 0.8180894 0.1819106  

z <- otu_table(DSFdata_Feed)
table(as.vector(z) > 0) / prod(dim(z))

#FALSE      TRUE 
#0.8180894 0.1819106

z <- otu_table(ZDSFdata_Feed)
table(as.vector(z) > 0) / prod(dim(z))

#FALSE      TRUE 
#0.8180894 0.1819106

##Standardize with CSS:
#https://forum.qiime2.org/t/css-normalization-of-sequence-counts-with-metagenomeseq-in-r/4288

#Convert file types to use in metagenomeSeq
MGS <- phyloseq_to_metagenomeSeq(Fdata_Feed) 

#Perform normalization following: https://bioconductor.org/packages/release/bioc/vignettes/metagenomeSeq/inst/doc/metagenomeSeq.pdf
CSSp <- cumNormStatFast(MGS)
CSSp
MGS <- cumNorm(MGS, p = CSSp)

#Store post-transformation abundance table
normmybiom <- MRcounts(MGS, norm = T)

#Create copy of filtered data to be modified
MGSFdata_Feed = Fdata_Feed
MGSFdata_Feed

#Switch out old ASV table for transformed one
otu_table(MGSFdata_Feed) <-otu_table(normmybiom, taxa_are_rows = TRUE)
MGSFdata_Feed

#Plot each of the different methods to see the differences#
plot_ordination(SRFdata_Feed, ordinate(SRFdata_Feed, "PCoA", "Bray_Feed"), color = "Treatment..", shape = "Treatment..") 
plot_ordination(ZDSFdata_Feed, ordinate(ZDSFdata_Feed, "PCoA", "Bray_Feed"), color = "Treatment..", shape = "Treatment..")
plot_ordination(MGSFdata_Feed, ordinate(MGSFdata_Feed, "PCoA", "Bray_Feed"), color = "Treatment..", shape = "Treatment..")

#Check your library size differences by dividing the sequences sum of the sample with the most sequences by the one with the least sequences typically want it to be <10x difference


max(sample_sums(Fdata_Feed))/min(sample_sums(Fdata_Feed))#8.68932
max(sample_sums(SRFdata_Feed))/min(sample_sums(SRFdata_Feed)) #5.767051
max(sample_sums(ZDSFdata_Feed))/min(sample_sums(ZDSFdata_Feed)) #3.762684
max(sample_sums(MGSFdata_Feed))/min(sample_sums(MGSFdata_Feed)) #2.093023

#Check to see if they are statistically normal, if >0.05 then it is likely normal
shapiro.test(sample_sums(Fdata_Feed)) #0.9403
shapiro.test(sample_sums(SRFdata_Feed)) #0.9983
shapiro.test(sample_sums(ZDSFdata_Feed)) #1
shapiro.test(sample_sums(MGSFdata_Feed)) #0.2568

#Visual representation of the sample sums
hist(sample_sums(Fdata_Feed))
hist(sample_sums(SRFdata_Feed))
hist(sample_sums(ZDSFdata_Feed))
hist(sample_sums(MGSFdata_Feed))

#ZDSFdata_Feed performed the best, having 2nd lowest ratio, and more obvious normally distributed histogram

###BETA DIVERSITY ANALYSIS SQUARE ROOT VERSION###

##Export ESV table for PRIMER7##

# Extract abundance matrix from the phyloseq object#
QIIME2_WS_Bio = as(otu_table(Fdata_Feed), "matrix")

# Coerce to data.frame#
QIIME2_WS_Bio = as.data.frame(QIIME2_WS_Bio)

#Make a txt file#
write.table(QIIME2_WS_Bio, file = 'QIIME2_WS_Bio.txt', sep = "\t")

##PERMANOVA Analysis##

#make a distance matrix#
Bray_Feed<-phyloseq::distance(ZDSFdata_Feed, "Bray")

#convert metadata of the phyloseq to a dataframe#
sampledf_Feed<-data.frame(sample_data(ZDSFdata_Feed))

#run a PERMANVOA using vegan's adonis function#
adonis2(Bray_Feed~Treatment.., data=sampledf_Feed, permutations=9999)

# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 9999
# 
# adonis2(formula = Bray_Feed ~ Treatment.., data = sampledf_Feed, permutations = 9999)
# Df SumOfSqs      R2      F Pr(>F)  
# Treatment..  3   1.0063 0.30353 1.1622 0.0392 *
#   Residual     8   2.3091 0.69647                
# Total       11   3.3154 1.00000                
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

##PCoA Analysis##

#Make a dissimilarity/similarity matrix, by default Bray_Feed Curtis is used#
ZDSFdata_FeedBrayPCoAO <- ordinate(ZDSFdata_Feed,"PCoA")

PCoA_Bray_Feed_Treatment <- plot_ordination(ZDSFdata_Feed, ZDSFdata_FeedBrayPCoAO, color="Treatment..", shape="Treatment..")+
  geom_point(size=7, alpha=0.75)+
  ggtitle("QIIME2 Workshop PCOA with Bray-Curtis")+
  scale_colour_manual(name="Treatment", values=c("blue", "rosybrown", "purple", "orange", "forestgreen")) +
  scale_shape_manual(name="Treatment", values =c(15, 16, 17, 18))
PCoA_Bray_Feed_Treatment

##Save Graph for Publication##
tiff('PCoA by Treatment.tiff', units="in", width=10, height=6, res=300)
PCoA_Bray_Feed_Treatment
dev.off()

###TAXONOMIC BAR PLOTS###

# get abundance in relative percentage#
phy_Feed<- transform_sample_counts(Fdata_Feed, function(x) 100*x/sum(x))

##Create table ready for making stacked bar graph for Phyla >1% by Site-Survey combo##

# agglomerate taxa such that the ASVs with the same upper taxonomic ranks are merged at whatever rank you choose#
glom <- tax_glom(phy_Feed, taxrank = 'Phylum')

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
Phylum_Fdata_Feed_Treatment <- rbind(dat, Abundance)

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
Genus_Fdata_Feed_Treatment <- rbind(dat, Abundance)


##Make a series of graphs using ggplot based on above tables##

#Make Treatment Phylum graph#

spatial_plot_Phylum_Fdata_Feed_Treatment <- ggplot(data=Phylum_Fdata_Feed_Treatment, aes(x=Treatment.., y=Abundance, fill=Phylum)) + 
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  #theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("Phylum", values=c("orange", "lightblue", "purple", "firebrick3", "forestgreen", "yellow"))+
  xlab("Treatment") +
  ylab("Percentage")
#ggtitle("Phyla by Treatment")+
#theme(plot.title = element_text(hjust = 0.5)) +
#scale_x_discrete(limits=c("FP-Pre-Hurricane", "FP-Post-Hurricane", "HT-Pre-Hurricane", "HT-Post-Hurricane"), labels=c("Fort Pierce before Irma", "Fort Pierce after Irma", "Harbortown Marina before Irma", "Harbortown Marina after Irma"))
spatial_plot_Phylum_Fdata_Feed_Treatment

##Save Graph for Publication##
tiff('Phylum by Treatment.tiff', units="in", width=10, height=6, res=300)
spatial_plot_Phylum_Fdata_Feed_Treatment
dev.off()

##Make Genus Graph##

spatial_plot_Genus_Fdata_Feed_Treatment <- ggplot(data=Genus_Fdata_Feed_Treatment, aes(x=Treatment.., y=Abundance, fill=Genus)) +
  geom_bar(aes(), stat="identity", position="stack", color="black") + 
  # theme(axis.text.x  = element_text(angle=45, hjust=1)) +
  scale_fill_manual("Family;Genus", values=c("deeppink", "midnightblue", "blue", "red", "yellow", "orange", "green", "maroon", "dodgerblue", "purple", "lightblue", "orangered", "firebrick3", "darkseagreen3", "navy", "black", "seagreen", "turquoise", "white", "rosybrown"))+
  #ggtitle("Genera by Treatment") +
  xlab("Treatment") +
  ylab("Percentage") 
#theme(plot.title = element_text(hjust = 0.5)) +
#scale_x_discrete(limits=c("FP-Pre-Hurricane", "FP-Post-Hurricane", "HT-Pre-Hurricane", "HT-Post-Hurricane"), labels=c("Fort Pierce before Irma", "Fort Pierce after Irma", "Harbortown Marina before Irma", "Harbortown Marina after Irma"))
spatial_plot_Genus_Fdata_Feed_Treatment

##Save Graph for Publication##
tiff('Family;Genera by Treatment.tiff', units="in", width=10, height=6, res=300)
spatial_plot_Genus_Fdata_Feed_Treatment
dev.off()

#You will have to probably adjust the plot area 

#Determine Genus Fold Changes

#Switch the structure of dataframe from long to wide
wide_Genus_Fdata_Feed_Treatment <- spread(Genus_Fdata_Feed_Treatment, Treatment.., Abundance) 
wide_Genus_Fdata_Feed_Treatment

#Add a new column that is the TM1 abundance of a particular genera divided by the TM4 abundance of that same genera
wide_Genus_Fdata_Feed_Treatment$TM1TM4 <- wide_Genus_Fdata_Feed_Treatment$TM1 / wide_Genus_Fdata_Feed_Treatment$TM4 

#Same thing except that its TM4 divided by TM1
wide_Genus_Fdata_Feed_Treatment$TM4TM1 <- wide_Genus_Fdata_Feed_Treatment$TM4 / wide_Genus_Fdata_Feed_Treatment$TM1

#TM1 divided by TM2
wide_Genus_Fdata_Feed_Treatment$TM1TM2 <- wide_Genus_Fdata_Feed_Treatment$TM1 / wide_Genus_Fdata_Feed_Treatment$TM2 

#TM2 divided by TM1
wide_Genus_Fdata_Feed_Treatment$TM2TM1 <- wide_Genus_Fdata_Feed_Treatment$TM2 / wide_Genus_Fdata_Feed_Treatment$TM1

#TM1 divided by TM3
wide_Genus_Fdata_Feed_Treatment$TM1TM3 <- wide_Genus_Fdata_Feed_Treatment$TM1 / wide_Genus_Fdata_Feed_Treatment$TM3 

#TM3 divided by TM1
wide_Genus_Fdata_Feed_Treatment$TM3TM1 <- wide_Genus_Fdata_Feed_Treatment$TM3 / wide_Genus_Fdata_Feed_Treatment$TM1

###ALPHA DIVERSITY###


##Calculate alpha diversity##

Falphadiv_Feed = estimate_richness(Fdata_Feed, split = TRUE)
rownames(Falphadiv_Feed) <- c("17-1",  "17-2",  "17-3", "18-1",  "18-2",  "18-3", "19-1",  "19-2",  "19-3", "20-1",  "20-2",  "20-3")
Falphadiv_Feed
write.csv(Falphadiv_Feed, file='Falphadiv_Feed.csv')

#cbind in the metadata from the filtered phyloseq#
Falphadiv_Feed_metadata = cbind(as(Falphadiv_Feed, "data.frame"), as(sample_data(Fdata_Feed)[rownames(Falphadiv_Feed), ], "data.frame"))

##Test for normality##
shapiro.test(Falphadiv_Feed_metadata$Shannon) #0.8998
shapiro.test(Falphadiv_Feed_metadata$Observed) #0.9367
shapiro.test(Falphadiv_Feed_metadata$Fisher) #0.9473
shapiro.test(Falphadiv_Feed_metadata$Simpson) #0.3316

#Shannon and Simpson are not normally distributed thus must use Pearson for correlation tests

##Correlation Tests##
cor.test(Falphadiv_Feed_metadata$Shannon, Falphadiv_Feed_metadata$Observed, method = "pearson", exact = FALSE) #6.773e-05
cor.test(Falphadiv_Feed_metadata$Shannon, Falphadiv_Feed_metadata$Fisher, method = "pearson", exact = FALSE) #3.062e-05
cor.test(Falphadiv_Feed_metadata$Shannon, Falphadiv_Feed_metadata$Simpson, method = "pearson", exact = FALSE) #0.0001626

##Multiple test adjustment##
#Make a list of all the p vales from the correlation tests, make it a dataframe, add column of adjusted p values#
alpha_div_p.value_Feed=c(6.773e-05,3.062e-05,0.0001626)
alpha_div_p.value_Feed=data.frame(alpha_div_p.value_Feed)
alpha_div_p.value_Feed$padj <- p.adjust(alpha_div_p.value_Feed$alpha_div_p.value_Feed, method = "BH")
alpha_div_p.value_Feed

# alpha_div_p.value_Feed        padj
# 1              6.773e-05 0.000101595
# 2              3.062e-05 0.000091860
# 3              1.626e-04 0.000162600

#Shannon is significantly correlated to other alpha diversity metrics

#Besides normality, another good thing to test for statistical analysis is heteroscedacity
bartlett.test(Shannon ~ Treatment.., Falphadiv_Feed_metadata) #0.8565

#Shannon is normal and homoscedastic, will use ANOVA and Tukey to test

##Basic visualization of Shannon statistics##
ddply(Falphadiv_Feed_metadata, .(Treatment..), summarise, mean=mean(Shannon), sd=sd(Shannon), median=median(Shannon), IQR=IQR(Shannon))

# Treatment..     mean        sd   median       IQR
# 1         TM1 3.416990 0.7986628 3.677094 0.3904523
# 2         TM2 2.984462 1.0495770 3.331221 0.8505816
# 3         TM3 2.677683 1.4891391 3.430333 1.0766987
# 4         TM4 3.512199 0.2211072 3.539258 0.1943155

#Make a simple boxplot to get a visual summary of these statistics#
ggplot(Falphadiv_Feed_metadata, aes(x=Treatment.., y=Shannon, color=Treatment..)) + 
  geom_boxplot() +
  theme(axis.text.x  = element_text(angle=45, hjust=1))

###TESTS FOR TOTAL AND PAIRWISE SIGNIFICANCE###
#First store results in a variable, and then show a summary of the results
res.aov_Feed <- aov(Shannon ~ Treatment.., data = Falphadiv_metadata)
anova(res.aov_Feed) #0.3257

#Alpha diversity is not significant, so will not test pair-wise Treatments

###COMBINING STATISTICAL AND ALPHA DIVERSITY TESTING###

##Making a more informative boxplot##

#Most of this you've seen before with bar plot#
Shannon_Treatment_Boxplots_Feed <- ggplot(Falphadiv_Feed_metadata, aes(x=Treatment.., y=Shannon)) + 
  geom_boxplot() +
  #theme(axis.text.x  = element_text(angle=45, hjust=1)) + 
  xlab("Treatment") +
  #scale_x_discrete(limits=c("FP-Pre-Hurricane", "FP-Post-Hurricane", "HT-Pre-Hurricane", "HT-Post-Hurricane")) +
  annotate("text", x = c(1:4) , y = c(2.8, 3.6, 3.25, 3.2), label = c("a", "a", "a", "a"), size=5)+
  theme(axis.text=element_text(size=13), axis.title=element_text(size=15), legend.text=element_text(size=11), legend.key.size = unit(0.7,"cm"), legend.title=element_text(size=13))
Shannon_Treatment_Boxplots_Feed

tiff('Shannon by Treatment.tiff', units="in", width=10, height=6, res=300)
Shannon_Treatment_Boxplots_Feed
dev.off()