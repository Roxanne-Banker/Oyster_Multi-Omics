#rm(list = ls())

setwd("/Users/roxannebanker/bioinformatics/projects/oyster-oa/microbiome/")


# February 2020
# Document to analyze oyster larvae/juvenile data from 
# Oyster multi-omics project run in Aug-Oct 2019
# Following this tutorial: https://benjjneb.github.io/dada2/tutorial.html
# and the Harbor_R_Notebook from Cassie Ettinger

# loading in required packages

library(phyloseq)
#packageVersion("phyloseq")
library(ggplot2)
packageVersion("ggplot2")
library(rmarkdown)
library(dada2)
#packageVersion("dada2")
library(decontam)
#packageVersion("decontam")
library(phangorn)
library(DECIPHER)
library(vegan)
#packageVersion("vegan")
library(dplyr)
library(RColorBrewer)
library(reshape)
library(coin)
library(FSA)
#packageVersion("FSA")
library(Biostrings)
library(DESeq2)
library(decontam)
library(gridExtra)
library(seqinr)
#salibrary(QsRutils)


# Set seed to make results reproducible
set.seed(5311)



# Using Dada2 to create an amplicon sequence variant (ASV) table
######################################################

# Setting the path for the fastq files
path <- "/Users/roxannebanker/bioinformatics/projects/oyster-oa/microbiome/raw_data/" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

# Sort and get sample names
fnFs <- sort(list.files(path, pattern="_R1_001.fastq"))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq"))
sample.names <- sapply(strsplit(fnFs, "_S"), `[`, 1)

# specify full paths to the data
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)

# Quality Inspection
# Looking at our data to try to see at what point the read quality drops below Q20

plotQualityProfile(fnFs[1:4])
# For the forward reads, the controls look a bit bad, but the other samples
# look good. How about the reverse reads?

plotQualityProfile(fnRs[1:4])
# These are maybe worse, but tutorials indicate that is okay for reverse reads

# Note: might have to be careful in subsequent trimming steps
# The tutorial this script is following indicates that the tutorial 
# was written for 2x250 V4 sequence data, which means that the forward and
# reverse reads almost completely overlap, and can use quality scores to 
# guide filtering. I am not sure how well our primer set, V4-V5, overlaps. 
# I am going to try what the tutorial uses, then adjust
# the truncLen from there.

# tutorial "your truncLen must be large enough to maintain 20 + 
# biological.length.variation nucleotides of overlap between them."
######################################################


# Quality Filtering
######################################################

# Assign the filenames for the filtered fastq.gz files.
# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(300,240),
                     maxN=0, maxEE=c(3,6), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE,
                     trimLeft = c(19, 20)) # On Windows set multithread=FALSE
head(out)

# this is the first attempt with the default truncLen values in the tutorial
# truncLen=c(240,160)
# reads.in reads.out
# b-A_S179_L001_R1_001.fastq       358       109
# b-B_S191_L001_R1_001.fastq       329        94
# C1-A_S332_L001_R1_001.fastq    61807     56665
# C1-B_S344_L001_R1_001.fastq    82331     75507
# C1-C_S108_L001_R1_001.fastq   103194     95612
# C2-A_S120_L001_R1_001.fastq   111945    103285

# truncLen=c(200,250)
# reads.in reads.out
# b-A_S179_L001_R1_001.fastq       358        72
# b-B_S191_L001_R1_001.fastq       329        66
# C1-A_S332_L001_R1_001.fastq    61807     50486
# C1-B_S344_L001_R1_001.fastq    82331     67597
# C1-C_S108_L001_R1_001.fastq   103194     88307
# C2-A_S120_L001_R1_001.fastq   111945     95590

# truncLen=c(300,300)
# reads.in reads.out
# b-A_S179_L001_R1_001.fastq       358        41
# b-B_S191_L001_R1_001.fastq       329        42
# C1-A_S332_L001_R1_001.fastq    61807     37846
# C1-B_S344_L001_R1_001.fastq    82331     52180
# C1-C_S108_L001_R1_001.fastq   103194     70509
# C2-A_S120_L001_R1_001.fastq   111945     74158

# Settled on these parameters for now:
# truncLen=c(300,240)
# maxEE=c(3,6)
# reads.in reads.out
# b-A_S179_L001_R1_001.fastq       358       129
# b-B_S191_L001_R1_001.fastq       329       111
# C1-A_S332_L001_R1_001.fastq    61807     57342
# C1-B_S344_L001_R1_001.fastq    82331     76414
# C1-C_S108_L001_R1_001.fastq   103194     96472
# C2-A_S120_L001_R1_001.fastq   111945    104226



# Learn the Error Rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# It is always worthwhile, as a sanity check if nothing else, to visualize the estimated error rates:
plotErrors(errF, nominalQ=TRUE)


# Sample Inference

# We are now ready to apply the core sample inference algorithm 
# to the filtered and trimmed sequence data.
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

# Inspecting the returned dada-class object:
dadaFs[[5]]

# Merge paired reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)


# Inspect the merger data.frame from the first sample
head(mergers[[4]])


# I am not sure if enough merged, but I can revisit after I ask David or Cassie
######################################################


# Construct sequence table
######################################################

# We can now construct an amplicon sequence variant table (ASV) table, 
# a higher-resolution version of the OTU table produced 
# by traditional methods.
seqtab <- makeSequenceTable(mergers)

dim(seqtab)
# [1]   38 7916


# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# 287  312  326  341  351  352  356  359  362  364  365  366  367  368  369  370 
# 1    1    1    3    1    1    1    1    4    2    6   12    7  107   63  623 
# 371  372  373  374  375  376  377  378  379  380  381  382  383  389  390  392 
# 145 2445 1639 2592  142   75    8    2    1    4    2    1    4    1    1    1 
# 394  396  398  399  400  401  402  403  406  408  440  449  459 
# 2    1    1    3    1    2    1    1    2    2    1    1    1 

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

dim(seqtab.nochim)
# [1]   38 2043

sum(seqtab.nochim)/sum(seqtab)
# [1] 0.8725215


# Track reads through the pipeline. As a final check of our progress, we’ll look at the number of 
# reads that made it through each step in the pipeline:
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
print(track)

# input filtered denoisedF denoisedR merged nonchim
# b-A      358      131        66        69     29      29
# b-B      329      112        64        59     29      29
# C1-A   61807    57394     56873     56814  40971   40188
# C1-B   82331    76506     75615     75635  57865   56057
# C1-C  103194    96545     95275     95332  58883   55748
# C2-A  111945   104313    103567    103599  52138   51056
# C2-B   48171    44209     43860     43851  25644   25501
# C2-C   52745    49208     48743     48827  38736   38329
# C3-A   72834    68929     67957     68029  57528   53841
# C3-B  181999   171222    169612    169923 155405  108690
# C3-C   89734    83726     82849     82819  54065   50487
# L-A    18247    16587     16370     16404   2882    2831
# L-B    29787    27191     26849     26911   5739    5669
# L-C    26955    24461     24299     24321   1899    1861
# OA1-A  51681    48228     47693     47660  41232   40496
# OA1-B  94592    89076     88636     88690  80121   79345
# OA1-C  62762    59207     58630     58647  49466   48652
# OA2-A  48235    45093     44815     44819  28676   28318
# OA2-B  98698    92921     92417     92447  63522   62162
# OA2-C  73463    69087     68635     68616  57009   55294
# OA3-A  67733    62698     62328     62325  23570   23389
# OA3-B  71386    66405     66068     66074  25761   25613
# OA3-C  87453    81346     80860     80898  35360   35041
# P-A   173770   161972    161394    161639 153915  112375
# P-B   194384   180920    180115    180450 171954  121895
# P-C   225269   210796    210218    210455 201179  142924
# SM1-A    153       59         6        17      0       0
# SM1-B  40809    38304     37988     38011  28854   28529
# SM1-C  85050    79772     79168     79197  62430   61080
# SM2-A    206       59        25        25     13      13
# SM2-B  52845    49388     48995     49079  34543   34096
# SM2-C  53553    50634     50267     50291  44195   43398
# SM3-A 102784    96386     94966     95067  72161   67314
# SM3-B 165704   156444    154557    154547 136051  122201
# SM3-C  92555    87310     85756     85805  69615   62767
# t-A      144       53        13        22     12      12
# t-B      265      102        34        37     20      20
# t-C      374      113         9        33      7       7


# Some improvement on sequence retention after adjusting Trunclen and maxEE parameters.
# But is this enough???
######################################################


# Assign taxonomy
######################################################

taxa <- assignTaxonomy(seqtab.nochim, "~/bioinformatics/projects/oyster-oa/microbiome/raw_data/tax/silva_nr_v132_train_set.fa.gz", multithread=TRUE)

# Assign species identity
taxa <- addSpecies(taxa, "~/bioinformatics/projects/oyster-oa/microbiome/raw_data/tax/silva_species_assignment_v132.fa.gz")

# Let’s inspect the taxonomic assignments:
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

######################################################


# Evaluate accuracy
######################################################

# The tutorial has mock communities to assess accuracy of tax assignment
# for my data set these are P-A, P-B, and P-C

# I am going to revisit how to use the samples later. For now, some fun stuff.
######################################################



# Build Phylogeny
# All from Ettinger Harbor Notebook
# March 26 2020: This is computationally expensive and I am thus going to skip right now.
# So far, it seems like the need for the phylogeny comes from calculating Unifrac distances
# I plan to use Bray-Curtis, so perhaps not needed right now
#################################################################


#get DNA sequences from chimera-removed ASV table
seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs # This propagates to the tip labels of the tree

#Align DNA sequences
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)

#Turn alignment into a matrix
phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
#calculate maximum liklihood values for alignment
dm <- dist.ml(phangAlign)


#build neighbor joining tree
#a GTR+G+I (Generalized time-reversible with Gamma rate variation)
#maximum likelihood tree using the neighbor-joining tree
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE,
                    optGamma=TRUE, rearrangement = "stochastic", control = pml.control(trace = 0))
# need to rerun, took too long ^

#saveRDS(fitGTR, "Harbor_seqtab2_tree.rds")
#Now to root the tree
fitGTR <- readRDS("Harbor_seqtab2_tree.rds")
#using most abundant archaea as outgroup
outgroupseq <- "CCGGTTCGTGCCCCTAGGCTTCGTCCCTGACCGTCGGATCCGTTCCAGTGTGGCGCCTTCGCAACTGGTGGTCCTCTAAGGATGACAGGATTTTACCCCTACCCCTAGAGTACCCCACACCTCTCCC"

phyloseq.tree <- root(fitGTR$tree, outgroup = outgroupseq, resolve.root = TRUE)

saveRDS(phyloseq.tree, "Harbor_seqtab2_tree_archaea_rooted.rds")

#################################################################







# Phyloseq


# Identify Contaminants
#################################################################
sample_names(ps) # bead blanks: 1,2; tube blanks: 36,37,38
vector_for_decontam <- c(rep(TRUE, 2), rep(FALSE, 33), rep(TRUE, 3))

contam_df <- isContaminant(ps, neg=vector_for_decontam)

table(contam_df$contaminant) # identified 1 as contaminant

#which ASV is contaminant
head(which(contam_df$contaminant))
#seq 1
contam_df[1,]
# TACGTAGGTGGCAAGCGTTGTCCGGAATTATTGGGCGTAAAGGGCTCGCAGGCGGTTCCTTAAGTCTGATGTGAAAGCCCCCGGCTCAACCGGGGAGGGTCATTGGAAACTGGGGAACTTGAGTGCAGAAGAGGAGAGTGGAATTCCACGTGTAGCGGTGAAATGCGTAGAGATGTGGAGGAACACCAGTGGCGAAGGCGACTCTCTGGTCTGTAACTGACGCTGAGGAGCGAAAGCGTGGGGAGCGAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGAGTGCTAAGTGTTAGGGGGTTTCCGCCCCTTAGTGCTGCAGCTAACGCATTAAGCACTCCGCCTGGGGAGTACGGTCGCAAGACTG 

# make a list of the contaminant
contaminant <- ("TACGTAGGTGGCAAGCGTTGTCCGGAATTATTGGGCGTAAAGGGCTCGCAGGCGGTTCCTTAAGTCTGATGTGAAAGCCCCCGGCTCAACCGGGGAGGGTCATTGGAAACTGGGGAACTTGAGTGCAGAAGAGGAGAGTGGAATTCCACGTGTAGCGGTGAAATGCGTAGAGATGTGGAGGAACACCAGTGGCGAAGGCGACTCTCTGGTCTGTAACTGACGCTGAGGAGCGAAAGCGTGGGGAGCGAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGAGTGCTAAGTGTTAGGGGGTTTCCGCCCCTTAGTGCTGCAGCTAACGCATTAAGCACTCCGCCTGGGGAGTACGGTCGCAAGACTG")
#################################################################



# Creating phyloseq object
#################################################################

# reading in mapping file
mapping <- read.csv("phyloseq-metadata.csv")

row.names(mapping) <- mapping$sample_names

# We now construct a phyloseq object directly from the dada2 outputs.
otu_table = otu_table(seqtab.nochim, taxa_are_rows=FALSE)
mapping_file = sample_data(mapping)
taxa_table = tax_table(taxa)

# creating phyloseq object
ps <- phyloseq(otu_table,mapping_file, taxa_table)
#################################################################



# Remove contaminant
#################################################################
#get all taxa names
allTaxa = taxa_names(ps)
#returns list of all taxa, except contaminants
allTaxa <- allTaxa[!(allTaxa %in% contaminant)]
#keeps taxa in allTaxa and gets rid of anthing else
ps_filt = prune_taxa(allTaxa, ps)
#################################################################



# Remove Chloroplasts, Mitochondria, and Animal sequences
#################################################################
# removing chloroplasts, mitochondria FOR SILVA 132 - silva 128 has chloro as a Class
ps_noChloro <- subset_taxa(ps_filt, Order !="o__Chloroplast")
ps_noChloro_noMito <- subset_taxa(ps_noChloro, Family != "f__Mitochondria")
pssub_noChloro_noMito_noMollusc <- subset_taxa(ps_noChloro_noMito, Phylum != "p__Mollusca")
ps_noChloro_noMito_noMollusc_noAni <- subset_taxa(ps_noChloro_noMito_noMollusc, Kingdom != "p__Animalia")
# checking to see how many sequences lost during these filtering steps
filter_seq_loss <- (1-(sample_sums(ps_noChloro_noMito_noMollusc_noAni)/sample_sums(ps)))*100
print(filter_seq_loss)
# need to remove samples that have a zero count, I think it is disrupting deseq2 down the line
data <- prune_samples( sample_sums(ps_noChloro_noMito_noMollusc_noAni) > 0, ps_noChloro_noMito_noMollusc_noAni )
#################################################################



# Summary Function
#################################################################
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
#################################################################



# Alpha diversity plots
#################################################################
# creating a list to subset phyloseq object to remove negative and positive controls
names <- sample_names(data)
names <- names[-c(1,2,24,25,26,29,35,36)]

ps.alpha <- prune_samples(names,data)

alpha.measure <- estimate_richness(ps.alpha, measures=c("Observed", "Shannon"))
alpha.measure$sample_name <- sample_names

mapping$sample_names <- as.character(mapping$sample_names)
alpha.measure.full <- full_join(alpha.measure, mapping, by = c("sample_name" = "sample_names"))
alpha.measure.full <- na.omit(alpha.measure.full)

observed <- summarySE(alpha.measure.full, measurevar="Observed", groupvars=c("bucket.id"))
shannon <- summarySE(alpha.measure.full, measurevar="Shannon", groupvars=c("bucket.id"))

alpha_treatment <- c(rep("control", 3), rep("larvae", 1), rep("low-ph", 3), rep("molybdate", 3))
shannon$treatment <- alpha_treatment
observed$treatment <- alpha_treatment



alpha.breaks <- as.character(observed$bucket.id)

alpha.ids <- c("C38-1",  "C38-2",  "C51",  "L",  
                "LP38-1", "LP28-2", "LP51", 
                "SM38-1", "SM38-2", "SM51")

alpha.color <-c("#8c510a", "#bf812d", "#dfc27d", "#636363",
                 "#800026", "#e31a1c", "#fd8d3c",
                 "#253494", "#1d91c0", "#7fcdbb")

# factor ordering so the x axis is in the order i want
observed$bucket.id <- factor(observed$bucket.id,levels = c("L", "C1", "C2", "C3", "OA1", "OA2", "OA3", "SM1", "SM2", "SM3"))


# Simpson
observed.plot <- ggplot(observed, aes(x=bucket.id, y=Observed)) +
  geom_point() +
  geom_errorbar(aes(ymin=Observed-se, ymax=Observed+se),
                width=.1, size=0.4,                    # Width of the error bars
                position=position_dodge(.9)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  labs(x = "", y="Observed ASVs") +
  scale_x_discrete(breaks=NULL,
                   labels=NULL) +
  theme( legend.position = "none" )

observed.plot

# factor ordering so the x axis is in the order i want
shannon$bucket.id <- factor(shannon$bucket.id,levels = c("L", "C1", "C2", "C3", "OA1", "OA2", "OA3", "SM1", "SM2", "SM3"))

# Shannon
shannon.plot <- ggplot(shannon, aes(x=bucket.id, y=Shannon)) +
  geom_point() +
  geom_errorbar(aes(ymin=Shannon-se, ymax=Shannon+se),
                width=.1, size=0.4,                    # Width of the error bars
                position=position_dodge(.9)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  labs(x = "Bucket ID", y="Shannon Index") +
  scale_x_discrete(breaks=alpha.breaks,
                   labels=alpha.ids) +
  theme( legend.position = "none",
         axis.text.x = element_text(face=2))

shannon.plot

grid.arrange(observed.plot, shannon.plot, nrow = 2)




# Alpha Diversity Statistics

# bucket ID
kruskal_test(Shannon ~ bucket.id, distribution = approximate(nresample = 9999), data=shannon)
# chi-squared = 9, p-value = 1

kruskal_test(Observed ~ bucket.id, distribution = approximate(nresample = 9999), data=observed)
# chi-squared = 9, p-value = 1


# need to make the treatment colums factors for this analysis
shannon$treatment <- as.factor(shannon$treatment)
observed$treatment <- as.factor(observed$treatment)

# treatment
kruskal_test(Shannon ~ treatment, distribution = approximate(nresample = 9999), data=shannon)
# chi-squared = 1.5818, p-value = 0.7597

kruskal_test(Observed ~ treatment, distribution = approximate(nresample = 9999), data=observed)
# chi-squared = 3.4, p-value = 0.3723

#################################################################





# Normalizing for sample depth via variance stabalizing transformation
#################################################################
# Will be using Deseq2 to variance stabilizing transformation, instead of the classic rarefaction or
# turning counts into proportions. This is recommended by the McMurdie and Holmes 2014 Plos Computational Biology paper, Waste not Want not.
# using this tutorial: https://astrobiomike.github.io/amplicon/dada2_workflow_ex#analysis-in-r

# first we need to make a DESeq2 object
# the design I set, bucket.id, will create groups based on both treatment and time (point sampled)
data_deseq <- phyloseq_to_deseq2(data, ~ bucket.id)

gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
geoMeans = apply(counts(data_deseq), 1, gm_mean)
data_deseq = estimateSizeFactors(data_deseq, geoMeans = geoMeans, type="poscount")
#data_deseq = DESeq(data_deseq, test="Wald", fitType="parametric")
# great this worked!

data_deseq_vst <- varianceStabilizingTransformation(data_deseq, blind = FALSE, fitType = "parametric")

# and here is pulling out our transformed table
vst_transform <- assay(data_deseq_vst)

# okay. the VST log-like transformation produces negative values
# for counts that are <1
# therefore I am going to set these to zero so that distance metrics can be applied
# https://www.bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-FAQ.html#negative-numbers-in-my-transformed-data-table

vst_transform[vst_transform < 0.0] <- 0.0

# and calculating our Euclidean distance matrix
euc_dist <- assay(data_deseq_vst)
euc_dist <- dist(t(vst_transform))

#################################################################




# Hierarchical clustering
#################################################################
# Astro Tutortial: Now that we have our Euclidean distance matrix, let’s make and plot a hierarchical clustering of our samples.

euc_clust <- hclust(euc_dist, method="ward.D2")

euc_dend <- as.dendrogram(euc_clust, hang=0.1)
dend_cols <- as.character(sample_info_tab$color[order.dendrogram(euc_dend)])
labels_colors(euc_dend) <- dend_cols

plot(euc_dend, ylab="VST Euc. dist.")

#################################################################




# New Phyloseq object with adjusted ASV counts via DESEQ2
#################################################################
# making our phyloseq object with transformed table
vst_transform_phy <- otu_table(vst_transform, taxa_are_rows=T)

# creating phyloseq object
ps.vst <- phyloseq(vst_transform_phy, mapping_file, taxa_table)

# Calculate RA per sample
ps.vst.RA = transform_sample_counts(ps.vst, function(x) 100 * x/sum(x))

#################################################################




# Ordinatation
#################################################################

# ordination using vst adjust sample counts: ps.vst

# filtering out the postitive and negative controls
ps.filt <- subset_samples( ps.vst.RA, sample_data(ps.vst.RA)$age >= 1)
# removing SM2A because it has low counts, throws off analyses
ps.all <- subset_samples( ps.filt, sample_names(ps.filt) != "SM2-A")



# Ordination WITH larval samples
# Bray-Curtis Distances PCoA with entire dataset
ord.pcoa.bray <- ordinate(ps.all, method="PCoA", distance="bray")

# plot
plot_ordination(ps.all, ord.pcoa.bray , shape="treatment", color="treatment", title="PCoA Bray ")

# interesting! now to spruce it up a bit

bucket.ids <- c("C1",  "C2",  "C3",  "L",  
                "OA1", "OA2", "OA3", 
                "SM1", "SM2", "SM3")

bucket.labs <- c("C38-1",  "C38-2",  "C51",  "L",  
               "LP38-1", "LP28-2", "LP51", 
               "SM38-1", "SM38-2", "SM51")

bucket.color <-c("#8c510a", "#bf812d", "#dfc27d", "#636363",
                 "#800026", "#e31a1c", "#fd8d3c",
                 "#253494", "#1d91c0", "#7fcdbb")


bucket.shapes <-c(16, 16, 16, 3,
                  15, 15, 15,
                  17, 17, 17)

bucket.fill <-c("#8c510a", "#bf812d", "#ffffff", "#636363",
                "#003c30", "#01665e", "#ffffff",
                "#762a83", "#9970ab", "#ffffff")

# might be a better way to do this, but given I want colors to be
# CB friendly, doing it this way

plot_ordination(physeq = ps.filt, ordination = ord.pcoa.bray, color="bucket.id", shape="bucket.id") +
  #stat_ellipse(aes(group=treatment), linetype = 2, type="t") + #t distribution at level 0.95
  #scale_colour_brewer(type = "qual", palette = "Set1") +
  geom_point(size = 4) +
  theme(legend.title=element_blank()) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(name="Bucket ID", labels = bucket.labs,
                     values= bucket.color) +
  scale_shape_manual(name="Bucket ID", values=bucket.shapes, guide=NULL) +
  scale_fill_manual(name="Bucket.ID", values=bucket.fill) +
  xlab("PCoA 1 [42.1% of variation]") + 
  ylab("PCoA 2 [16.8% of variation]")



# Ordination WITHOUT larval samples
# filtering out the larval samples
ps.juv <- subset_samples( ps.all, sample_data(ps.filt)$age >= 5)


# Bray-Curtis Distances PCoA with entire dataset
juv.ord.pcoa.bray <- ordinate(ps.juv, method="PCoA", distance="bray")

# plot
plot_ordination(ps.juv, juv.ord.pcoa.bray , shape="treatment", color="treatment", title="PCoA Bray ")


bucket.ids2 <- c("C1",  "C2",  "C3",  
                "OA1", "OA2", "OA3", 
                "SM1", "SM2", "SM3")

bucket.color2 <-c("#8c510a", "#bf812d", "#dfc27d",
                  "#800026", "#e31a1c", "#fd8d3c",
                  "#253494", "#1d91c0", "#7fcdbb")

bucket.shapes2 <-c(16, 16, 16,
                  15, 15, 15,
                  17, 17, 17)


plot_ordination(physeq = ps.juv, ordination = juv.ord.pcoa.bray, color="bucket.id", shape="bucket.id") +
  #stat_ellipse(aes(group=treatment), linetype = 2, type="t") + #t distribution at level 0.95
  #scale_colour_brewer(type = "qual", palette = "Set1") +
  geom_point(size = 4) +
  theme(legend.title=element_blank()) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(name="Bucket ID", labels = bucket.labs,
                     values= bucket.color2) +
  scale_shape_manual(name="Bucket ID", values=bucket.shapes2, guide= NULL) +
  xlab("PCoA 1 [30.6% of variation]") + 
  ylab("PCoA 2 [19.0% of variation]")



# July 20, 2020
# Okay, now trying to produce a PCA instead

otu <- veganotu(ps.all)
otu.h <- decostand(otu, "hellinger")
otu.h.d <- vegdist(otu.h, "euclidean")

#################################################################




# Ordination Statistics (PERMANOVA, Levene's Test)

# all treatments, all time points within treatments
#################################################################


# the initial distance call
Dist.treat = phyloseq::distance(ps.all, method = "bray", type="samples")
# make a data frame from the sample_data
sample_data_filt <- data.frame(sample_data(ps.all))


# Treatments Comparison
# betadisper: dispersion test
treatment.beta <- betadisper(Dist.treat, sample_data_filt$treatment)
permutest(treatment.beta)
#            Df   Sum Sq   Mean Sq     F N.Perm Pr(>F)
# Groups     3 0.005837 0.0019458 0.692    999  0.566
# Residuals 24 0.067484 0.0028118  

# Adonis test
adonis(Dist.treat ~ treatment, data = sample_data_filt)
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# treatment  3    3.0171 1.00569  14.261 0.64063  0.001 ***
# Residuals 24    1.6924 0.07052         0.35937           
# Total     27    4.7095                 1.00000  


# Treatment-Time Comparison
# betadisper: dispersion test
treatment.time.beta <- betadisper(Dist.treat, sample_data_filt$group)
permutest(treatment.time.beta)
#            Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)  
# Groups     6 0.037320 0.0062199 3.2811    999  0.023 *
# Residuals 21 0.039809 0.0018957 

# Adonis test
adonis(Dist.treat ~ group, data = sample_data_filt)
#            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# group      6    3.8293 0.63822  13.647 0.79588  0.001 ***
# Residuals 21    0.9821 0.04677         0.20412           
# Total     27    4.8114                 1.00000   

# applying bonferrni correction to p values from both tests

# permanova
levene.p <- c(0.566, 0.023) 

levene.p.bonf <- p.adjust(levene.p, method="bonferroni")
levene.p.bonf

# permanova
permanova.p <- c(0.001, 0.001) 

permanova.p.bonf <- p.adjust(permanova.p, method="bonferroni")
permanova.p.bonf

#################################################################


# separating ps objects for treatment pairwise comparisons
#################################################################

ps.larv <- subset_samples(ps.all, treatment=="larvae")

ps.control <- subset_samples(ps.all, treatment=="control")

ps.lp <- subset_samples(ps.all, treatment=="low-ph")

ps.sm <- subset_samples(ps.all, treatment=="molybdate")

#################################################################


# pairwise between treatments
#################################################################

# Larva:Control
ps.larv.cntrl <- merge_phyloseq(ps.larv,ps.control)
# make a data frame from the sample_data
data.larv.cntrl <- data.frame(sample_data(ps.larv.cntrl))

Dist.larv.cntrl = phyloseq::distance(ps.larv.cntrl, method = "bray", type="samples")

# betadisper: dispersion test
larv.cntrl.beta <- betadisper(Dist.larv.cntrl, data.larv.cntrl$treatment)
permutest(larv.cntrl.beta)
#            Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
# Groups     1 0.002672 0.0026715 0.5692    999    0.5
# Residuals 10 0.046933 0.0046933      

# Adonis test
adonis(Dist.larv.cntrl ~ treatment, data = data.larv.cntrl)
#            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
# treatment  1    1.7923 1.79229  24.769 0.71239  0.003 **
# Residuals 10    0.7236 0.07236         0.28761          
# Total     11    2.5159                 1.00000 



# Larva:Low-pH 
ps.larv.oa <- merge_phyloseq(ps.larv,ps.lp)
# make a data frame from the sample_data
data.larv.oa <- data.frame(sample_data(ps.larv.oa))
Dist.larv.oa = phyloseq::distance(ps.larv.oa, method = "bray", type="samples")

# betadisper: dispersion test
larv.oa.beta <- betadisper(Dist.larv.oa, data.larv.oa$treatment)
permutest(larv.oa.beta)
#            Df    Sum Sq   Mean Sq     F N.Perm Pr(>F)
# Groups     1 0.0057708 0.0057708 2.749    999  0.132
# Residuals 10 0.0209927 0.0020993           

# Adonis test
adonis(Dist.larv.oa ~ treatment, data = data.larv.oa)
#            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
# treatment  1   1.77364 1.77364  22.662 0.69383  0.006 **
# Residuals 10   0.78265 0.07827         0.30617          
# Total     11   2.55629                 1.00000



# Larva:Molybdate
ps.larv.sm <- merge_phyloseq(ps.larv,ps.sm)
# make a data frame from the sample_data
data.larv.sm <- data.frame(sample_data(ps.larv.sm))
Dist.larv.sm = phyloseq::distance(ps.larv.sm, method = "bray", type="samples")

# betadisper: dispersion test
larv.sm.beta <- betadisper(Dist.larv.sm, data.larv.sm$treatment)
permutest(larv.sm.beta)
#            Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
# Groups     1 0.002687 0.0026870 0.9568    999  0.347
# Residuals  8 0.022467 0.0028083         

# Adonis test
adonis(Dist.larv.sm ~ treatment, data = data.larv.sm)
#            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# treatment  1   1.66530 1.66530  22.741 0.73976  0.014 *
# Residuals  8   0.58584 0.07323         0.26024         
# Total      9   2.25114                 1.00000 




# Control - OA
ps.cntrl.oa <- merge_phyloseq(ps.control,ps.lp)
# make a data frame from the sample_data
data.cntrl.oa <- data.frame(sample_data(ps.cntrl.oa))

Dist.cntrl.oa = phyloseq::distance(ps.cntrl.oa, method = "bray", type="samples")

# betadisper: dispersion test
cntrl.oa.beta <- betadisper(Dist.cntrl.oa, data.cntrl.oa$age)
permutest(cntrl.oa.beta)
#            Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)  
# Groups     1 0.007009 0.0070090 3.4164    999   0.08 .
# Residuals 16 0.032825 0.0020516            

# Adonis test
adonis(Dist.cntrl.oa ~ age, data = data.cntrl.oa)
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
# age        1    0.3562  0.3562  3.5909 0.18329   0.01 **
# Residuals 16    1.5871  0.0992         0.81671          
# Total     17    1.9433                 1.00000  



# Control - SM 
ps.cntrl.sm <- merge_phyloseq(ps.control,ps.sm)
# make a data frame from the sample_data
data.cntrl.sm <- data.frame(sample_data(ps.cntrl.sm))

Dist.cntrl.sm = phyloseq::distance(ps.cntrl.sm, method = "bray", type="samples")

# betadisper: dispersion test
cntrl.sm.beta <- betadisper(Dist.cntrl.sm, data.cntrl.sm$age)
permutest(cntrl.sm.beta)
#            Df    Sum Sq   Mean Sq      F N.Perm Pr(>F)
# Groups     1 0.0039854 0.0039854 2.1628    999  0.164
# Residuals 14 0.0257978 0.0018427             

# Adonis test
adonis(Dist.cntrl.sm ~ age, data = data.cntrl.sm)
#            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
# age        1   0.40936 0.40936  6.8104 0.32726  0.002 **
# Residuals 14   0.84150 0.06011         0.67274          
# Total     15   1.25086                 1.00000   



# OA - SM
ps.oa.sm <- merge_phyloseq(ps.lp,ps.sm)
# make a data frame from the sample_data
data.oa.sm <- data.frame(sample_data(ps.oa.sm))

Dist.oa.sm = phyloseq::distance(ps.oa.sm, method = "bray", type="samples")

# betadisper: dispersion test
oa.sm.beta <- betadisper(Dist.oa.sm, data.oa.sm$age)
permutest(oa.sm.beta)
#            Df    Sum Sq    Mean Sq     F N.Perm Pr(>F)  
# Groups     1 0.0026709 0.00267092 4.988    999  0.023 *
# Residuals 14 0.0074966 0.00053547        

# Adonis test
adonis(Dist.oa.sm ~ age, data = data.oa.sm)
#            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
# age        1   0.39272 0.39272  4.6116 0.24778  0.002 **
# Residuals 14   1.19223 0.08516         0.75222          
# Total     15   1.58495                 1.00000      


# applying bonferrni correction to p values from both tests

# permanova
levene.p <- c(0.5, 0.132, 0.347, 0.08, 0.164, 0.023) 

levene.p.bonf <- p.adjust(levene.p, method="bonferroni")
levene.p.bonf

# permanova
permanova.p <- c(0.003, 0.006, 0.014, 0.01, 0.002, 0.002) 

permanova.p.bonf <- p.adjust(permanova.p, method="bonferroni")
permanova.p.bonf

#################################################################


# separating ps objects for all time points
#################################################################

ps.l <- subset_samples(ps.all, group=="L")

ps.c38 <- subset_samples(ps.all, group=="C38")

ps.c51 <- subset_samples(ps.all, group=="C51")

ps.lp38 <- subset_samples(ps.all, group=="LP38")

ps.lp51 <- subset_samples(ps.all, group=="LP51")

ps.sm38 <- subset_samples(ps.all, group=="SM38")

ps.sm51 <- subset_samples(ps.all, group=="SM51")

#################################################################


# pairwise comparisons for all time points
#################################################################


# Larva:Control 38
ps.larv.c38 <- merge_phyloseq(ps.larv,ps.c38)
# make a data frame from the sample_data
data.larv.c38 <- data.frame(sample_data(ps.larv.c38))
Dist.larv.c38 = phyloseq::distance(ps.larv.c38, method = "bray", type="samples")

# betadisper: dispersion test
larv.c38.beta <- betadisper(Dist.larv.c38, data.larv.c38$treatment)
permutest(larv.c38.beta)
#            Df   Sum Sq  Mean Sq      F N.Perm Pr(>F)
# Groups     1 0.001359 0.001359 0.6833    999   0.44
# Residuals  7 0.013923 0.001989        

# Adonis test
adonis(Dist.larv.c38 ~ treatment, data = data.larv.c38)
#            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# treatment  1    1.6400 1.63999  32.111 0.82102  0.011 *
# Residuals  7    0.3575 0.05107         0.17898         
# Total      8    1.9975                 1.00000  


# Larva:Control 51
ps.larv.c51 <- merge_phyloseq(ps.larv,ps.c51)
# make a data frame from the sample_data
data.larv.c51 <- data.frame(sample_data(ps.larv.c51))
Dist.larv.c51 = phyloseq::distance(ps.larv.c51, method = "bray", type="samples")

# betadisper: dispersion test
larv.c51.beta <- betadisper(Dist.larv.c51, data.larv.c51$treatment)
permutest(larv.c51.beta)
#           Df   Sum Sq   Mean Sq     F N.Perm Pr(>F)
# Groups     1 0.003848 0.0038480 0.588    719 0.4014
# Residuals  4 0.026179 0.0065448            

# Adonis test
adonis(Dist.larv.c51 ~ treatment, data = data.larv.c51)
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# treatment  1    1.2639  1.2639  21.905 0.84559    0.1
# Residuals  4    0.2308  0.0577         0.15441       
# Total      5    1.4947                 1.00000  


# Larva:LP38
ps.larv.lp38 <- merge_phyloseq(ps.larv,ps.lp38)
# make a data frame from the sample_data
data.larv.lp38 <- data.frame(sample_data(ps.larv.lp38))
Dist.larv.lp38 = phyloseq::distance(ps.larv.lp38, method = "bray", type="samples")

# betadisper: dispersion test
larv.lp38.beta <- betadisper(Dist.larv.lp38, data.larv.lp38$treatment)
permutest(larv.lp38.beta)
#            Df    Sum Sq    Mean Sq      F N.Perm Pr(>F)
# Groups     1 0.0006144 0.00061439 0.2127    999  0.633
# Residuals  7 0.0202219 0.00288884             

# Adonis test
adonis(Dist.larv.lp38 ~ treatment, data = data.larv.lp38)
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
# treatment  1   1.59567 1.59567  23.619 0.77138  0.009 **
# Residuals  7   0.47292 0.06756         0.22862          
# Total      8   2.06859                 1.00000  


# Larva:LP51
ps.larv.lp51 <- merge_phyloseq(ps.larv,ps.lp51)
# make a data frame from the sample_data
data.larv.lp51 <- data.frame(sample_data(ps.larv.lp51))
Dist.larv.lp51 = phyloseq::distance(ps.larv.lp51, method = "bray", type="samples")

# betadisper: dispersion test
larv.lp51.beta <- betadisper(Dist.larv.lp51, data.larv.lp51$treatment)
permutest(larv.lp51.beta)
#           Df   Sum Sq   Mean Sq     F N.Perm Pr(>F)
# Groups     1 0.011755 0.0117548 4.026    719 0.1014
# Residuals  4 0.011679 0.0029197             

# Adonis test
adonis(Dist.larv.lp51 ~ treatment, data = data.larv.lp51)
#            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# treatment  1    1.2849 1.28493  26.952 0.87077    0.1
# Residuals  4    0.1907 0.04768         0.12923       
# Total      5    1.4756                 1.00000       


# Larva:SM38 
ps.larv.sm38 <- merge_phyloseq(ps.larv,ps.sm38)
# make a data frame from the sample_data
data.larv.sm38 <- data.frame(sample_data(ps.larv.sm38))
Dist.larv.sm38 = phyloseq::distance(ps.larv.sm38, method = "bray", type="samples")

# betadisper: dispersion test
larv.sm38.beta <- betadisper(Dist.larv.sm38, data.larv.sm38$treatment)
permutest(larv.sm38.beta)
#            Df    Sum Sq   Mean Sq      F N.Perm Pr(>F)
# Groups     1 0.0045646 0.0045646 1.8843    999  0.255
# Residuals  5 0.0121124 0.0024225             

# Adonis test
adonis(Dist.larv.sm38 ~ treatment, data = data.larv.sm38)
#            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# treatment  1   1.39143 1.39143   27.91 0.84807  0.021 *
# Residuals  5   0.24927 0.04985         0.15193         
# Total      6   1.64070                 1.00000       


# Larva:SM51
ps.larv.sm51 <- merge_phyloseq(ps.larv,ps.sm51)
# make a data frame from the sample_data
data.larv.sm51 <- data.frame(sample_data(ps.larv.sm51))
Dist.larv.sm51 = phyloseq::distance(ps.larv.sm51, method = "bray", type="samples")

# betadisper: dispersion test
larv.sm51.beta <- betadisper(Dist.larv.sm51, data.larv.sm51$treatment)
permutest(larv.sm51.beta)
#            Df    Sum Sq   Mean Sq      F N.Perm Pr(>F)
# Groups     1 0.0098148 0.0098148 3.0288    719 0.2014
# Residuals  4 0.0129618 0.0032405            

# Adonis test
adonis(Dist.larv.sm51 ~ treatment, data = data.larv.sm51)
#            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# treatment  1   1.29345 1.29345  26.227 0.86767    0.1
# Residuals  4   0.19727 0.04932         0.13233       
# Total      5   1.49072                 1.00000 


# Larva:SM51
ps.larv.sm51 <- merge_phyloseq(ps.larv,ps.sm51)
# make a data frame from the sample_data
data.larv.sm51 <- data.frame(sample_data(ps.larv.sm51))
Dist.larv.sm51 = phyloseq::distance(ps.larv.sm51, method = "bray", type="samples")

# betadisper: dispersion test
larv.sm51.beta <- betadisper(Dist.larv.sm51, data.larv.sm51$treatment)
permutest(larv.sm51.beta)
#            Df    Sum Sq   Mean Sq      F N.Perm Pr(>F)
# Groups     1 0.0098148 0.0098148 3.0288    719 0.2014
# Residuals  4 0.0129618 0.0032405            

# Adonis test
adonis(Dist.larv.sm51 ~ treatment, data = data.larv.sm51)
#            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# treatment  1   1.29345 1.29345  26.227 0.86767    0.1
# Residuals  4   0.19727 0.04932         0.13233       
# Total      5   1.49072                 1.00000 


# C38:C51
ps.c38.c51 <- merge_phyloseq(ps.c38,ps.c51)
# make a data frame from the sample_data
data.c38.c51 <- data.frame(sample_data(ps.c38.c51))
Dist.c38.c51 = phyloseq::distance(ps.c38.c51, method = "bray", type="samples")

# betadisper: dispersion test
c38.c51.beta <- betadisper(Dist.c38.c51, data.c38.c51$group)
permutest(c38.c51.beta)
#            Df    Sum Sq   Mean Sq     F N.Perm Pr(>F)
# Groups     1 0.0012086 0.0012086 0.492    999  0.494
# Residuals  7 0.0171944 0.0024563         

# Adonis test
adonis(Dist.c38.c51 ~ group, data = data.c38.c51)
#            Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)   
# group      1   0.27857 0.278570  6.4619 0.48002  0.007 **
# Residuals  7   0.30177 0.043109         0.51998          
# Total      8   0.58034                  1.00000  


# C38:C51
ps.c38.lp38 <- merge_phyloseq(ps.c38,ps.lp38)
# make a data frame from the sample_data
data.c38.lp38 <- data.frame(sample_data(ps.c38.lp38))
Dist.c38.lp38 = phyloseq::distance(ps.c38.lp38, method = "bray", type="samples")

# betadisper: dispersion test
c38.lp38.beta <- betadisper(Dist.c38.lp38, data.c38.lp38$group)
permutest(c38.lp38.beta)
#            Df    Sum Sq   Mean Sq      F N.Perm Pr(>F)  
# Groups     1 0.0057013 0.0057013 5.0745    999  0.048 *
# Residuals 10 0.0112352 0.0011235         

# Adonis test
adonis(Dist.c38.lp38 ~ group, data = data.c38.lp38)
#            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
# group      1   0.61271 0.61271  11.265 0.52975  0.003 **
# Residuals 10   0.54389 0.05439         0.47025          
# Total     11   1.15660                 1.00000  


# C38:LP51
ps.c38.lp51 <- merge_phyloseq(ps.c38,ps.lp51)
# make a data frame from the sample_data
data.c38.lp51 <- data.frame(sample_data(ps.c38.lp51))
Dist.c38.lp51 = phyloseq::distance(ps.c38.lp51, method = "bray", type="samples")

# betadisper: dispersion test
c38.lp51.beta <- betadisper(Dist.c38.lp51, data.c38.lp51$group)
permutest(c38.lp51.beta)
#            Df    Sum Sq   Mean Sq      F N.Perm Pr(>F)    
# Groups     1 0.0078017 0.0078017 20.271    999  0.001 ***
# Residuals  7 0.0026941 0.0003849   

# Adonis test
adonis(Dist.c38.lp51 ~ group, data = data.c38.lp51)
#            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# group      1   0.57392 0.57392  15.353 0.68684  0.013 *
# Residuals  7   0.26167 0.03738         0.31316         
# Total      8   0.83559                 1.00000  


# C38:SM38
ps.c38.sm38 <- merge_phyloseq(ps.c38,ps.sm38)
# make a data frame from the sample_data
data.c38.sm38 <- data.frame(sample_data(ps.c38.sm38))
Dist.c38.sm38 = phyloseq::distance(ps.c38.sm38, method = "bray", type="samples")

# betadisper: dispersion test
c38.sm38.beta <- betadisper(Dist.c38.sm38, data.c38.sm38$group)
permutest(c38.sm38.beta)
#            Df    Sum Sq    Mean Sq      F N.Perm Pr(>F)  
# Groups     1 0.0015648 0.00156480 4.0024    999  0.066 .
# Residuals  8 0.0031277 0.00039097   

# Adonis test
adonis(Dist.c38.sm38 ~ group, data = data.c38.sm38)
#            Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)   
# group      1   0.25207 0.252071  6.2971 0.44045  0.005 **
# Residuals  8   0.32024 0.040029         0.55955          
# Total      9   0.57231                  1.00000 


# C38:SM51
ps.c38.sm51 <- merge_phyloseq(ps.c38,ps.sm51)
# make a data frame from the sample_data
data.c38.sm51 <- data.frame(sample_data(ps.c38.sm51))
Dist.c38.sm51 = phyloseq::distance(ps.c38.sm51, method = "bray", type="samples")

# betadisper: dispersion test
c38.sm51.beta <- betadisper(Dist.c38.sm51, data.c38.sm51$group)
permutest(c38.sm51.beta)
#            Df    Sum Sq   Mean Sq     F N.Perm Pr(>F)  
# Groups     1 0.0060111 0.0060111 10.58    999  0.017 *
# Residuals  7 0.0039771 0.0005682 

# Adonis test
adonis(Dist.c38.sm51 ~ group, data = data.c38.sm51)
#            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
# group      1   0.34591 0.34591  9.0269 0.56323   0.01 **
# Residuals  7   0.26824 0.03832         0.43677          
# Total      8   0.61415                 1.00000  


# C51:LP38 
ps.c51.lp38 <- merge_phyloseq(ps.c51,ps.lp38)
# make a data frame from the sample_data
data.c51.lp38 <- data.frame(sample_data(ps.c51.lp38))
Dist.c51.lp38 = phyloseq::distance(ps.c51.lp38, method = "bray", type="samples")

# betadisper: dispersion test
c51.lp38.beta <- betadisper(Dist.c51.lp38, data.c51.lp38$group)
permutest(c51.lp38.beta)
#           Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
# Groups     1 0.009296 0.0092960 2.7697    999  0.111
# Residuals  7 0.023494 0.0033563   

# Adonis test
adonis(Dist.c51.lp38 ~ group, data = data.c51.lp38)
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# group      1   0.40020  0.4002  6.7151 0.48961  0.016 *
# Residuals  7   0.41718  0.0596         0.51039         
# Total      8   0.81738                 1.00000  


# C51:LP51
ps.c51.lp51 <- merge_phyloseq(ps.c51,ps.lp51)
# make a data frame from the sample_data
data.c51.lp51 <- data.frame(sample_data(ps.c51.lp51))
Dist.c51.lp51 = phyloseq::distance(ps.c51.lp51, method = "bray", type="samples")

# betadisper: dispersion test
c51.lp51.beta <- betadisper(Dist.c51.lp51, data.c51.lp51$group)
permutest(c51.lp51.beta)
#           Df    Sum Sq   Mean Sq      F N.Perm Pr(>F)
# Groups     1 0.0021518 0.0021518 0.5757    719 0.6014
# Residuals  4 0.0149509 0.0037377  

# Adonis test
adonis(Dist.c51.lp51 ~ group, data = data.c51.lp51)
#            Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
# group      1   0.29558 0.295581  8.7603 0.68653    0.1
# Residuals  4   0.13496 0.033741         0.31347       
# Total      5   0.43055                  1.00000       


# C51:SM51
ps.c51.sm51 <- merge_phyloseq(ps.c51,ps.sm51)
# make a data frame from the sample_data
data.c51.sm51 <- data.frame(sample_data(ps.c51.sm51))
Dist.c51.sm51 = phyloseq::distance(ps.c51.sm51, method = "bray", type="samples")

# betadisper: dispersion test
c51.sm51.beta <- betadisper(Dist.c51.sm51, data.c51.sm51$group)
permutest(c51.sm51.beta)
#           Df    Sum Sq   Mean Sq     F N.Perm Pr(>F)
# Groups     1 0.0013718 0.0013718 0.338    719 0.7014
# Residuals  4 0.0162339 0.0040585

# Adonis test
adonis(Dist.c51.sm51 ~ group, data = data.c51.sm51)
#           Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
# group      1   0.12766 0.127661   3.608 0.47423    0.1
# Residuals  4   0.14153 0.035383         0.52577       
# Total      5   0.26919                  1.00000 


# LP38:LP51 
ps.lp38.lp51 <- merge_phyloseq(ps.lp38,ps.lp51)
# make a data frame from the sample_data
data.lp38.lp51 <- data.frame(sample_data(ps.lp38.lp51))
Dist.lp38.lp51 = phyloseq::distance(ps.lp38.lp51, method = "bray", type="samples")

# betadisper: dispersion test
lp38.lp51.beta <- betadisper(Dist.lp38.lp51, data.lp38.lp51$group)
permutest(lp38.lp51.beta)
#           Df    Sum Sq   Mean Sq      F N.Perm Pr(>F)    
# Groups     1 0.0224936 0.0224936 17.507    999  0.001 ***
# Residuals  7 0.0089938 0.0012848 

# Adonis test
adonis(Dist.lp38.lp51 ~ group, data = data.lp38.lp51)
#           Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
# group      1   0.26230 0.262299  4.8692 0.41024  0.014 *
# Residuals  7   0.37708 0.053869         0.58976         
# Total      8   0.63938                  1.00000 


# LP38:SM38
ps.lp38.sm38 <- merge_phyloseq(ps.lp38,ps.sm38)
# make a data frame from the sample_data
data.lp38.sm38 <- data.frame(sample_data(ps.lp38.sm38))
Dist.lp38.sm38 = phyloseq::distance(ps.lp38.sm38, method = "bray", type="samples")

# betadisper: dispersion test
lp38.sm38.beta <- betadisper(Dist.lp38.sm38, data.lp38.sm38$group)
permutest(lp38.sm38.beta)
#            Df    Sum Sq   Mean Sq      F N.Perm Pr(>F)    
# Groups     1 0.0114690 0.0114690 9.7324    999  0.001 ***
# Residuals  8 0.0094274 0.0011784 

# Adonis test
adonis(Dist.lp38.sm38 ~ group, data = data.lp38.sm38)
#            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
# group      1   0.36470 0.36470  6.6971 0.45568  0.004 **
# Residuals  8   0.43565 0.05446         0.54432          
# Total      9   0.80035                 1.00000   


# LP38:SM51 
ps.lp38.sm51 <- merge_phyloseq(ps.lp38,ps.sm51)
# make a data frame from the sample_data
data.lp38.sm51 <- data.frame(sample_data(ps.lp38.sm51))
Dist.lp38.sm51 = phyloseq::distance(ps.lp38.sm51, method = "bray", type="samples")

# betadisper: dispersion test
lp38.sm51.beta <- betadisper(Dist.lp38.sm51, data.lp38.sm51$group)
permutest(lp38.sm51.beta)
#           Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)    
# Groups     1 0.019372 0.0193719 13.195    999  0.001 ***
# Residuals  7 0.010277 0.0014681  

# Adonis test
adonis(Dist.lp38.sm51 ~ group, data = data.lp38.sm51)
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# group      1   0.45357 0.45357  8.2757 0.54176  0.012 *
# Residuals  7   0.38365 0.05481         0.45824         
# Total      8   0.83722                 1.00000  


# LP51:SM38 
ps.lp51.sm38 <- merge_phyloseq(ps.lp51,ps.sm38)
# make a data frame from the sample_data
data.lp51.sm38 <- data.frame(sample_data(ps.lp51.sm38))
Dist.lp51.sm38 = phyloseq::distance(ps.lp51.sm38, method = "bray", type="samples")

# betadisper: dispersion test
lp51.sm38.beta <- betadisper(Dist.lp51.sm38, data.lp51.sm38$group)
permutest(lp51.sm38.beta)
#            Df     Sum Sq    Mean Sq      F N.Perm Pr(>F)  
# Groups     1 0.00233705 0.00233705 13.216    999  0.032 *
# Residuals  5 0.00088421 0.00017684  

# Adonis test
adonis(Dist.lp51.sm38 ~ group, data = data.lp51.sm38)
#            Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
# group      1   0.42219 0.42219  13.758 0.73345  0.028 *
# Residuals  5   0.15343 0.03069         0.26655         
# Total      6   0.57562                 1.00000 


# LP51:SM51 
ps.lp51.sm51 <- merge_phyloseq(ps.lp51,ps.sm51)
# make a data frame from the sample_data
data.lp51.sm51 <- data.frame(sample_data(ps.lp51.sm51))
Dist.lp51.sm51 = phyloseq::distance(ps.lp51.sm51, method = "bray", type="samples")

# betadisper: dispersion test
lp51.sm51.beta <- betadisper(Dist.lp51.sm51, data.lp51.sm51$group)
permutest(lp51.sm51.beta)
#            Df     Sum Sq   Mean Sq      F N.Perm Pr(>F)
# Groups     1 0.00008741 8.741e-05 0.2017    719 0.7014
# Residuals  4 0.00173362 4.334e-04   

# Adonis test
adonis(Dist.lp51.sm51 ~ group, data = data.lp51.sm51)
#            Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)
# group      1   0.29045 0.290447  11.453 0.74115    0.1
# Residuals  4   0.10144 0.025359         0.25885       
# Total      5   0.39188                  1.00000  


# SM38:SM51  
ps.sm38.sm51 <- merge_phyloseq(ps.sm38,ps.sm51)
# make a data frame from the sample_data
data.sm38.sm51 <- data.frame(sample_data(ps.sm38.sm51))
Dist.sm38.sm51 = phyloseq::distance(ps.sm38.sm51, method = "bray", type="samples")

# betadisper: dispersion test
sm38.sm51.beta <- betadisper(Dist.sm38.sm51, data.sm38.sm51$group)
permutest(sm38.sm51.beta)
#            Df    Sum Sq    Mean Sq      F N.Perm Pr(>F)
# Groups     1 0.0014706 0.00147057 3.3928    999  0.159
# Residuals  5 0.0021672 0.00043344   

# Adonis test
adonis(Dist.sm38.sm51 ~ group, data = data.sm38.sm51)
#            Df SumsOfSqs  MeanSqs F.Model      R2 Pr(>F)  
# group      1   0.28257 0.282573  8.8303 0.63847  0.041 *
# Residuals  5   0.16000 0.032001         0.36153         
# Total      6   0.44258                  1.00000  




# applying bonferrni correction to p values from both tests
# a
# levene
levene.p <- c(0.44, 0.4014, 0.633, 0.1014, 0.255, 0.2014, 0.494, 0.048, 0.001, 0.066, 0.017, 0.111, 0.6014, 0.971,
              0.7014, 0.001, 0.001, 0.001, 0.032, 0.7014, 0.159) 

levene.p.bonf <- p.adjust(levene.p, method="bonferroni")
levene.p.bonf

# permanova
permanova.p <- c(0.011, 0.1, 0.009, 0.1, 0.021, 0.1, 0.007, 0.003, 0.013, 0.005, 0.01, 0.016, 0.1, 0.032, 0.1, 0.014, 
                 0.004, 0.012, 0.028, 0.1, 0.041) 

permanova.p.bonf <- p.adjust(permanova.p, method="bonferroni")
permanova.p.bonf


# b
# levene
levene.p <- c(0.44, 0.4014, 0.633, 0.1014, 0.255, 0.2014, 0.494, 0.048, 0.066, 0.6014, 0.7014, 0.001, 0.159) 

levene.p.bonf <- p.adjust(levene.p, method="bonferroni")
levene.p.bonf

# permanova
permanova.p <- c(0.011, 0.1, 0.009, 0.1, 0.021, 0.1, 0.007, 0.003, 0.005, 0.1, 0.1, 0.014, 0.041) 

permanova.p.bonf <- p.adjust(permanova.p, method="bonferroni")
permanova.p.bonf




#################################################################








# Taxonomy tests

# All Treatments
# Comparison on means, KW, Dunn 
#################################################################

#define standard error function
se <- function(x) sqrt(var(x)/length(x))

#collapse ASVs at family level
AvgRA_o = tax_glom(ps.all, taxrank="Family", NArm = FALSE)

#filter out taxa that don't vary > 0.2 %
AvgRA99_O = filter_taxa(AvgRA_o, function(x) var(x) > .002, TRUE)
df_o <- psmelt(AvgRA99_O)
#group and calculate mean, sd and se for different taxonomic levels
grouped_o <- group_by(df_o, treatment, Family, Order, Class, Phylum)
avgs_o <- summarise(grouped_o, mean=mean(Abundance), sd=sd(Abundance), se=se(Abundance))

avgs_filt <- filter(avgs_o, mean > 0)
# I want to read out this table for sup data
write.csv(avgs_filt, "treatment_fam_avg_sd_se.csv")


avgs_filt$ST <- as.factor(avgs_filt$treatment)
avgs_filt_st <- melt(data.frame(avgs_filt), id.vars=c("Family", "treatment"),
                  measure.vars=c("mean"))
kruskal_test(mean ~ ST, distribution=approximate(nresample=9999), data=avgs_filt)

for (i in 1:length(levels(df_o$Family))) {
  new_df <- subset(df_o, df_o$Family == levels(df_o$Family)[i])
  new_df$ST <- as.factor(new_df$treatment)
  print(new_df$Genus[1])
  print(kruskal_test(Abundance ~ ST, distribution=approximate(nresample=9999), data=new_df))
  print( dunnTest(Abundance ~ ST, data=new_df, method ="bonferroni"))}

#start
chisq = NULL
pvals = NULL
listofcats_sig = NULL
DataSet = df_o
taxa = as.factor(unique(df_o$Family))


for (cat in taxa) {
  new_df <- subset(DataSet, DataSet$Family == cat)
  new_df$ST <- as.factor(new_df$treatment)
  kw = kruskal.test(Abundance ~ ST, data=new_df)
  chisq = c(kw$statistic, chisq)
  pvals = c(kw$p.value, pvals)
  if (kw$p.value <= 0.05) {
    listofcats_sig = c(cat, listofcats_sig)
  }
}

pvals.bonf = p.adjust(pvals, method="bonferroni")
df_taxa = data.frame(rev(taxa), chisq, pvals, pvals.bonf)
write.table(df_taxa, 'Treatment_VST_Mean_FAM_KW.txt', sep="\t")

pvals.dunn = NULL
pvals.dunn.bonf = NULL
Zsc.dunn = NULL
comparison.dunn = NULL
cats.dunn = NULL

for (cat in listofcats_sig){
  new_df <- subset(DataSet, DataSet$Family == cat)
  new_df$ST <- as.factor(new_df$treatment)
  dT = dunnTest(Abundance ~ ST, data = new_df, method = "bonferroni")
  for (i in 1:length(dT$res$Comparison)) {
    print(cat)
    pvals.dunn = c(dT$res$P.unadj[i], pvals.dunn)
    pvals.dunn.bonf = c(dT$res$P.adj[i], pvals.dunn.bonf)
    Zsc.dunn = c(dT$res$Z[i], Zsc.dunn)
    comparison.dunn = c(levels(dT$res$Comparison)[i], comparison.dunn)
    cats.dunn = c(cat, cats.dunn)
  }
}

df_taxa = data.frame(cats.dunn, comparison.dunn, Zsc.dunn, pvals.dunn, pvals.dunn.bonf)
tax_oi <- filter(df_taxa, pvals.dunn.bonf < 0.05)

write.csv(df_taxa, 'Treatment_VST_Mean_FAM_KW_Dunn_OI.csv')


#################################################################





# Control 38:52
# Comparison on means, KW, Dunn 
#################################################################

# ps.control is a phyloseq object with only control samples
# was first subset in the Ordination statistics section

#collapse ASVs at family level
AvgRA_o = tax_glom(ps.larv.cntrl, taxrank="Family", NArm = FALSE)

#filter out taxa that don't vary > 0.2 %
AvgRA99_O = filter_taxa(AvgRA_o, function(x) var(x) > .002, TRUE)
df_o <- psmelt(AvgRA99_O)
#group and calculate mean, sd and se for different taxonomic levels
grouped_o <- group_by(df_o, age, Family, Order, Class, Phylum)
#avgs_o <- summarise(grouped_o, mean=100*mean(Abundance), sd=100*sd(Abundance), se=100*se(Abundance))
avgs_o <- summarise(grouped_o, mean=mean(Abundance), sd=sd(Abundance), se=se(Abundance))

write.csv(avgs_o, 'Larv-Cntrl-Age/L-Cntrl_Age_Mean_sd_se.csv')


avgs_o$ST <- as.factor(avgs_o$age)
# avgs_o_st <- melt(data.frame(avgs_o), id.vars=c("Order", "Substrate"),
#                   measure.vars=c("mean"))
kruskal_test(mean ~ ST, distribution=approximate(nresample=9999), data=avgs_o)

for (i in 1:length(levels(df_o$Family))) {
  new_df <- subset(df_o, df_o$Family == levels(df_o$Family)[i])
  new_df$ST <- as.factor(new_df$age) # make sure to change the variable being tested
  print(new_df$Order[1])
  print(kruskal_test(Abundance ~ ST, distribution=approximate(nresample=9999), data=new_df))
  print( dunnTest(Abundance ~ ST, data=new_df, method ="bonferroni"))}

#start
chisq = NULL
pvals = NULL
listofcats_sig = NULL
DataSet = df_o
taxa = unique(df_o$Family)

for (cat in taxa) {
  new_df <- subset(DataSet, DataSet$Family == cat)
  new_df$ST <- as.factor(new_df$age)
  kw = kruskal.test(Abundance ~ ST, data=new_df)
  chisq = c(kw$statistic, chisq)
  pvals = c(kw$p.value, pvals)
  if (kw$p.value <= 0.05) {
    16
    listofcats_sig = c(cat, listofcats_sig)
  }
}

pvals.bonf = p.adjust(pvals, method="bonferroni")
df_taxa = data.frame(rev(taxa), chisq, pvals, pvals.bonf)
write.csv(df_taxa, 'Larv-Cntrl-Age/L-Cntrl_Age_Mean_Fam_KW_002var.csv')

pvals.dunn = NULL
pvals.dunn.bonf = NULL
Zsc.dunn = NULL
comparison.dunn = NULL
cats.dunn = NULL

for (cat in listofcats_sig){
  new_df <- subset(DataSet, DataSet$Family == cat)
  new_df$ST <- as.factor(new_df$age)
  dT = dunnTest(Abundance ~ ST, data = new_df, method = "bonferroni")
  for (i in 1:length(dT$res$Comparison)) {
    print(cat)
    pvals.dunn = c(dT$res$P.unadj[i], pvals.dunn)
    pvals.dunn.bonf = c(dT$res$P.adj[i], pvals.dunn.bonf)
    Zsc.dunn = c(dT$res$Z[i], Zsc.dunn)
    comparison.dunn = c(levels(dT$res$Comparison)[i], comparison.dunn)
    cats.dunn = c(cat, cats.dunn)
  }
}

#
df_taxa = data.frame(cats.dunn, comparison.dunn, Zsc.dunn, pvals.dunn, pvals.dunn.bonf)
write.csv(df_taxa, 'Larv-Cntrl-Age/L-Cntrl_Age_Mean_Order_KW_002var_Dunn.csv')

#################################################################



# Low-pH 38:52
# Comparison on means, KW, Dunn 
#################################################################

# ps.larv.oa is a phyloseq object with only control samples
# was first subset in the Ordination statistics section

#collapse ASVs at family level
AvgRA_o = tax_glom(ps.larv.oa, taxrank="Family", NArm = FALSE)

#filter out taxa that don't vary > 0.2 %
AvgRA99_O = filter_taxa(AvgRA_o, function(x) var(x) > .002, TRUE)
df_o <- psmelt(AvgRA99_O)
#group and calculate mean, sd and se for different taxonomic levels
grouped_o <- group_by(df_o, age, Family, Order, Class, Phylum)
#avgs_o <- summarise(grouped_o, mean=100*mean(Abundance), sd=100*sd(Abundance), se=100*se(Abundance))
avgs_o <- summarise(grouped_o, mean=mean(Abundance), sd=sd(Abundance), se=se(Abundance))

write.csv(avgs_o, 'Larv-pH-Age/L-pH_Age_Mean_sd_se.csv')


avgs_o$ST <- as.factor(avgs_o$age)
# avgs_o_st <- melt(data.frame(avgs_o), id.vars=c("Order", "Substrate"),
#                   measure.vars=c("mean"))
kruskal_test(mean ~ ST, distribution=approximate(nresample=9999), data=avgs_o)

for (i in 1:length(levels(df_o$Family))) {
  new_df <- subset(df_o, df_o$Family == levels(df_o$Family)[i])
  new_df$ST <- as.factor(new_df$age) # make sure to change the variable being tested
  print(new_df$Order[1])
  print(kruskal_test(Abundance ~ ST, distribution=approximate(nresample=9999), data=new_df))
  print( dunnTest(Abundance ~ ST, data=new_df, method ="bonferroni"))}

#start
chisq = NULL
pvals = NULL
listofcats_sig = NULL
DataSet = df_o
taxa = unique(df_o$Family)

for (cat in taxa) {
  new_df <- subset(DataSet, DataSet$Family == cat)
  new_df$ST <- as.factor(new_df$age)
  kw = kruskal.test(Abundance ~ ST, data=new_df)
  chisq = c(kw$statistic, chisq)
  pvals = c(kw$p.value, pvals)
  if (kw$p.value <= 0.05) {
    16
    listofcats_sig = c(cat, listofcats_sig)
  }
}

pvals.bonf = p.adjust(pvals, method="bonferroni")
df_taxa = data.frame(rev(taxa), chisq, pvals, pvals.bonf)
write.csv(df_taxa, 'Larv-pH-Age/L-pH_Age_Mean_Fam_KW_002var.csv')

pvals.dunn = NULL
pvals.dunn.bonf = NULL
Zsc.dunn = NULL
comparison.dunn = NULL
cats.dunn = NULL

for (cat in listofcats_sig){
  new_df <- subset(DataSet, DataSet$Family == cat)
  new_df$ST <- as.factor(new_df$age)
  dT = dunnTest(Abundance ~ ST, data = new_df, method = "bonferroni")
  for (i in 1:length(dT$res$Comparison)) {
    print(cat)
    pvals.dunn = c(dT$res$P.unadj[i], pvals.dunn)
    pvals.dunn.bonf = c(dT$res$P.adj[i], pvals.dunn.bonf)
    Zsc.dunn = c(dT$res$Z[i], Zsc.dunn)
    comparison.dunn = c(levels(dT$res$Comparison)[i], comparison.dunn)
    cats.dunn = c(cat, cats.dunn)
  }
}

#
df_taxa = data.frame(cats.dunn, comparison.dunn, Zsc.dunn, pvals.dunn, pvals.dunn.bonf)
write.csv(df_taxa, 'Larv-pH-Age/L-pH_Age_Mean_Order_KW_002var_Dunn.csv')

#################################################################




# Molybdate 38:52
# Comparison on means, KW, Dunn 
#################################################################

# ps.arv.sm is a phyloseq object with only control samples
# was first subset in the Ordination statistics section

#collapse ASVs at family level
AvgRA_o = tax_glom(ps.larv.sm, taxrank="Family", NArm = FALSE)

#filter out taxa that don't vary > 0.2 %
AvgRA99_O = filter_taxa(AvgRA_o, function(x) var(x) > .002, TRUE)
df_o <- psmelt(AvgRA99_O)
#group and calculate mean, sd and se for different taxonomic levels
grouped_o <- group_by(df_o, age, Family, Order, Class, Phylum)
#avgs_o <- summarise(grouped_o, mean=100*mean(Abundance), sd=100*sd(Abundance), se=100*se(Abundance))
avgs_o <- summarise(grouped_o, mean=mean(Abundance), sd=sd(Abundance), se=se(Abundance))

write.csv(avgs_o, 'Larv-SM-Age/L-SM_Age_Mean_sd_se.csv')


avgs_o$ST <- as.factor(avgs_o$age)
# avgs_o_st <- melt(data.frame(avgs_o), id.vars=c("Order", "Substrate"),
#                   measure.vars=c("mean"))
kruskal_test(mean ~ ST, distribution=approximate(nresample=9999), data=avgs_o)

for (i in 1:length(levels(df_o$Family))) {
  new_df <- subset(df_o, df_o$Family == levels(df_o$Family)[i])
  new_df$ST <- as.factor(new_df$age) # make sure to change the variable being tested
  print(new_df$Order[1])
  print(kruskal_test(Abundance ~ ST, distribution=approximate(nresample=9999), data=new_df))
  print( dunnTest(Abundance ~ ST, data=new_df, method ="bonferroni"))}

#start
chisq = NULL
pvals = NULL
listofcats_sig = NULL
DataSet = df_o
taxa = unique(df_o$Family)

for (cat in taxa) {
  new_df <- subset(DataSet, DataSet$Family == cat)
  new_df$ST <- as.factor(new_df$age)
  kw = kruskal.test(Abundance ~ ST, data=new_df)
  chisq = c(kw$statistic, chisq)
  pvals = c(kw$p.value, pvals)
  if (kw$p.value <= 0.05) {
    16
    listofcats_sig = c(cat, listofcats_sig)
  }
}

pvals.bonf = p.adjust(pvals, method="bonferroni")
df_taxa = data.frame(rev(taxa), chisq, pvals, pvals.bonf)
write.csv(df_taxa, 'Larv-SM-Age/L-SM_Age_Mean_Fam_KW_002var.csv')

pvals.dunn = NULL
pvals.dunn.bonf = NULL
Zsc.dunn = NULL
comparison.dunn = NULL
cats.dunn = NULL

for (cat in listofcats_sig){
  new_df <- subset(DataSet, DataSet$Family == cat)
  new_df$ST <- as.factor(new_df$age)
  dT = dunnTest(Abundance ~ ST, data = new_df, method = "bonferroni")
  for (i in 1:length(dT$res$Comparison)) {
    print(cat)
    pvals.dunn = c(dT$res$P.unadj[i], pvals.dunn)
    pvals.dunn.bonf = c(dT$res$P.adj[i], pvals.dunn.bonf)
    Zsc.dunn = c(dT$res$Z[i], Zsc.dunn)
    comparison.dunn = c(levels(dT$res$Comparison)[i], comparison.dunn)
    cats.dunn = c(cat, cats.dunn)
  }
}

#
df_taxa = data.frame(cats.dunn, comparison.dunn, Zsc.dunn, pvals.dunn, pvals.dunn.bonf)
write.csv(df_taxa, 'Larv-SM-Age/L-SM_Age_Mean_Order_KW_002var_Dunn.csv')

#################################################################




# Control : Low-pH
# Comparison on means, KW, Dunn 
#################################################################

var <- paste( sample_data(ps.cntrl.oa)$treatment, sample_data(ps.cntrl.oa)$age, sep = "")

sample_data(ps.cntrl.oa)$treatage <- var


#define standard error function
se <- function(x) sqrt(var(x)/length(x))

#collapse ASVs at family level
AvgRA_o = tax_glom(ps.cntrl.oa, taxrank="Family", NArm = FALSE)

#filter out taxa that don't vary > 0.2 %
AvgRA99_O = filter_taxa(AvgRA_o, function(x) var(x) > .002, TRUE)
df_o <- psmelt(AvgRA99_O)
#group and calculate mean, sd and se for different taxonomic levels
grouped_o <- group_by(df_o, treatage, Family, Order, Class, Phylum)
avgs_o <- summarise(grouped_o, mean=mean(Abundance), sd=sd(Abundance), se=se(Abundance))

avgs_filt <- filter(avgs_o, mean > 0)
# I want to read out this table for sup data
write.csv(avgs_filt, "cntrl-lp_fam_avg_sd_se .csv")


avgs_filt$ST <- as.factor(avgs_filt$treatage)
avgs_filt_st <- melt(data.frame(avgs_filt), id.vars=c("Family", "treatage"),
                     measure.vars=c("mean"))
kruskal_test(mean ~ ST, distribution=approximate(nresample=9999), data=avgs_filt)

for (i in 1:length(levels(df_o$Family))) {
  new_df <- subset(df_o, df_o$Family == levels(df_o$Family)[i])
  new_df$ST <- as.factor(new_df$treatage)
  print(new_df$Genus[1])
  print(kruskal_test(Abundance ~ ST, distribution=approximate(nresample=9999), data=new_df))
  print( dunnTest(Abundance ~ ST, data=new_df, method ="bonferroni"))}

#start
chisq = NULL
pvals = NULL
listofcats_sig = NULL
DataSet = df_o
taxa = as.factor(unique(df_o$Family))


for (cat in taxa) {
  new_df <- subset(DataSet, DataSet$Family == cat)
  new_df$ST <- as.factor(new_df$treatage)
  kw = kruskal.test(Abundance ~ ST, data=new_df)
  chisq = c(kw$statistic, chisq)
  pvals = c(kw$p.value, pvals)
  if (kw$p.value <= 0.05) {
    listofcats_sig = c(cat, listofcats_sig)
  }
}

pvals.bonf = p.adjust(pvals, method="bonferroni")
df_taxa = data.frame(rev(taxa), chisq, pvals, pvals.bonf)
write.csv(df_taxa, 'Cntrl-LP_VST_Mean_FAM_KW.csv')

pvals.dunn = NULL
pvals.dunn.bonf = NULL
Zsc.dunn = NULL
comparison.dunn = NULL
cats.dunn = NULL

for (cat in listofcats_sig){
  new_df <- subset(DataSet, DataSet$Family == cat)
  new_df$ST <- as.factor(new_df$treatage)
  dT = dunnTest(Abundance ~ ST, data = new_df, method = "bonferroni")
  for (i in 1:length(dT$res$Comparison)) {
    print(cat)
    pvals.dunn = c(dT$res$P.unadj[i], pvals.dunn)
    pvals.dunn.bonf = c(dT$res$P.adj[i], pvals.dunn.bonf)
    Zsc.dunn = c(dT$res$Z[i], Zsc.dunn)
    comparison.dunn = c(levels(dT$res$Comparison)[i], comparison.dunn)
    cats.dunn = c(cat, cats.dunn)
  }
}

df_taxa = data.frame(cats.dunn, comparison.dunn, Zsc.dunn, pvals.dunn, pvals.dunn.bonf)
tax_oi <- filter(df_taxa, pvals.dunn.bonf < 0.05)

write.csv(df_taxa, 'Cntrl-LP_VST_Mean_FAM_KW_Dunn_OI.csv')

#################################################################




# Control : Sodium Molybdate
# Comparison on means, KW, Dunn 
#################################################################

var <- paste( sample_data(ps.cntrl.sm)$treatment, sample_data(ps.cntrl.sm)$age, sep = "")

sample_data(ps.cntrl.sm)$treatage <- var


#define standard error function
se <- function(x) sqrt(var(x)/length(x))

#collapse ASVs at family level
AvgRA_o = tax_glom(ps.cntrl.sm, taxrank="Family", NArm = FALSE)

#filter out taxa that don't vary > 0.2 %
AvgRA99_O = filter_taxa(AvgRA_o, function(x) var(x) > .002, TRUE)
df_o <- psmelt(AvgRA99_O)
#group and calculate mean, sd and se for different taxonomic levels
grouped_o <- group_by(df_o, treatage, Family, Order, Class, Phylum)
avgs_o <- summarise(grouped_o, mean=mean(Abundance), sd=sd(Abundance), se=se(Abundance))

avgs_filt <- filter(avgs_o, mean > 0)
# I want to read out this table for sup data
write.csv(avgs_filt, "cntrl-sm_fam_avg_sd_se .csv")


avgs_filt$ST <- as.factor(avgs_filt$treatage)
avgs_filt_st <- melt(data.frame(avgs_filt), id.vars=c("Family", "treatage"),
                     measure.vars=c("mean"))
kruskal_test(mean ~ ST, distribution=approximate(nresample=9999), data=avgs_filt)

for (i in 1:length(levels(df_o$Family))) {
  new_df <- subset(df_o, df_o$Family == levels(df_o$Family)[i])
  new_df$ST <- as.factor(new_df$treatage)
  print(new_df$Genus[1])
  print(kruskal_test(Abundance ~ ST, distribution=approximate(nresample=9999), data=new_df))
  print( dunnTest(Abundance ~ ST, data=new_df, method ="bonferroni"))}

#start
chisq = NULL
pvals = NULL
listofcats_sig = NULL
DataSet = df_o
taxa = as.factor(unique(df_o$Family))


for (cat in taxa) {
  new_df <- subset(DataSet, DataSet$Family == cat)
  new_df$ST <- as.factor(new_df$treatage)
  kw = kruskal.test(Abundance ~ ST, data=new_df)
  chisq = c(kw$statistic, chisq)
  pvals = c(kw$p.value, pvals)
  if (kw$p.value <= 0.05) {
    listofcats_sig = c(cat, listofcats_sig)
  }
}

pvals.bonf = p.adjust(pvals, method="bonferroni")
df_taxa = data.frame(rev(taxa), chisq, pvals, pvals.bonf)
write.csv(df_taxa, 'Cntrl-SM_VST_Mean_FAM_KW.csv')

pvals.dunn = NULL
pvals.dunn.bonf = NULL
Zsc.dunn = NULL
comparison.dunn = NULL
cats.dunn = NULL

for (cat in listofcats_sig){
  new_df <- subset(DataSet, DataSet$Family == cat)
  new_df$ST <- as.factor(new_df$treatage)
  dT = dunnTest(Abundance ~ ST, data = new_df, method = "bonferroni")
  for (i in 1:length(dT$res$Comparison)) {
    print(cat)
    pvals.dunn = c(dT$res$P.unadj[i], pvals.dunn)
    pvals.dunn.bonf = c(dT$res$P.adj[i], pvals.dunn.bonf)
    Zsc.dunn = c(dT$res$Z[i], Zsc.dunn)
    comparison.dunn = c(levels(dT$res$Comparison)[i], comparison.dunn)
    cats.dunn = c(cat, cats.dunn)
  }
}

df_taxa = data.frame(cats.dunn, comparison.dunn, Zsc.dunn, pvals.dunn, pvals.dunn.bonf)
tax_oi <- filter(df_taxa, pvals.dunn.bonf < 0.05)

write.csv(df_taxa, 'Cntrl-SM_VST_Mean_FAM_KW_Dunn_OI.csv')

#################################################################






# Taxonomic Bar Charts
#################################################################


# new df, change families
avgs_f <- avgs_filt
avgs_f$Family <- factor(avgs_f$Family,
                        levels=families)
avgs_f <- na.omit(avgs_f)
# remove families with %<5
avgs_5 <- filter(avgs_f, mean > 5)
avgs_3 <- filter(avgs_f, mean > 3)
avgs_1 <- filter(avgs_f, mean > 1)

# New facet label names for treatment variable
ad_labeller <- as_labeller( c("control" = "Control",
                              "larvae" = "Larva",
                              "low-ph" = "Low-pH",
                              "molybdate" = "Sodium Molybdate") )

# 12 colors
col_12 <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928')

alpha.ids <- c("C38-1",  "C38-2",  "C51",  "L",  
               "LP38-1", "LP28-2", "LP51", 
               "SM38-1", "SM38-2", "SM51")

# All treatments
ggplot(avgs_3, aes(x=bucket.id, y=mean, fill=Family))+
  geom_bar(stat="identity", position = "dodge") +
  geom_errorbar(aes(ymin=(mean-se), ymax=(mean+se)),
                width=.4, position=position_dodge(.9)) +
  theme(legend.title=element_blank()) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  facet_grid( . ~ factor(treatment,levels=c("larvae","control","low-ph","molybdate")), scales="free_x", space="free", labeller=ad_labeller) +
  xlab("Bucket ID") + ylab("Mean Percent Abundance") +
  scale_fill_manual(values = col_12) +
  coord_cartesian(ylim = c(1.6,32)) +
  scale_x_discrete(breaks=alpha.breaks,
                   labels=alpha.ids)

# avgs_3 families
# [1] Flavobacteriaceae      Rhodobacteraceae       Rubinisphaeraceae     
# [4] Saprospiraceae         Rhizobiaceae           A4b                   
# [7] Phycisphaeraceae       Crocinitomicaceae      Cryomorphaceae        
# [10] Cyanobiaceae           Pseudoalteromonadaceae Rubritaleaceae  

# acgs_5 families
# Flavobacteriaceae Rhodobacteraceae  Saprospiraceae    Crocinitomicaceae
# Rubritaleaceae    Phycisphaeraceae 

# change x axis labels??????????



#################################################################




# PICRUST2
#################################################################
# July 28 2020
# exporting data from previous analyses to give to David as inputs for PICRUST2

library(ShortRead)
library(biom)
setwd("/Users/roxannebanker/bioinformatics/projects/oyster-oa/microbiome/picrust2/")

# to make biom
picrust.seqtab <- otu_table(ps.all)

asv <- as(picrust.seqtab,"matrix")
asv_biom <- make_biom(data=asv)
write_biom(asv_biom, "oyster_seqs.biom")



# make fasta 
# make tax_table from ps.all 
picrust.taxtab <- tax_table(ps.all)

# pull out sequences
seqs_study <- getSequences(picrust.taxtab)

# export, then i subsequently change suffix to .fna
writeFasta(seqs_study, file = "oyster_seqs.fasta")




#################################################################





