# Week 6 Linkage Analysis
# Sungmin Lee

# Objectives
# Describe the concept of linkage analysis
# Demonstrate how to install the GenABEL package from BioConductor
# Highlight the different ways to analyze data and to generate various plots like the Manhattan Plot
# Discuss the Bonferroni correction and demonstrates how to can use genotype and phenotype data to look at Expression Quantitative Traits Loci(eQTLs)

# Linkage Analysis
# - one of the most powerful gene-discovey tools to be developed
# - Linkage analysis uses SNPs (markers) and traits (phenotype) and does a massive correlation analysis across a population
#   (trait-phenotype association)
# - This is data and time intensive
# - two types of traits can be examined
# 1. continuous trait (height and such)
# 2. binary traits (dominant and recessive) which can be muddled by penetrance

########
# GWAS #
########
# Early
# - Genome-wide association studies (GWAS) focused on binary conditions and SNP profiles (SNP for possible heterozygous states)
# Simple analysis with GenABEL package

library(GenABEL)
data <- load.gwaa.data(phenofile="phenotype.dat",genofile="genotype.raw")

# access the genotype at marker 3 through 5
gtdata(data)[,3:5]
# or if you want the name of marker (SNP) 17000
snpnames(gtdata(data)[,17000])
# or retrieve the actual genotypes at specific marker for specific individuals
as.character(gtdata(data)[1:4,1177])
# to see the reference ("wild type") allele
refallele(gtdata(data)[,1777])

# Hardy-Weinberg equilibrium states that the population is in equilibrium for the distribution of alleles
summary.snp.data(gtdata(data))
# p.11,p.12,p.22 give actual number for genotype, Fmax column shows the deviation from HWE

# distribution of continuous trait
ct <- phdata(data)$ct
hist(ct,col="slateblue",main="Distribution of ct")

# check if contribution trait depends on other vars
boxplot(ct~phdata(data)$sex, col=c("blue","red"),xlab="sex",names=c("F","M"),main="ct by sex")
# what does ct~ do? -> depends on, check whether continuous trait depends on sex column

# Quality control analysis
?check.marker

qc <- check.marker(data, call = 0.95, perid.call=0.95, maf=1e-08, p.lev=1e-08)
# call, marker variable set to
# perid, individual called marker 
# maf, minor allele frequency
# p.lev, pvalue for HWE

# Create new dataset containing individuals and snps that are ok
data.qc <- data[qc$idok, qc$snpok]

# can test each marker one at a time using Cochran-Armitage Trend Test
# if we want to see if ct is associated with any of the markers in our dataset.
an <- qtscore(ct~1, data=data.qc, trait="gaussian")
plot(an,col=c("olivedrab","slateblue"),cex=.5, main="Manhattan plot")
# maybe strong correlation with ct with chromosome 2 - big spike

# before we need to see if p-value for each of many tests are significant
# looking at "lambda" or genomic inflation factor 
# can result in lower than expected p values due to cryptic relatedness or population structure
estlambda(an[,"P1df"],plot=T)
# estimate = 1.182951 above 1 -> small amount of genomic inflation

# Pc1df contains genomic control information
estlambda(an[,"Pc1df"], plot=T)
# running lambda -> much closer to 1
# runnning manhattan plot with Pc1df
plot(an, col=c("olivedrab","slateblue"), cex=.5, main="Manhattan plot", df="Pc1df")

# there is a possibility that p-value is significant by chance
# normally we can increase p-value to compensate
# but Bonferroni correlation does divide each p-value for each test by the number of tests
pval.threshold <- 0.05
bonferroni <- -log10(pval.threshold/nids(data.qc))

# let's add this to manhattan plot
plot(an, col=c("olivedrab","slateblue"), cex=.5, main="Manhattan plot", df="Pc1df")
abline(h=bonferroni, lty=3, color="red")


########
# eQTL #
########
# In addition, genotype and phenotype data to look at Expression Quantitative Traits Loci, or eQTLs
# eQTLs are regions of genome containing DNA sequence variants that influence expression level of one or more genes. 
library(MatrixEQTL)

useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# set significant threshold
pvOutputThreshold = 1e-2 # only SNPs that have an effect on gene expression with a significance of 10-2 or greater will be retained

# covariance matrix
errorCovariance = numeric()

# Load data
snps = SlicedData$new() # initialize a new and blank dataframe called snps
snps$fileDelimiter = "\t" # tell MatrixEQTL that values are separated by tabs
snps$fileOmitCharacters = "NA" # tells MatrixEQTL to ignore values that are NA
snps$fileSkipRows = 1 # MatrixEQTL that the first row is a label row
snps$gileSkipsColumns = 1 # tells MatrixEQTL that the first column is a label.
snps$fileSliceSize = 2000 # MatrixEQTL to read chunks of 2000 rows
snps$LoadFile("SNP.txt") # MatrixEQTL what our filename is

gene = SlicedData$new()
gene$fileDelimiter = "\t"      
gene$fileOmitCharacters = "NA" 
gene$fileSkipRows = 1          
gene$fileSkipColumns = 1       
gene$fileSliceSize = 2000      
gene$LoadFile("GE.txt")

cvrt = SlicedData$new()
cvrt$fileDelimiter = "\t"      
cvrt$fileOmitCharacters = "NA" 
cvrt$fileSkipRows = 1          
cvrt$fileSkipColumns = 1       
cvrt$LoadFile("Covariates.txt")

# Run analysis
me = Matrix_eQTL_engine(snps=snps,gene=gene,
                        cvrt=cvrt, output_file_name = "outfile.txt",
                        pvOutputThreshold = pvOutputThreshold, useModel = useModel,
                        errorCovariance = errorCovariance, verbose=TRUE,
                        pvalue.hist=FALSE, min.pv.by.genesnp=FALSE, noFDRsaveMemory=FALSE)

# output file is a list of SNP-gene combinations that exceeded our threshold
# indicating those SNPs somehow affecting the expression levels of the indicated genes. 
# For example, SNP_11 changes the expression level of Gene_06

# Results:
cat('Analysis done in: ', me$time.in.sec, 'seconds', '\n')
cat('Detected eQTLs', '\n')
show(me$all$eqtls)

# see if there is a difference between SNPs that are close to our gene and those are distant
# temp files to store cis and trans file
output_file_name_cis = tempfile()
output_file_name_tra = tempfile()

pvOutputThreshold_cis = 1e-2
pvOutputThreshold_tra = 1e-2

# establish maximum cis distance
# tell MatrixEQTL that gene/SNP pairs that are further part than this should essentially be considered trans
cisDist = 1e6 # 1 million base pairs

snpspos = read.table("snpsloc.txt",header=TRUE, stringsAsFactors = FALSE)
genepos = read.table("geneloc.txt",header=TRUE, stringsAsFactors = FALSE)

# Run analysis
me = Matrix_eQTL_main(snps=snps, gene = gene, cvrt = cvrt,
                      output_file_name = output_file_name_tra, pvOutputThreshold = pvOutputThreshold_tra, 
                      useModel = useModel, errorCovariance = errorCovariance, verbose = TRUE, 
                      output_file_name.cis = output_file_name_cis, pvOutputThreshold.cis = pvOutputThreshold_cis,
                      snpspos = snpspos,genepos=genepos, cisDist = cisDist, pvalue.hist = "qqplot", 
                      min.pv.by.genesnp = FALSE, noFDRsaveMemory = FALSE)
# min.pv.by.genesnp tells MatrixEQTL to save all p values even if they do not meet the threshold if set for TRUE.
# noFDRsaveMemory tells MatrixEQTL to write the output directly to file when set for TRUE – this saves memory.

# Result
show(me$cis$eqtls)
cat('Detected distant eQTLs: ', '\n')
show(me$trans$eqtls)

# Q-Q plot or quantile-quantile plot
# a graphical technique for determining if two data sets come from populations with a common distribution
# plot of quantiles of the first data against the quantiles of the second data set. 
# If not, symmetric distribution, Q-Q plot suggests that they don't both come from same population
# skewed towards a local distribution, meaning that SNP/Gene pairs that are closer tend to show more association

###########
# Summary #
###########
# The concept of linkage analysis
# The installation of GenABEL package from BioConductor
# The different ways to analyze data and to generate various plots like the Manhattan Plot
# The concept of Bonferroni correction and how to can use genotype and phenotype data to look at Expression Quantitative Traits Loci(eQTLs)