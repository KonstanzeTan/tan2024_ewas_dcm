#' -----------------------------------------------------------------------------
#' cis-eQTM analysis (MatrixeQTL)
#' 
#' @description This script performs cis-eQTM (expression quantitative trait 
#' methylation) analysis to investigate relationships between sentinel CpG methylation
#' and proximal gene expression (<1Mb based on gene transcriptional start site).
#'
#' @author Konstanze Tan <konstanz001@e.ntu.edu.sg) 
#'
#' @date Friday Dec 13 2024
#' -----------------------------------------------------------------------------
# load betas
phe_magnet <- read.delim("phe_expr.txt")
beta_magnet  <- readRDS("beta_QN_detP001_marker095_rmch_magnet.rds")
beta_magnet <- as.data.frame(beta_magnet)
beta_magnet <- beta_magnet[,colnames(beta_magnet) %in% phe_magnet$BeadChip]

# load covariates
cov_magnet <-read.delim("cvrt_306_AgeGenderRaceRINPF7891213", check.names=F) #rownames: covariates (factor encoded), colnames: samples
cov_magnet <-as.data.frame(t(cov_magnet))

# read in gene expression file
gene_expr<-read.delim("magnet_expr_filt_norm", check.names=F) #rownames: ENSGID, colnames: BeadChip

# check if column order is standardized
gene_expr <- gene_expr[, colnames(beta_magnet)]
colnames(beta_magnet)==colnames(gene_expr)

cov_magnet <- cov_magnet[, colnames(beta_magnet)]
colnames(beta_magnet)==colnames(cov_magnet)

# read in CpG and Gene position information

cpgspos <- read.delim("cpgspos_147033_sentinels_backg") # 0-start to correspond with gene
cpgspos <- cpgspos[cpgspos$CpG %in% sentinels,]
genepos <- readRDS("/genepos_16465_auto_hg19.rds")#geneid, chr, left, right

# prepare SlicedData object for gene expression

gene = SlicedData$new()
gene$CreateFromMatrix(as.matrix(gene_expr))
gene$fileSliceSize = 2000# read file in pieces of 2,000 rows

# prepare SlicedData object for covariates
cvrt = SlicedData$new(); # leave as this if not running with covariates 
cvrt$CreateFromMatrix(as.matrix(cov_magnet))

# prepare SlicedData object for cpgs (matrix, rownames=CpGs, colnames=samples)
cpgs = SlicedData$new()
cpgs$CreateFromMatrix(as.matrix(beta_magnet))
cpgs$fileSliceSize = 2000 #read file in pieces of 2,000 rows (too big= faster to compute but slower to load; balance!)

# Run eQTM analysis

output_file_name = c("eQTM_AgeGenderRaceRINPF7891213_cis_1Mb_sentinel_magnet")
  
  me = Matrix_eQTL_main(
    snps = cpgs, #sliced data object 
    gene = gene,
    cvrt = cvrt,
    output_file_name.cis = output_file_name,
    pvOutputThreshold.cis = 1,
    pvOutputThreshold = 0, # set=0 if analysing cis-pairs
    snpspos=cpgspos,
    genepos =genepos,
    useModel = modelLINEAR, # standard additive linear model rather than dominant or recessive linear model (applies to genetic data)
    errorCovariance = numeric(), # independence for gene expression
    verbose = TRUE,
    cisDist=1e6, # distance for local gene-CpG pairs. Make sure its TSS-centered!
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE); # calculate FDR?
