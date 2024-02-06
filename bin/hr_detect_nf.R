#!/usr/local/bin/Rscript
library(VariantAnnotation)
library(signature.tools.lib)
args = commandArgs(trailingOnly=TRUE)

sample <- args[1]
snvexp_path <- args[2]
svexp_path <- args[3]
cnv_path <- args[4]
indel_path <- args[5]

genomev = 'hg38'
##make the empty input matrix
col_hrdetect <- c("del.mh.prop", "SNV3", "SV3", "SV5", "hrd", "SNV8")
input_matrix <- matrix(NA,nrow = 1,ncol = length(col_hrdetect),dimnames = list(sample,col_hrdetect))
####################################################################
####indels##########################################################
####################################################################
indeltab <- c(indel_path)
names(indeltab) <- sample
####################################################################
####snvs############################################################
####################################################################

#snv catalogues
snv_full <- read.csv(snvexp_path)
snv <- subset(snv_full, sample==sample)
rownames(snv) <- 1:length(rownames(snv))
input_matrix[sample, 'SNV3'] <- sv[1, 'SNV3']
input_matrix[sample, 'SNV8'] <- sv[1, 'SNV8']

####################################################################
####svs############################################################
####################################################################
sv_full <- read.csv(svexp_path)
sv <- subset(sv_full, sample==sample)
rownames(sv) <- 1:length(rownames(sv))
input_matrix[sample, 'SV3'] <- sv[1, 'SV3']
input_matrix[sample, 'SV5'] <- sv[1, 'SV5']


####################################################################
####cnvs############################################################
####################################################################
###make cnv input through modification of [sample]_CNVs.tsv MTR input
cnvs <- c(cnv_path)
names(cnvs) <- sample

####################################################################
####run HRDetect####################################################
####################################################################
res <- HRDetect_pipeline(input_matrix,
                         genome.v = genomev,
                         Indels_vcf_files = indeltab,
                         CNV_tab_files = cnvs)

df <- as.data.frame(res$hrdetect_output)
df['sample'] = sample
write.table(df, paste0(sample, '_hr_detect.tsv'),sep='\t', row.names = F,quote=F)

