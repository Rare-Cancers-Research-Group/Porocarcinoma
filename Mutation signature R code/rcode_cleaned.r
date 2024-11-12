# prefiltering done in bash due to signeRs' use of VariantAnnotation module doesn't do filtering
# within R
# variant tables are dragens .hard-filtered.vcf.gz
# for i in *.vcf.gz;do n=${i%.vcf.gz}".som.vcf";zcat $i | perl -lane 'm/^#/ && print && next; m/PASS/ || next; m/SOMATIC/ || next; length($F[3])==1 && length($F[4])==1 || next; print;' > $n;done

library(dplyr)
library(reshape2)
library(kableExtra)
library(ggplot2)
library(gridExtra)
library(BSgenome.Hsapiens.UCSC.hg38)
library(mutSignatures)

hg38 = BSgenome.Hsapiens.UCSC.hg38

vcf_files=c("P10.hard-filtered.som.vcf","P14.hard-filtered.som.vcf","P23.hard-filtered.som.vcf","P32.hard-filtered.som.vcf","P5.hard-filtered.som.vcf","P11.hard-filtered.som.vcf","P15.hard-filtered.som.vcf","P24.hard-filtered.som.vcf","P38.hard-filtered.som.vcf","P6.hard-filtered.som.vcf","P12.hard-filtered.som.vcf","P16.hard-filtered.som.vcf","P26.hard-filtered.som.vcf","P3.hard-filtered.som.vcf","P7.hard-filtered.som.vcf","P13.hard-filtered.som.vcf","P17.hard-filtered.som.vcf","P2.hard-filtered.som.vcf","P4.hard-filtered.som.vcf","P9.hard-filtered.som.vcf")
sample_names=c("P10","P14","P23","P32","P5","P11","P15","P24","P38","P6","P12","P16","P26","P3","P7","P13","P17","P2","P4","P9")
poro = importVCFfiles(vcfFiles=vcf_files, sampleNames=sample_names)
poroData = processVCFdata(poro, BSGenomeDb=hg38, sample_colName="SAMPLEID")
poro.counts = countMutTypes(mutTable = poroData, mutType_colName = "mutType", sample_colName = "SAMPLEID")

# number of signatures from signeR
num.sign = 5
poro.params = mutSignatures::setMutClusterParams(num_processesToExtract = num.sign, num_totIterations = 20, num_parallelCores = 8)
poro.analysis = decipherMutationalProcesses(input = poro.counts, params = poro.params)
poro.sig = poro.analysis$Results$signatures
poro.exp = poro.analysis$Results$exposures
# Plot signature 1 (standard barplot, you can pass extra args such as ylim)
msigPlot(poro.sig, signature = 2, ylim = c(0, 0.10))
msigPlot(poro.exp, main = "samples")
png(file="5_signatures.png",
width=1200, height=1500)
par(mfrow=c(3,2))
msigPlot(poro.sig, signature = 1, ylim = c(0, 0.10))
msigPlot(poro.sig, signature = 2, ylim = c(0, 0.10))
msigPlot(poro.sig, signature = 3, ylim = c(0, 0.10))
msigPlot(poro.sig, signature = 4, ylim = c(0, 0.10))
msigPlot(poro.sig, signature = 5, ylim = c(0, 0.10))
dev.off()

msigPlot(poro.exp, main = "samples")
xprt = coerceObj(x = poro.sig, to = "data.frame")
write.csv(xprt, "5_signatures.csv")
xprt = coerceObj(x = poro.exp, to = "data.frame")
write.csv(xprt, "5_signatures_exposures.csv")

# Compare de-novo signatures with selected COSMIC signatures

cosmic = mutSignatures::as.mutation.signatures(read.delim("COSMIC_v3.3.1_set1.txt", header = TRUE, row.names=1))
match_sig = matchSignatures(mutSign = poro.sig, reference = cosmic, threshold = 0.45, plot = TRUE)
match_sig$plot + ggtitle("Match to COSMIC signatures.")


# signeR part
library(VariantAnnotation)
library(rtracklayer)
library(signeR)
library(BSgenome.Hsapiens.UCSC.hg38)

hg38 = BSgenome.Hsapiens.UCSC.hg38
vcf_files=c("P10.hard-filtered.som.vcf","P14.hard-filtered.som.vcf","P23.hard-filtered.som.vcf","P32.hard-filtered.som.vcf","P5.hard-filtered.som.vcf","P11.hard-filtered.som.vcf","P15.hard-filtered.som.vcf","P24.hard-filtered.som.vcf","P38.hard-filtered.som.vcf","P6.hard-filtered.som.vcf","P12.hard-filtered.som.vcf","P16.hard-filtered.som.vcf","P26.hard-filtered.som.vcf","P3.hard-filtered.som.vcf","P7.hard-filtered.som.vcf","P13.hard-filtered.som.vcf","P17.hard-filtered.som.vcf","P2.hard-filtered.som.vcf","P4.hard-filtered.som.vcf","P9.hard-filtered.som.vcf")
mut = matrix(ncol=96,nrow=0)
for(i in vcf_files) {
    vo = readVcf(i, "hg38")
    # sample name (should pick up from the vcf automatically if available)
    # colnames(vo) = i
    m0 = genCountMatrixFromVcf(hg38, vo)
    mut = rbind(mut, m0)
    print(paste(i, "done,", date()))
}
# dim(mut) # show mutation matrix dimensions
# [1] 20 96

target_regions = import(con="Twist_Core_Exome_RefSeq_targets_hg38.bed", format="bed")
opp = genOpportunityFromGenome(hg38, target_regions, nsamples=nrow(mut))

signatures = signeR(M=mut, Opport=opp)
# The optimal number of signatures is 5.
# Running  Gibbs sampler for 5 signatures...Done.
