### load packages
library(BiocManager)
library(AnnotationHub)
library(rtracklayer)
library(Gviz)
library(Rsamtools)
library(biomaRt)
library(GenomicRanges)
library(GenomicFeatures)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(IRanges)
library(tidyverse)

### set working directory
setwd("D:/R-scripts/ChIPseq")

### import bigWig file (downloaded from ENCODE project db)
ChIP_bigWig <- import.bw("D:/01_Private Dateien/Biochemie Studium/Masterstudium/01 - Molecular Biosciences - Major Cancer Biology - RKU Heidelberg/4th semester - master thesis/ENCODE_ChIPseq/ENCFF796GHK.bigWig", as="GRanges")
ChIP_bigWig2 <- import.bw("D:/01_Private Dateien/Biochemie Studium/Masterstudium/01 - Molecular Biosciences - Major Cancer Biology - RKU Heidelberg/4th semester - master thesis/ENCODE_ChIPseq/ENCFF677COF.bigWig", as="GRanges")

### import narrow peak BED file
extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer")

# MYC
ChIP_bed <- import("D:/01_Private Dateien/Biochemie Studium/Masterstudium/01 - Molecular Biosciences - Major Cancer Biology - RKU Heidelberg/4th semester - master thesis/ENCODE_ChIPseq/ENCFF838NIW.bed", format = "BED", extraCols = extraCols_narrowPeak)

# MAX
ChIP_bed2 <- import("D:/01_Private Dateien/Biochemie Studium/Masterstudium/01 - Molecular Biosciences - Major Cancer Biology - RKU Heidelberg/4th semester - master thesis/ENCODE_ChIPseq/ENCFF481AOS.bed", format = "BED", extraCols = extraCols_narrowPeak)

## set genome to "hg19" if not already set to this particular reference genome
genome(ChIP_bigWig) = "hg38" 
gen <- genome(ChIP_bigWig)
gen <- ChIP_bigWig@seqinfo@genome
gen

### subset bigWig file to reduce size for subsequent analysis select your chromosome of interest etc.
ChIPseq_chr <- keepSeqlevels(ChIP_bigWig, paste0("chr","16"), pruning.mode="coarse")
ChIPseq2_chr <- keepSeqlevels(ChIP_bigWig2, paste0("chr","16"), pruning.mode="coarse")

### save bigWig subset of your chromosome of interest 
#export(K562_KLF1_chr8, "K562_KLF1_chr1.bigWig")

### create GenomeAxisTrack object, which acts as coordinate axis for the analysed genomic regions
GT = GenomeAxisTrack(exponent=3) # gives distance in kilobases kb

## select genomic ranges of chr (hg38.p13 chr lenght)
chr1 <- 248956422
chr2 <- 242193529
chr3 <- 198295559
chr4 <- 190214555
chr5 <- 181538259
chr6 <- 170805979
chr7 <- 159345973
chr8 <- 145138636
chr9 <- 138394717
chr10 <- 133797422
chr11 <- 135086622
chr12 <- 133275309
chr13 <- 114364328
chr14 <- 107043718
chr15 <- 101991189
chr16 <- 90338345
chr17 <- 83257441
chr18 <- 80373285
chr20 <- 58617616
chr19 <- 64444167
chr22 <- 46709983
chr21 <- 50818468
chrX <- 156040895
chrY <- 57227415

### create gene-region track from hg38.p13 reference genome biomart
hg38 <- useMart("ensembl", "hsapiens_gene_ensembl")
gTrack <- BiomartGeneRegionTrack(biomart = hg38, chromosome = "chr16", start = 1, end = chr12, name = "Transcripts", shape="arrow")

## chromosome name
chr <- as.character(unique(seqnames(ChIPseq_chr)))
chr

### Ideogram track for visualization of Chromosomes etc.
IT = IdeogramTrack(genome = "hg38", chromosome = chr)
IT

### Datatrack ok K562 chr-subset
DT1 <- DataTrack(range = ChIPseq_chr, genome = "hg38", type = c("polygon"), name = "Max\n", lwd.mountain=2,
                col.mountain = "dodgerblue3", fill.mountain = c("black","dodgerblue3"), ylim=c(0,160))

DT2 <- DataTrack(range = ChIPseq2_chr, genome = "hg38", type = c("polygon"), name = "c-Myc\n", lwd.mountain=2,
                col.mountain = "firebrick3", fill.mountain = c("black","firebrick3"), ylim = c(0,160))

### Annotationtrack of ChIPseq peaks
AT1 = AnnotationTrack(range= ChIP_bed, genome = "hg38", chromosome = "chr16", name = "peaks", featureAnnotation = "id",
                      fill = "firebrick3", cex.feature = 0, strands = c("+","-","*"), shape = "fixedArrow")
AT2 = AnnotationTrack(range= ChIP_bed2, genome = "hg38", chromosome = "chr16", name = "peaks", featureAnnotation = "id",
                      fill = "dodgerblue3", cex.feature = 0, strands = c("+","-","*"), shape = "fixedArrow")

# plot Annotation tracks alone without signal curves from bigWig
track.list=list(IT, GT, gTrack, AT1, AT2)

plotTracks(track.list,from = 4257186, to = 4273075, transcriptAnnotation="symbol", window="auto",
           cex.title=1, cex.bands = 0.4, showBandId=TRUE, fontsize=13, fontface.main = 1,
           sizes = c(0.2,0.2,0.5,0.15,0.15))

### plot tracks of ChIPseq data (select your relevant datatracks)
track.list=list(IT, GT, gTrack, AT1, DT2, AT2, DT1)

plotTracks(track.list,from = 4257000, to = 4274000, transcriptAnnotation="symbol", window="auto",
           cex.title=1, cex.bands = 0.4, showBandId=TRUE, fontsize=13, fontface.main = 1,
           sizes = c(0.2,0.2,0.5,0.2,0.4,0.2,0.4))



### plot your favourite gene locus without tracks from bidWig files
# full hg38.p13 gene region track and ideogram
gTrack <- BiomartGeneRegionTrack(biomart = hg38, name = "Transcripts", shape="arrow", chromosome = "chr11")
IT = IdeogramTrack(genome = "hg38", chromosome = "chr11")

track.list=list(IT, GT, gTrack)

plotTracks(track.list,from = 5210000, to = 5520000, transcriptAnnotation="symbol", window="auto",
           cex.title=1, cex.bands = 0.4, showBandId=TRUE, fontsize=14, fontface.main = 1, sizes = c(0.10,0.1,0.8))
