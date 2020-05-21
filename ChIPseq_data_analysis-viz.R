#### Analysis and visualization of ChIPseq data (retrieved from the ENCODE project database)

### load necessary packages
library(BiocManager)
library(tidyverse)
library(Rsamtools)
library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
library(ggraph)
library(clusterProfiler)
library(ChIPpeakAnno)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(reshape2)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TFBSTools)
library(Biostrings)
library(ChIPseeker)
library(biomaRt)
library(phylotools)
library(msigdbr)
library(motifRG)
library(motifStack)
library(JASPAR2018)
library(AnnotationHub)
library(Gviz)

### set wd
setwd("path to your working directory")

# load txdb gene model and shorten name for convenience
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

### load files (narrowpeak = BED format)
extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric", qValue = "numeric", peak = "integer")

# Snyder Stanford convervative IDR peaks
ChIP_bed <- import("~/ENCFF875JHB.bed", format = "BED", extraCols = extraCols_narrowPeak)

### annotate the peaks with precompiled ensembl annotation
data(TSS.human.GRCh38)
ucsc.hg38.knownGene <- genes(TxDb.Hsapiens.UCSC.hg38.knownGene)

# use biomaRt for annotation with ensmbl
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# annotate GRanges
ChIP.anno <- annotatePeakInBatch(ChIP_bed, AnnotationData = ucsc.hg38.knownGene)
ChIP.anno <- addGeneIDs(annotatedPeak = ChIP.anno, orgAnn = "org.Hs.eg.db", feature_id_type = "entrez_id", IDs2Add = c("ensembl","symbol"), mart = ensembl)

## calculate percentage of peaks in promoters etc.
ChIP_aCR <- assignChromosomeRegion(ChIP.anno, nucleotideLevel=FALSE, precedence=c("Promoters", "immediateDownstream", 
                                   "fiveUTRs", "threeUTRs", "Exons", "Introns"), TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene)

# create dataframes
ChIP_aCR_p <- ChIP_aCR$percentage
ChIP_aCR_p <- as.data.frame(ChIP_aCR_p)

# plot bargraphs
ggplot(ChIP_aCR_p, aes(x=`subjectHits`, y=`Freq`)) + geom_bar(position= "dodge", colour="black", width=.8, stat = "identity",
  show.legend = TRUE) + ylab("percentage [%]\n") + xlab("\ngenomic elements") + theme_bw() + labs(title = "ChIPseq peaks") +
  theme(plot.title = element_text(color="black", size=16, face= "bold"),
        axis.title.x = element_text(color="black", size=16, face= "bold"),
        axis.text.x = element_text(angle = -90, vjust = 0.5, hjust = 0.0, color="black", size=12, face= "bold"),
        axis.text.y = element_text(color="black", size=12, face = "bold"),
        axis.title.y = element_text(color="black", size=16, face="bold"))

# peak heatmap
peakHeatmap(ChIP.anno, TxDb = txdb, upstream = 2000, downstream = 2000, color = "firebrick3", title = "peak heatmap")

# density plot
plotAvgProf2(ChIP.anno, TxDb = txdb, upstream = 2000, downstream = 2000, conf = 0.95, resample = 1000,
             xlab = "Genomic Region (5'->3')", ylab = "Read Count Frequency ChIPseq peaks")

### get promoters via from reference genome hg38
genes <- genes(txdb)
promoters <- promoters(genes, upstream = 2000, downstream = 200) #  define promoter range around TSS
promoters

# subset annotaded ChIPseq peaks within promoters
ChIP_prom <- subsetByOverlaps(ChIP.anno, promoters)
ChIP_prom

### get dataframe of GRanges object of ChIPseq data
ChIP_prom <- getAllPeakSequence(ChIP_prom, genome = Hsapiens)

ChIP_prom_df <- as.data.frame(ChIP_prom)
ChIP_prom_df <- tibble::rownames_to_column(ChIP_prom_df, "Ranges")

# create fasta
ChIP_seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, ChIP_prom)

# export fasta
writeXStringSet(ChIP_seq, "./ChIP_seq.fasta")

# import fasta sequences
ChIP_seq <- readDNAStringSet("ChIP_seq.fasta")

#count G and C nucleotides in sequence and divide by all nucleotides etc.
ChIP_prom_df_G <- str_count(ChIP_prom_df$sequence, "G")
ChIP_prom_df_C <- str_count(ChIP_prom_df$sequence, "C")
ChIP_prom_df_GCAT <- str_count(ChIP_prom_df$sequence, "")
ChIP_prom_df_GC <- ((ChIP_prom_df_C + ChIP_prom_df_G) / ChIP_prom_df_GCAT) * 100                                

# add % GC-content as column to dataframe
ChIP_prom_df["GC_content"] <- ChIP_prom_df_GC
colnames(ChIP_prom_df)

# save as txt
write_tsv(ChIP_prom_df, "ChIPseq_prom_df.txt")

### GO enrichment analysis with Clusterprofiler
ChIP_prom_df <- read_tsv(file = "~/ChIPseq_prom_df.txt") # import table if needed

## perform enrichment analysis
# use known genes from hg38 reference genome as background genes -> annotate via org.Hs.eg.db to get gene symbols
backgroundgenes <- ucsc.hg38.knownGene$gene_id
backgroundgenes.anno <- addGeneIDs(annotatedPeak=backgroundgenes, 
                                   orgAnn="org.Hs.eg.db", feature_id_type="entrez_id", c("ensembl","symbol"), mart = ensembl)

# create term2gene dataframe , category = "H" == hallmark gene sets , "C2" == curated gene sets , "C3" == motif genesets
m_t2g <- msigdbr(species = "Homo sapiens") %>% dplyr::select(gs_name, gene_symbol)
head(m_t2g)

em <- enricher(ChIP_prom_df$symbol, TERM2GENE=m_t2g, pvalueCutoff = 0.05, pAdjustMethod = "BH",
               universe=backgroundgenes.anno$symbol, minGSSize = 10, maxGSSize = 500, qvalueCutoff = 0.05)
enricher_result <- em@result

# save txt file
write_tsv(enricher_result, "ChIPseq_enricher_result.txt")

### plot GO enrich results
enricher_result = enricher_result[order(enricher_result$p.adjust), c(1:9)]
enricher_result <- enricher_result[c(1:8),]
nLog10_adjP <- -log10(enricher_result$p.adjust)
enricher_result["log10(p.adjust)"] <- nLog10_adjP
ordeR <- c(8:1)
enricher_result["order"] <- ordeR

ggplot(enricher_result, aes(x=reorder(`Description`, `order`), y=`log10(p.adjust)`)) +
  geom_bar(position= "dodge", fill="dodgerblue3", width=.8, stat = "identity", show.legend = TRUE, colour="black") + 
  theme_bw() + labs(title = "") + ylab("-log10(adj. p-value)\n") + xlab("") +
  scale_y_continuous(expand = c(0,0), limits = c()) + 
  geom_text(aes(y=`p.adjust`, label=round(`log10(p.adjust)`, digits = 2)), vjust=0.5, hjust=-0.5, size=5.5,
            position=position_stack(vjust = 0.5), colour="white") + 
  theme(plot.title = element_text(color="firebrick3", size=18, face= "bold", hjust=1.0),
        axis.title.x = element_text(color="black", size=18, face= "bold", hjust=0.5),
        axis.text.x = element_text(angle = 0, vjust = 1, hjust = 0, color="black", size=18),
        axis.text.y = element_text(color="black", size=14, face = "bold"),
        axis.title.y = element_text(color="black", size=18, face="bold")) + coord_flip()

### motif analysis
ChIP_prom_f <- ChIP_prom[order(ChIP_prom$signalValue, decreasing = TRUE)]

top500_peaks <- head(ChIP_prom_f, n=500)

# reduce peaks to merge nearby located peaks
#reduced_peaks <- reduce(top500_peaks)

# resize peak sequences (get only +/- 100 bp from center)
resized_peaks <- resize(top500_peaks, width=200, fix="center")

# get peak sequences
peak_seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, resized_peaks)

set.seed(1) # shuffle sequences of peak_seq and use as background for analysis
control_seq <- DNAStringSet(lapply(peak_seq, sample))

### motif analysis
set.seed(1)

#calculate motif enrichments
motifs <- findMotifFgBg(fg.seq = peak_seq, bg.seq = control_seq, max.motif = 1, enriched.only = TRUE, max.width = 15, both.strand = TRUE)

# refine motif edges
refined_motifs <- lapply(motifs$motifs, function(x){motifRG::refinePWMMotifExtend(motifs = x@match$pattern, seqs = peak_seq, max.width = 15)})

# plot seqlogo
seqLogo(refined_motifs[[1]]$model$prob)
plotMotifLogo(refined_motifs[[1]]$model$prob, motifName = "", font = "Helvetica-Bold", xlcex = 2.9, ylcex = 2.9, ncex = , xlab = "\nPosition", ylab = "Bits\n")

### motif annotation
unknown_motif <- refined_motifs[[1]]$model$prob

# create PWM (position weight matrix)
unknown_pwm <- PWMatrix(ID = "unk", profileMatrix = unknown_motif)

# get PWM library from JASPAR core motif database
pwm_library <- getMatrixSet(JASPAR2018, opts = list(collection = "CORE", species = "Homo sapiens", matrixtype = "PWM"))

# find motifs that are simmilar to the motif (pearson correlation)
pwm_sim <- PWMSimilarity(pwm_library, unknown_pwm, method = "Pearson")

pwm_library_list = lapply(pwm_library, function(x){data.frame(ID = ID(x), name = name(x))})
pwm_library_dt = dplyr::bind_rows(pwm_library_list)
pwm_library_dt$similarity = pwm_sim[pwm_library_dt$ID]
pwm_library_dt = pwm_library_dt[order(-pwm_library_dt$similarity),]
pwm_library_dt$order = order(pwm_library_dt$similarity)
head(pwm_library_dt)

ggplot(pwm_library_dt[1:8,], aes(x=reorder(`name`, `order`), y=`similarity`)) + 
  geom_bar(position= "dodge", fill="navy", width=.8, stat = "identity", show.legend = TRUE, colour="black") + 
  theme_bw() + labs(title = "Motif alignment") + ylab("\nPWM similarity (pearson)") + xlab("") +
  scale_y_continuous(expand = c(0,0), limits = c(0,0.9)) + 
  geom_text(aes(y=`similarity`, label=round(`similarity`, digits = 2)), vjust=0.5, hjust=-0.5, size=5.5, position=position_stack(vjust = 0.5), colour="white") + 
  theme(plot.title = element_text(color="black", size=20, face= "bold", hjust=0.5),
        axis.title.x = element_text(color="black", size=18, face= "bold", hjust=0.5),
        axis.text.x = element_text(angle = 0, vjust = 1, hjust = 0, color="black", size=18),
        axis.text.y = element_text(color="black", size=14, face = "bold"),
        axis.title.y = element_text(color="black", size=18, face="bold")) + coord_flip()


### plot a specific gene locus using Gviz
# create genome axis track
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

# BiomartGeneRegionTrack
hg38 <- useMart("ensembl", "hsapiens_gene_ensembl") # enter chr for gene region track
gTrack <- BiomartGeneRegionTrack(biomart = hg38, chromosome = "chr12", start = 1, end = chr12, name = "Transcripts", shape="arrow")

# Ideogramtrack
IT = IdeogramTrack(genome = "hg38", chromosome = "chr12")

# modify ChIP_prom to have arrows in annotation track (feature stand + - * vs strand *)
colnames(ChIP_prom_df)[colnames(ChIP_prom_df) == "feature_strand"] <- "strand"
ChIP_prom_df <- ChIP_prom_df[,c(1:5,7:25)] # remove original strand info

ChIP_bed <- makeGRangesFromDataFrame(ChIP_prom_df, keep.extra.columns=TRUE, ignore.strand=FALSE, seqinfo=NULL,
                                      seqnames.field=c("seqnames","seqname","chromosome","chrom","chr", "chromosome_name","seqid"),
                                      start.field="start", end.field=c("end", "stop"), strand.field="strand", starts.in.df.are.0based=FALSE)

# Annotationtrack for ChIPseq peaks
AT = AnnotationTrack(range= ChIP_bed, genome = "hg38", chromosome = "chr12", name = "peaks", featureAnnotation = "id",
                      fill = "firebrick3", cex.feature = 0, feature_strand = c("+","-","*"), shape = "fixedArrow", feature = ChIP_bed@strand)

### import bigWig file (downloaded from ENCODE project db) -> eg. GATA1 in K562 cell line
ChIP_bigWig <- import.bw("~/ENCFF932GTA.bigWig", as="GRanges")

# subset bigWIg for relevant chromosome
ChIPseq_chr <- keepSeqlevels(ChIP_bigWig2, paste0("chr","12"), pruning.mode="coarse")

DT <- DataTrack(range = ChIPseq_chr, genome = "hg38", type = c("polygon"), name = "GATA1", lwd.mountain=1,
                 col.mountain = "darkgrey", fill.mountain = c("black","dodgerblue3"), ylim=c())

# plot Annotation tracks alone without signal curves from bigWig
track.list=list(IT, GT, gTrack, AT, DT)

plotTracks(track.list,from = 119665000, to = 119682000, transcriptAnnotation = "symbol", window = "auto",
           cex.title = 1, cex.bands = 0.4, showBandId = TRUE, fontsize = 13, fontface.main = 1,
           sizes = c(0.15,0.2,0.65,0.25,0.4))
