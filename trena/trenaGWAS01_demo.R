# demo.R
# Assess GWAS SNPs In the Context of a TRENA Gene Model, based on
# data from a 2014 Nature article:
#
#  Biological insights from 108 schizophrenia-associated genetic loci
#  Schizophrenia Working Group of the Psychiatric Genomics Consortium
#  http://www.nature.com/nature/journal/v511/n7510/full/nature13595.html
#
# In this small demo we focus on just two SNPs in locus 5 of those 108 loci
# All the data needed for the analysis is contained within the R package
# used here, "trenaGWAS01"
#

library(trenaGWAS01)
stopifnot(packageVersion("trenaGWAS01") >= "0.99.6")

stopifnot(packageVersion("MotifDb") >= "1.19.17")
stopifnot(packageVersion("BrowserViz") >= "1.9.15")
stopifnot(packageVersion("httpuv") >= "1.3.5")
stopifnot(packageVersion("trena") >= "0.99.187")
stopifnot(packageVersion("trenaViz") >= "0.99.21")

# chunk 1: initialize the package, establish connection to trenaViz in   your web browser
# "tgwas" is created by the trenaGWAS01 package.  it provides relatively high-level
# function calls to trena, trenaViz, and MotifDb packages.  It also
# includes pre-calculated gene regulatory models for TSNARE1, along with variant
# data for locus 5.  The hg38 genome is used.

tgwas <- trenaGWAS01(gwasLocusNumber=5, targetGene="TSNARE1", targetGene.tss=142354831, quiet=TRUE)


# chunk 2: direct the igv.js genome browser to display a small (4kb)  portion of the GWAS locus 5
# igv's genome browser can be controlled by R function calls.  and most immediately
# relevant here: the region on display, set programmatically (see below) or intereactively
# with your mouse, defines the genomic region used in subsequent footprint and lookup and
# calculatesions

setRegion(tgwas, "chr8:142,230,914-142,234,913") # rs49776977 is at 142,232,793

# chunk 3:  load DNase footprints for the currently displayed region,
# from four databases hosted by BDDS (UChicago).   These four are the
# result processing the same ENCODE brain (mixed tissue) footprint
# reads, using two different methods (HINT and Wellington) and two
# different alignment seeds sizes.  This step currently takes a few
# minutes, longer that we like.  Future refinement to this package
# will provide choice of database, choice of alternate open chromatin
# data (ATAC-seq, for instance), and much improved speeds
# this takes about 30 seconds: 4 footprint databases, 4000 bases

tbl.fp <- loadFootprints(tgwas)

# chunk 4: load SNPs provided by the GWAS study.  the observed "index SNP" is in one track
# 10 additional imputed snps ("credible SNPS") in purple are in their own track

tbl.snps <- loadSNPs(tgwas)

# chunk 5: add one more track: imputed snps which fall in the previously
# loaded footprints

tbl.snpsInFp <- findSNPsInFootprints(tgwas, tbl.snps, tbl.fp)


# chunk 6: load a previously calculated trena model,  built upon
# RNA-seq data from the cerebellum.  (In another document we give
# several examples of gene regulatory model construction.)  The
# "assess..." function reports any transcription factor in the model
# which is gained or lost as a consequence of each SNP.
#
# each snp (currently 2) is assessed in these ways:
#
# 1) the wild-type sequence at snp location +/- "shoulder" is matched against about 500
#    human transcription factor binding motifs
#
# 2) Those motifs which match the wild-type sequence with a score >= 90% are then
#    associated (conservatively, using MotifDb) with their cognate transcription factors
#
# 3) We then ask: do any of these >= 90% match TFs appear in the previously calculated
#    trena model?
#
# 4) for those TFs which DO fall in the model, their rank in that model is reported
#
load(system.file(package="trenaGWAS01", "extdata", "tbl.geneModel.cer.RData"))
targetGene <- "TSNARE1"
assessSnpsInContextOfGeneModel(tgwas, tbl.snpsInFp, tbl.geneModel, targetGene, matchThreshold=90, shoulder=8)

# chunk 7: we also have a more detailed method for assessing SNPs.  This reports all we know about
# motifs & TFs observed, gained and lost in the area of interest, in the wild-type sequence and
# after injecting the SNP
# first collect a broad range of plausible motifs

jaspar2016.pfms <- query(MotifDb, "jaspar2016")
human.pfms <- query(jaspar2016.pfms, "sapiens")
mouse.pfms <- query(jaspar2016.pfms, "mus")
rat.pfms   <- query(jaspar2016.pfms, "rnorvegicus")
pfms <- as.list(c(human.pfms, mouse.pfms, rat.pfms))   # 619

# request the detailed assesment.  set the shoulder (as above), to 8 bases up and downstream.
# relax the match score to 80 so that you can see how MA0671.1 (NFIX) scores high in wt,
#  > 10% worse in the variant sequence

tbl.snpInfo <- assessSnp(trena, pfms, "rs4976977", 8, 80)

# one of many possible manipulations of this data.frame (transposed for easier reading)
# notice how the motifRelativeScore is 0.9675 in wt, 0.8619 in mutant, and the
# conserved second base of the motif (reverse complement) is T -> A

as.data.frame(t(tbl.snpInfo[grep("NFIX", tbl.snpInfo$motifName),]))

# motifName                      Hsapiens-jaspar2016-NFIX-MA0671.1             Hsapiens-jaspar2016-NFIX-MA0671.1
# status                                                        wt                                           mut
# assessed                                                 in.both                                       in.both
# motifRelativeScore                                     0.9675002                                     0.8618954
# delta                                                          0                                             0
# signature          Hsapiens-jaspar2016-NFIX-MA0671.1;142232792;- Hsapiens-jaspar2016-NFIX-MA0671.1;142232792;-
# chrom                                                       chr8                                          chr8
# motifStart                                             142232792                                     142232792
# motifEnd                                               142232800                                     142232800
# strand                                                         -                                             -
# match                                                  GTTGGCAGA                                     GATGGCAGA
# variant                                                rs4976977                                     rs4976977
#
# this is the reverse complement PCM for MA0671.1.  note 23864 vs 2399 at position 2:
#
# A  [ 6439  2399   110     0    29  2539 32864  8043  7378 ]
# C  [ 9135 12688   180     0  3257 32864  5911 10855  8725 ]
# G  [ 9039  1356  2230 32864 32864   186  9695  6134  8810 ]
# T  [ 8251 32864 32864     0     3   132 11279  7832  7951 ]
