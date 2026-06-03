#!/usr/bin/env Rscript
# Bambu with an explicit NDR (novel discovery rate) so it does not invoke the online
# NDR-recommendation step (which needs BiocManager/network). NDR=0.5 is a permissive
# discovery setting, recorded in the protocol for reproducibility. Writes
# <outdir>/extended_annotations.gtf (reference + novel BambuTx transcripts); the
# protocol then keeps only the novel BambuTx rows.
#
# Usage: run_bambu_ndr.R <bam> <ref_gtf> <genome_fa> <outdir> [ncore]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("usage: run_bambu_ndr.R <bam> <ref_gtf> <genome_fa> <outdir> [ncore]")
bam <- args[1]; ref_gtf <- args[2]; genome <- args[3]; outdir <- args[4]
ncore <- if (length(args) >= 5) as.integer(args[5]) else 1L

suppressMessages(library(bambu))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
ann <- prepareAnnotations(ref_gtf)
se <- bambu(reads = bam, annotations = ann, genome = genome, ncore = ncore,
            NDR = 0.5, quant = FALSE)
writeToGTF(se, file.path(outdir, "extended_annotations.gtf"))
cat("bambu wrote", length(se), "transcripts\n")
