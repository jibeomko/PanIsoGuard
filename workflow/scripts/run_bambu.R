#!/usr/bin/env Rscript
# Run Bambu over a pre-computed (shared) BAM and write extended_annotations.gtf.
# Usage: run_bambu.R <bam> <ref_gtf> <genome_fa> <outdir> [ncore]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("usage: run_bambu.R <bam> <ref_gtf> <genome_fa> <outdir> [ncore]")
bam <- args[1]; ref_gtf <- args[2]; genome <- args[3]; outdir <- args[4]
ncore <- if (length(args) >= 5) as.integer(args[5]) else 1L

suppressMessages(library(bambu))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

annotations <- prepareAnnotations(ref_gtf)
se <- bambu(reads = bam, annotations = annotations, genome = genome, ncore = ncore)
writeBambuOutput(se, path = outdir)   # writes <outdir>/extended_annotations.gtf (+ counts)
