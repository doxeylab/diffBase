#!/usr/bin/env Rscript

# update.R
# Last modified: 2020-06-11 00:09:43 (CEST)
# BJM Tremblay

# possible tree building code:
# library(ape)
# aln <- read.alignment("file.fa", "fasta")
# tree <- bionj(dist.alignment(aln, matrix = "similarity"))
# tree <- ladderize(midpoint.root(tree))

UPDATE_DATE <- as.character(Sys.Date())

message("--- Running Diff-base update script (", UPDATE_DATE, ") ---")
message("")

make_names <- function(x, y, z) {
  names(z) <- paste0(x, y, ".", seq_len(length(z)))
  z
}

fix_update_date <- function(f) {
  l <- readr::read_lines(f)
  l[grep("^LAST_UPDATE_DATE <- function()", l)] <- paste0(
    "LAST_UPDATE_DATE <- function() \"",
    UPDATE_DATE, "\""
  )
  readr::write_lines(l, f)
}

fix_seqs <- function(x) {
  AAStringSet(toupper(gsub("-", "", as.character(x))))
}

suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(seqinr))
suppressPackageStartupMessages(library(phytools))

#-------------------------------------------------------------------------------

update_app <- function(app) {

  group <- paste0("groups-", app)
  let <- app
  app <- paste0("app-", let)
  meta <- paste0("metadata-", let)

  filesA <- list.files(group, full.names = TRUE)
  fileNamesA <- list.files(group)
  fileNamesA <- vapply(strsplit(fileNamesA, ".", fixed = TRUE), function(x) x[1], character(1))

  groupsA <- lapply(filesA, readAAStringSet)

  names(groupsA) <- fileNamesA

  namesA <- lapply(groupsA, names)
  names(namesA) <- fileNamesA

  for (i in seq_along(namesA)) {
    namesA[[i]] <- structure(paste0(names(namesA)[i], ".", seq_along(namesA[[i]])), names = namesA[[i]])
    names(groupsA[[i]]) <- unname(namesA[[i]])
  }

  groupsAaligned <- groupsA

  groupsA <- lapply(groupsA, fix_seqs)

  saveRDS(groupsA, paste0(app, "/data/ALL-sequences.RDS"))

  # namesA <- lapply(groupsA, names)

  MDa <- structure(lapply(
    list.files(meta, full.names = TRUE),
    function(x) suppressMessages(readr::read_tsv(x, comment=">"))
  ), names = list.files(meta))

  saveRDS(MDa, paste0(app, "/data/metadata.RDS"))

  groupsAflat <- groupsA[[1]]
  for (i in seq_along(groupsA)[-1]) groupsAflat <- c(groupsAflat, groupsA[[i]])

  groupsAalignedflat <- groupsAaligned[[1]]
  for (i in seq_along(groupsAaligned)[-1]) groupsAalignedflat <- c(groupsAalignedflat, groupsAaligned[[i]])

  writeXStringSet(groupsAalignedflat, paste0(app, "/data/ALL-sequences-aligned.fa"))

  writeXStringSet(groupsAflat, paste0(app, "/downloads/ALL-sequences.fa"))
  writeXStringSet(groupsAflat, paste0(app, "/blastdb/ALL-sequences.fa"))

  setwd(paste0(app, "/blastdb"))
  system("makeblastdb -in ALL-sequences.fa -dbtype prot -blastdb_version 4",
    ignore.stdout = TRUE, ignore.stderr = TRUE)
  setwd("../..")

  alnA <- read.alignment(paste0(app, "/data/ALL-sequences-aligned.fa"), "fasta")
  treeA <- dist.alignment(alnA)
  treeA <- hclust(treeA, method = "average")
  treeA <- as.phylo(treeA)

  saveRDS(treeA, paste0(app, "/data/tree.RDS"))

  repSeqs <- readr::read_tsv(paste0("rep-sequences/repSequences", let, ".txt"))
  readr::write_tsv(repSeqs, paste0(app, "/data/repseqs.tsv"))

  saveRDS(namesA, paste0(app, "/data/ALL-names.RDS"))

  fix_update_date(paste0(app, "/global.R"))

  unlink(list.files(paste0("app-", let, "/downloads"), full.names = TRUE))

  meta2acc <- suppressMessages(readr::read_delim(
    paste0("classifications/tcd", tolower(let), "-allBLASThits.classified.txt"),
    " ", col_names = FALSE
  ))
  colnames(meta2acc) <- c("Acc", "Subtype", "Identity", "Length", "Coverage")
  meta2acc <- meta2acc[meta2acc$Identity == 100 & meta2acc$Coverage == "COMPLETE", ]
  saveRDS(meta2acc, paste0(app, "/data/metadata2acc.RDS"))

}

message("Updating app-A..")
update_app("A")
message("Updating app-B..")
update_app("B")

