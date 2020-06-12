#!/usr/bin/env Rscript

# update.R
# Last modified: 2020-06-12 10:46:37 (CEST)
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

gather_metadata <- function(subtype, META2ACC, METADATA) {
  ACCs <- META2ACC$Acc[META2ACC$Subtype == subtype]
  ACCs <- vapply(strsplit(ACCs, ".", fixed = TRUE), function(x) x[1], character(1))
  ACCs <- unique(ACCs)
  if (any(!ACCs %in% names(METADATA))) {
    ACCs2 <- ACCs[!ACCs %in% names(METADATA)]
    ACCs <- ACCs[ACCs %in% names(METADATA)]
    ACCs2 <- tibble(
      Id = NA, Source = NA,
      `Nucleotide Accession` = NA,
      Start = NA, Stop = NA, Strand = NA, Protein = ACCs2,
      `Protein Name` = NA, Organism = NA, Strain = NA, Assembly = NA
    )
    if (!length(ACCs))
      out <- ACCs2
    else
      out <- rbind(do.call(rbind, METADATA[ACCs]), ACCs2)
  } else {
    x <- METADATA[ACCs]
    if (length(x) && sum(vapply(x, nrow, integer(1)))) {
      out <- do.call(rbind, x)
    } else {
      out <- tibble(
        Id = NA, Source = NA,
        `Nucleotide Accession` = NA,
        Start = NA, Stop = NA, Strand = NA, Protein = names(x),
        `Protein Name` = NA, Organism = NA, Strain = NA, Assembly = NA
      )
    }
  }
  out[!is.na(out$Protein), ]
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

  repSeqs <- suppressMessages(readr::read_tsv(paste0("rep-sequences/repSequences", let, ".txt")))
  readr::write_tsv(repSeqs, paste0(app, "/data/repseqs.tsv"))

  MDa <- structure(lapply(
    list.files(meta, full.names = TRUE),
    function(x) suppressMessages(readr::read_tsv(x, comment=">"))
  ), names = list.files(meta))

  saveRDS(MDa, paste0(app, "/data/metadata.RDS"))

  meta2acc <- suppressMessages(readr::read_delim(
    paste0("classifications/tcd", tolower(let), "-allBLASThits.classified.txt"),
    " ", col_names = FALSE
  ))
  colnames(meta2acc) <- c("Acc", "Subtype", "Identity", "Length", "Coverage")
  meta2acc <- meta2acc[meta2acc$Identity == 100 & meta2acc$Coverage == "COMPLETE", ]
  saveRDS(meta2acc, paste0(app, "/data/metadata2acc.RDS"))

  groupSubLens <- lapply(namesA, function(x) unlist(lapply(x, function(y) nrow(gather_metadata(y, meta2acc, MDa)))))

  for (i in seq_along(groupSubLens)) {
    names(groupSubLens[[i]]) <- namesA[[i]]
    groupSubLens[[i]] <- sort(groupSubLens[[i]], decreasing = TRUE)
    groupSubLens[[i]] <- c(
      groupSubLens[[i]][repSeqs$Rep_identifier[repSeqs$Subtype == names(groupSubLens)[i]]],
      groupSubLens[[i]][!names(groupSubLens[[i]]) %in% repSeqs$Rep_identifier[repSeqs$Subtype == names(groupSubLens)[i]]]
    )
  }

  for (i in seq_along(groupsA)) {
    groupsA[[i]] <- groupsA[[i]][names(groupSubLens[[i]])]
  }

  for (i in seq_along(namesA)) {
    namesA[[i]] <- names(groupSubLens[[i]])
  }

  for (i in seq_along(namesA)) {
    namesA[[i]] <- structure(paste0(names(namesA)[i], ".", seq_along(namesA[[i]])), names = namesA[[i]])
    names(groupsA[[i]]) <- unname(namesA[[i]])
  }

  groupsAaligned <- groupsA

  groupsA <- lapply(groupsA, fix_seqs)

  metaAll <- lapply(namesA, function(x) lapply(names(x), function(y) gather_metadata(y, meta2acc, MDa)))
  for (i in seq_along(metaAll)) {
    names(metaAll[[i]]) <- names(namesA[[i]])
  }
  for (i in seq_along(metaAll)) {
    for (j in seq_along(metaAll[[i]])) {
      metaAll[[i]][[j]]$SeqId <- names(metaAll[[i]])[j]
    }
  }
  for (i in seq_along(metaAll)) {
    metaAll[[i]] <- do.call(rbind, metaAll[[i]])
  }
  metaAll <- do.call(rbind, metaAll)
  saveRDS(metaAll, paste0(app, "/data/ALL-metadata.RDS"))

  saveRDS(groupsA, paste0(app, "/data/ALL-sequences.RDS"))

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

  saveRDS(namesA, paste0(app, "/data/ALL-names.RDS"))

  fix_update_date(paste0(app, "/global.R"))

  unlink(list.files(paste0("app-", let, "/downloads"), full.names = TRUE))

}

message("Updating app-A..")
update_app("A")
message("Updating app-B..")
update_app("B")

