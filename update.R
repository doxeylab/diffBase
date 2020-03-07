#!/usr/bin/env Rscript

UPDATE_DATE <- as.character(Sys.Date())

message("--- Running Diff-base update script (", UPDATE_DATE, ") ---")
message("")

if (as.logical(toupper(Sys.getenv("FORCE_UPDATE")))) {
  message("Detected FORCE_UPDATE, deleting existing MD5 checksums.")
  message("")
  unlink(".diff-base.groups.md5.RDS")
}

if (!all(c("app-A", "app-B") %in% list.files())) {
  stop(
    "Couldn't find app directories. Make sure you run this script ",
    "from within the diff-base directory."
  )
}

UPDATE_A <- FALSE
UPDATE_B <- FALSE

message("Checking whether apps need updating..")

MD5_A <- tools::md5sum(list.files("groups-A", full.names = TRUE))
MD5_B <- tools::md5sum(list.files("groups-B", full.names = TRUE))
MD5_ALL <- list(A = MD5_A, B = MD5_B)
if (!file.exists(".diff-base.groups.md5.RDS")) {
  UPDATE_A <- TRUE
  UPDATE_B <- TRUE
  message("  Couldn't find old MD5 checksums. Forcing update.")
  message("")
} else {
  MD5_OLD <- readRDS(".diff-base.groups.md5.RDS")
  if (length(MD5_ALL$A) != length(MD5_OLD$A)) {
    message("")
    stop(
      "The number of groups in groups-A appears to have changed.\n",
      "The update script cannot handle this kind of update: please\n",
      "update the source code manually and delete the\n",
      ".diff-base.groups.md5.RDS file before rerunning this script."
    )
  }
  if (length(MD5_ALL$B) != length(MD5_OLD$B)) {
    message("")
    stop(
      "The number of groups in groups-B appears to have changed.\n",
      "The update script cannot handle this kind of update: please\n",
      "update the source code manually and delete the\n",
      ".diff-base.groups.md5.RDS file before rerunning this script."
    )
  }
  if (any(MD5_ALL$A != MD5_OLD$A)) {
    message("  Detected changes in groups-A.")
    UPDATE_A <- TRUE
  }
  if (any(MD5_ALL$B != MD5_OLD$B)) {
    message("  Detected changes in groups-B.")
    UPDATE_B <- TRUE
  }
  message("")
}

if (!UPDATE_A && !UPDATE_B) {
  message("")
  message("Everything appears to be up-to-date.")
  q()
}

make_names_A <- function(x) {
  structure(names(Alist[[x]]), names = paste0(x, ".", 1:length(Alist[[x]])))
}

fix_update_date <- function(f) {
  l <- readr::read_lines(f)
  l[grep("^LAST_UPDATE_DATE <- function()", l)] <- paste0(
    "LAST_UPDATE_DATE <- function() \"",
    UPDATE_DATE, "\""
  )
  readr::write_lines(l, f)
}

#-------------------------------------------------------------------------------

if (UPDATE_A) {

  message("Updating app-A..")

  suppressPackageStartupMessages(library(Biostrings))

  gA <- readAAStringSet("groups-A/groupA.fa")
  gB <- readAAStringSet("groups-A/groupB.fa")
  gC <- readAAStringSet("groups-A/groupC.fa")
  gD <- readAAStringSet("groups-A/groupD.fa")
  gE <- readAAStringSet("groups-A/groupE.fa")
  gF <- readAAStringSet("groups-A/groupF.fa")

  Alist <- list(
    A = gA, B = gB, C = gC, D = gD, E = gE, F = gF
  )

  Anames <- list(
    A = make_names_A("A"),
    B = make_names_A("B"),
    C = make_names_A("C"),
    D = make_names_A("D"),
    E = make_names_A("E"),
    F = make_names_A("F")
  )

  MD <- structure(lapply(
    list.files("metadata-A", full.names = TRUE),
    function(x) suppressMessages(readr::read_tsv(x, skip = 2))
  ), names = list.files("metadata-A"))

  MD_NAMES <- names(MD)
  AnamesFlat <- unlist(Anames, use.names = FALSE)
  # if (length(MD_NAMES) != length(AnamesFlat)) {
  #   stop(
  #     "The number of metadata files does not match the number\n",
  #     "of group-A sequences."
  #   )
  # }
  # if (anyNA(pmatch(MD_NAMES, AnamesFlat, nomatch = NA))) {
  #   stop(
  #     "Found mismatches between metadata names and group-A\n",
  #     "sequences names."
  #   )
  # }

  saveRDS(MD, "app-A/data/metadata.RDS")

  saveRDS(Anames, "app-A/data/ALL-names.RDS")

  for (i in seq_along(Anames)) {
    saveRDS(
      Anames[[i]],
      paste0("app-A/data/", names(Anames)[i], "-names.RDS")
    )
    saveRDS(
      Alist[[i]],
      paste0("app-A/data/", names(Anames)[i], "-sequences.RDS")
    )
  }

  for (i in seq_along(Anames)) {
    names(Alist[[i]]) <- names(Anames[[i]])
  }

  Avec <- c(Alist$A, Alist$B, Alist$C, Alist$D, Alist$E, Alist$F)
  Avec2 <- AAStringSet(gsub("-", "", as.character(Avec), fixed = TRUE))

  writeXStringSet(Avec2, "app-A/data/ALL-sequences.fa")
  writeXStringSet(Avec2, "app-A/downloads/ALL-sequences.fa")
  writeXStringSet(Avec2, "app-A/blastdb/ALL-sequences.fa")
  saveRDS(Avec2, "app-A/data/ALL-sequences-AAStringSet.RDS")
  saveRDS(as.character(Avec2), "app-A/data/ALL-sequences.RDS")

  Alist2 <- lapply(
    Alist,
    function(x) AAStringSet(gsub("-", "", as.character(x), fixed = TRUE))
  )
  for (i in seq_along(Anames)) {
    writeXStringSet(
      Alist2[[i]],
      paste0("app-A/downloads/", names(Anames)[i], "-sequences.fa")
    )
  }

  if (as.logical(toupper(Sys.getenv("NO_UPDATE_TREE")))) {
    message("  Detected NO_UPDATE_TREE")
  } else {
    suppressPackageStartupMessages(library(phangorn))
    Amat <- as.matrix(dist.ml(phyDat(as.matrix(Avec), type = "AA")))
    Atree <- ape::as.phylo(hclust(as.dist(Amat)))
    saveRDS(Atree, "app-A/data/tree.RDS")
  }

  setwd("app-A/blastdb")
  if (system(
        "makeblastdb -blastdb_version 4 -dbtype prot -in ALL-sequences.fa",
        ignore.stdout = TRUE, ignore.stderr = TRUE
      )) {
    setwd("../../")
    stop("makeblastdb command failed")
  }
  setwd("../../")

  fix_update_date("app-A/global.R")

  message("")

}

make_names_B <- function(x) {
  structure(names(Blist[[x]]), names = paste0(x, ".", 1:length(Blist[[x]])))
}

#-------------------------------------------------------------------------------

if (UPDATE_B) {

  message("Updating app-B..")

  suppressPackageStartupMessages(library(Biostrings))

  gA <- readAAStringSet("groups-B/groupA.fa")
  gB <- readAAStringSet("groups-B/groupB.fa")
  gC <- readAAStringSet("groups-B/groupC.fa")
  gD <- readAAStringSet("groups-B/groupD.fa")
  gE <- readAAStringSet("groups-B/groupE.fa")
  gF <- readAAStringSet("groups-B/groupF.fa")
  gG <- readAAStringSet("groups-B/groupG.fa")
  gH <- readAAStringSet("groups-B/groupH.fa")
  gI <- readAAStringSet("groups-B/groupI.fa")
  gJ <- readAAStringSet("groups-B/groupJ.fa")
  gK <- readAAStringSet("groups-B/groupK.fa")
  gL <- readAAStringSet("groups-B/groupL.fa")

  Blist <- list(
    A = gA, B = gB, C = gC, D = gD, E = gE, F = gF,
    G = gG, H = gH, I = gI, J = gJ, K = gK, L = gL
  )

  Bnames <- list(
    A = make_names_B("A"),
    B = make_names_B("B"),
    C = make_names_B("C"),
    D = make_names_B("D"),
    E = make_names_B("E"),
    F = make_names_B("F"),
    G = make_names_B("G"),
    H = make_names_B("H"),
    I = make_names_B("I"),
    J = make_names_B("J"),
    K = make_names_B("K"),
    L = make_names_B("L")
  )

  MD <- structure(lapply(
    list.files("metadata-B", full.names = TRUE),
    function(x) suppressMessages(readr::read_tsv(x, skip = 2))
  ), names = list.files("metadata-B"))

  MD_NAMES <- names(MD)
  BnamesFlat <- unlist(Bnames, use.names = FALSE)
  # if (length(MD_NAMES) != length(BnamesFlat)) {
  #   stop(
  #     "The number of metadata files does not match the number\n",
  #     "of group-A sequences."
  #   )
  # }
  # if (anyNA(pmatch(MD_NAMES, BnamesFlat, nomatch = NA))) {
  #   stop(
  #     "Found mismatches between metadata names and group-B\n",
  #     "sequences names."
  #   )
  # }

  saveRDS(MD, "app-B/data/metadata.RDS")

  saveRDS(Bnames, "app-B/data/ALL-names.RDS")

  for (i in seq_along(Bnames)) {
    saveRDS(
      Bnames[[i]],
      paste0("app-B/data/", names(Bnames)[i], "-names.RDS")
    )
    saveRDS(
      Blist[[i]],
      paste0("app-B/data/", names(Bnames)[i], "-sequences.RDS")
    )
  }

  for (i in seq_along(Bnames)) {
    names(Blist[[i]]) <- names(Bnames[[i]])
  }

  Bvec <- c(
    Blist$A, Blist$B, Blist$C, Blist$D, Blist$E, Blist$F,
    Blist$G, Blist$H, Blist$I, Blist$J, Blist$K, Blist$L
  )
  Bvec2 <- AAStringSet(gsub("-", "", as.character(Bvec), fixed = TRUE))

  writeXStringSet(Bvec2, "app-B/data/ALL-sequences.fa")
  writeXStringSet(Bvec2, "app-B/downloads/ALL-sequences.fa")
  writeXStringSet(Bvec2, "app-B/blastdb/ALL-sequences.fa")
  saveRDS(Bvec2, "app-B/data/ALL-sequences-AAStringSet.RDS")
  saveRDS(as.character(Bvec2), "app-B/data/ALL-sequences.RDS")

  Blist2 <- lapply(
    Blist,
    function(x) AAStringSet(gsub("-", "", as.character(x), fixed = TRUE))
  )
  for (i in seq_along(Bnames)) {
    writeXStringSet(
      Blist2[[i]],
      paste0("app-B/downloads/", names(Bnames)[i], "-sequences.fa")
    )
  }

  if (as.logical(toupper(Sys.getenv("NO_UPDATE_TREE")))) {
    message("  Detected NO_UPDATE_TREE")
  } else {
    suppressPackageStartupMessages(library(phangorn))
    Bmat <- as.matrix(dist.ml(phyDat(as.matrix(Bvec), type = "AA")))
    Btree <- ape::as.phylo(hclust(as.dist(Bmat)))
    saveRDS(Btree, "app-B/data/tree.RDS")
  }

  setwd("app-B/blastdb")
  if (system(
        "makeblastdb -blastdb_version 4 -dbtype prot -in ALL-sequences.fa",
        ignore.stdout = TRUE, ignore.stderr = TRUE
      )) {
    setwd("../../")
    stop("makeblastdb command failed")
  }
  setwd("../../")

  fix_update_date("app-B/global.R")

  message("")

}

#-------------------------------------------------------------------------------

message("All done.")
saveRDS(MD5_ALL, ".diff-base.groups.md5.RDS")
