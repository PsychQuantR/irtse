posix_local_path_pattern <- function() {
  roots <- c("Users", "home", "tmp", "var", "private", "opt", "Volumes",
             "srv", "etc", "mnt", "usr", "root")
  paste0("(^|[^A-Za-z0-9._/:-])", "/(", paste(roots, collapse = "|"), ")(/|$)")
}

# This scan is a heuristic aid, not a proof that local paths are absent.
# For first folders with spaces it misses names beginning with anything other
# than a letter or digit, and paths whose spaces are tabs or non-ASCII spaces.
# It also misses paths broken across lines and relative paths, and, in the branch
# for folders with spaces, a path that directly follows "<": that branch does not
# start at "</", so that a closing HTML tag followed by text and another closing
# tag (as in the help index of an installed package) is not taken for a path.
# Two divisions with a space before each slash and none after are taken for a
# path (a false alarm); respacing the expression avoids it.
# A path whose second component begins with a space is not found.
additional_local_path_patterns <- function() {
  c(posix_local_path_pattern(),
    "(^|[^A-Za-z0-9._/:-])/[A-Za-z0-9._~-]+/[A-Za-z0-9._~/-]+",
    "(^|[^A-Za-z0-9._/:<-])/[A-Za-z0-9][^/\n]* [^/\n]*/[^[:space:]]",
    paste0("(^|[^A-Za-z])[A-Za-z]", ":", "/"))
}

# Installed lazy-load databases and RDS files are compressed, and can match
# a short path pattern by chance. Source files are scanned regardless of suffix.
# The default preserves installed-package callers; source scans opt in explicitly.
scannable_files <- function(paths, installed = TRUE) {
  serialized <- installed & grepl("[.](rdb|rdx|rds)$", paths, ignore.case = TRUE)
  version_control <- grepl("(^|/)[.]git(/|$)", paths)
  paths[!serialized & !version_control]
}

# References to the private tracker, written in pieces so that this file does
# not match its own pattern.
private_tracker_pattern <- function() {
  paste0("Psych", "Quant/[^\r\n]*", "#", "[0-9]+")
}

# The patterns of the scan of the shipped package.
shipped_scan_patterns <- function() {
  c(posix_local_path_pattern(), private_tracker_pattern())
}

# The patterns of the scan of the source tree.
source_scan_patterns <- function() {
  c(paste0("[A-Za-z]", ":", "\\\\"), private_tracker_pattern(),
    paste0("/", "issues", "/", "[0-9]+"), additional_local_path_patterns())
}

# Every file below `root`, hidden files included.
tree_files <- function(root) {
  list.files(root, recursive = TRUE, full.names = TRUE, all.files = TRUE)
}

# The files read by the scan of the shipped package: every file of the
# installed (or loaded) package at `package_root` and, when `source_root` is a
# source tree (it has an R folder), the files of that folder.
shipped_scan_files <- function(package_root, source_root) {
  source_r <- file.path(source_root, "R")
  in_source <- dir.exists(source_r)
  paths <- tree_files(package_root)
  if (in_source) paths <- c(paths, list.files(source_r, full.names = TRUE))
  unique(scannable_files(paths, installed = !in_source))
}

# The files read by the scan of the source tree at `source_root`. Where there
# is no source tree, they are the files of the installed package at
# `package_root`.
source_scan_files <- function(package_root, source_root) {
  in_source <- dir.exists(file.path(source_root, "R"))
  paths <- if (in_source) {
    c(file.path(source_root, c("DESCRIPTION", "README.md", "NEWS.md")),
      unlist(lapply(c("R", "man", "tests"), function(dir) {
        tree_files(file.path(source_root, dir))
      }), use.names = FALSE))
  } else {
    tree_files(package_root)
  }
  scannable_files(paths, installed = !in_source)
}

# The files among `paths` whose bytes match one of `patterns`.
files_matching <- function(paths, patterns) {
  hit <- vapply(paths, function(path) {
    bytes <- readBin(path, "raw", n = file.info(path)$size)
    any(vapply(patterns, function(pattern) length(grepRaw(pattern, bytes)) > 0L,
               logical(1)))
  }, logical(1))
  unname(paths[hit])
}
