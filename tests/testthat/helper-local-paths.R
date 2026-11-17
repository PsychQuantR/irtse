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
# A path whose first folder is not one of the twelve roots listed above and
# whose second component begins with a space is not found; with a listed root
# the first pattern finds it. Forms that occur in public packages and are not
# found (the test "Path forms that the scan does not find" runs each through
# the patterns): a file URL (file: and three slashes), a path written directly
# after a compiler flag (-I or -L), an entry of a search path that follows a
# colon after a name of two or more letters (bin: followed by the next path;
# after a single letter the drive pattern finds it), a path in the home folder
# written with a
# tilde and a single file name (~ and the name), and one written with a tilde
# and a user name.
additional_local_path_patterns <- function() {
  c(posix_local_path_pattern(),
    "(^|[^A-Za-z0-9._/:-])/[A-Za-z0-9._~-]+/[A-Za-z0-9._~/-]+",
    "(^|[^A-Za-z0-9._/:<-])/[A-Za-z0-9][^/\n]* [^/\n]*/[^[:space:]]",
    paste0("(^|[^A-Za-z])[A-Za-z]", ":", "/"))
}

# `paths` with the leading `root` removed. A path that is not below `root` is
# returned as it is.
relative_to_root <- function(paths, root) {
  prefix <- paste0(sub("/+$", "", root), "/")
  below <- startsWith(paths, prefix)
  paths[below] <- substring(paths[below], nchar(prefix) + 1L)
  paths
}

# Installed lazy-load databases and RDS files are compressed, and can match
# a short path pattern by chance. Source files are scanned regardless of suffix.
# The default preserves installed-package callers; source scans opt in explicitly.
# With `root`, what belongs to version control is decided from the path below
# that root, so that the name of a folder above the tree does not count.
# testthat writes tests/testthat/testthat-problems.rds and the folder
# tests/testthat/_problems after a failing run. They hold the output of the
# failures, which can name paths of the machine the tests ran on, and the next
# run would report them; they are not package files, .Rbuildignore keeps them
# out of the built package, and the scans skip them.
scannable_files <- function(paths, installed = TRUE, root = NULL) {
  inside <- if (is.null(root)) paths else relative_to_root(paths, root)
  serialized <- installed & grepl("[.](rdb|rdx|rds)$", inside, ignore.case = TRUE)
  version_control <- grepl("(^|/)[.]git(/|$)", inside)
  leftover <- grepl("(^|/)tests/testthat/(testthat-problems[.]rds$|_problems(/|$))", inside)
  paths[!serialized & !version_control & !leftover]
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
  paths <- scannable_files(tree_files(package_root), installed = !in_source,
                           root = package_root)
  # The plain listing of the R folder has no hidden entry and nothing below a
  # folder, so there is nothing for the selection to remove from it.
  if (in_source) paths <- c(paths, list.files(source_r, full.names = TRUE))
  unique(paths)
}

# The files read by the scan of the source tree at `source_root`. Where there
# is no source tree, they are the files of the installed package at
# `package_root`.
source_scan_files <- function(package_root, source_root) {
  if (!dir.exists(file.path(source_root, "R"))) {
    return(scannable_files(tree_files(package_root), installed = TRUE,
                           root = package_root))
  }
  paths <- c(file.path(source_root, c("DESCRIPTION", "README.md", "NEWS.md")),
             unlist(lapply(c("R", "man", "tests"), function(dir) {
               tree_files(file.path(source_root, dir))
             }), use.names = FALSE))
  scannable_files(paths, installed = FALSE, root = source_root)
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
