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
