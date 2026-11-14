posix_local_path_pattern <- function() {
  roots <- c("Users", "home", "tmp", "var", "private", "opt", "Volumes",
             "srv", "etc", "mnt", "usr", "root")
  paste0("(^|[^A-Za-z0-9._/:-])", "/(", paste(roots, collapse = "|"), ")(/|$)")
}

additional_local_path_patterns <- function() {
  c(posix_local_path_pattern(),
    "(^|[^A-Za-z0-9._/:-])/[A-Za-z0-9._~-]+/[A-Za-z0-9._~/-]+",
    paste0("(^|[^A-Za-z0-9._/:-])/[A-Za-z][A-Za-z0-9._~-]*",
           "( [A-Za-z0-9._~-]+){1,3}/[^[:space:]]"),
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
