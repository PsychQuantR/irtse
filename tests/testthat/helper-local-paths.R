posix_local_path_pattern <- function() {
  roots <- c("Users", "home", "tmp", "var", "private", "opt", "Volumes",
             "srv", "etc", "mnt", "usr", "root")
  paste0("(^|[^A-Za-z0-9._/:-])", "/(", paste(roots, collapse = "|"), ")(/|$)")
}

additional_local_path_patterns <- function() {
  c(posix_local_path_pattern(),
    paste0("(^|[^A-Za-z])[A-Za-z]", ":", "/"))
}

# Files the scan for local paths can read as text. Serialized R files (the
# lazy-load databases and .rds files of an installed package) are compressed,
# so a pattern of three or four bytes meets them by chance; a version-control
# folder is not part of the package and holds paths of the machine it is on.
scannable_files <- function(paths) {
  serialized <- grepl("[.](rdb|rdx|rds)$", paths, ignore.case = TRUE)
  version_control <- grepl("(^|/)[.]git(/|$)", paths)
  paths[!serialized & !version_control]
}
