posix_local_path_pattern <- function() {
  roots <- c("Users", "home", "tmp", "var", "private", "opt", "Volumes",
             "srv", "etc", "mnt", "usr", "root")
  paste0("(^|[^A-Za-z0-9._/:-])", "/(", paste(roots, collapse = "|"), ")(/|$)")
}

additional_local_path_patterns <- function() {
  c(posix_local_path_pattern(),
    paste0("(^|[^A-Za-z])[A-Za-z]", ":", "/"))
}
