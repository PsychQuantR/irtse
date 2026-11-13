additional_local_path_patterns <- function() {
  c(paste0("/", c("tmp", "var", "private", "opt", "Volumes"), "/"),
    paste0("(^|[^A-Za-z])[A-Za-z]", ":", "/"))
}
