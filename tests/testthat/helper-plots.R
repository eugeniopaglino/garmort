get_test_plot_dir <- function() {
  env_dir <- Sys.getenv("GARMORT_TEST_PLOTS_DIR", unset = "")
  if (nzchar(env_dir)) {
    dir.create(env_dir, recursive = TRUE, showWarnings = FALSE)
    return(env_dir)
  }
  
  # When running devtools::test() locally, working directory is tests/testthat
  if (file.exists("../../DESCRIPTION")) {
    local_dir <- normalizePath("../../test-diagnostics", winslash = "/", mustWork = FALSE)
    dir.create(local_dir, recursive = TRUE, showWarnings = FALSE)
    return(local_dir)
  }
  
  # Fallback during R CMD check so the check directory stays clean
  tmp_dir <- file.path(tempdir(), "garmort-diagnostics")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  tmp_dir
}