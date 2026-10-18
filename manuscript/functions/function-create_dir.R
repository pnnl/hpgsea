
create_data_dir <- function() {
  data_dir <- "manuscript/data"

  if (!dir.exists(data_dir)) {
    invisible(dir.create(data_dir))
  }

  return(data_dir)
}

create_figures_dir <- function() {
  figures_dir <- "manuscript/figures"

  if (!dir.exists(figures_dir)) {
    invisible(dir.create(figures_dir))
  }

  return(figures_dir)
}
