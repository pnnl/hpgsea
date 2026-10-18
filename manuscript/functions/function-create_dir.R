
create_data_dir <- function() {
  data_dir <- "simulation/data"

  if (!dir.exists(data_dir)) {
    invisible(dir.create(data_dir))
  }

  return(data_dir)
}

create_figures_dir <- function() {
  figures_dir <- "simulation/figures"

  if (!dir.exists(figures_dir)) {
    invisible(dir.create(figures_dir))
  }

  return(figures_dir)
}
