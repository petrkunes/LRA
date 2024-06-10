.onLoad <- function(libname, pkgname) {
  # Define a function to install packages from URLs if they are not already installed
  install_from_url <- function(package_name, url) {
    if (!requireNamespace(package_name, quietly = TRUE)) {
      install.packages(url, repos = NULL, type = "source")
    }
  }

  # Install required packages from custom URLs
  install_from_url("disqover", "https://github.com/MartinTheuerkauf/disqover/raw/8c9bd9bc08515e3f1a8b63464951c728cc673d0b/disqover_0.9.13.tar.gz")
}
