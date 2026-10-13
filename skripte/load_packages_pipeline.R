# =============================================================
# load_packages_pipeline.R
# Speicherort: Bachelorarbeit Pipeline/skripte/
# Installiert fehlende Pakete und lädt alle, die die Pipeline braucht.
# =============================================================

load_packages_pipeline <- function(){
  
  #### Paketlisten ####
  
  # CRAN
  cran_packages <- c(
    "Seurat",      # Visium-Objekte, SpatialFeaturePlot, AggregateExpression
    "hdf5r",       # Load10X_Spatial (.h5)
    "Matrix",      # colSums auf Sparse-Matrizen
    "FNN",         # get.knn (Zonierung)
    "dbscan",      # frNN (Zonierung)
    "readxl",      # read_xlsx (Talamanca-Zeittabelle)
    "CePa",        # read.gct (GTEx)
    "MetaCycle",   # meta2d (über tauFisher::run_metacycle)
    "nnet",        # multinom-Modell in tauFisher
    "boot",        # Bootstrap (Statistik)
    "ggplot2",
    "patchwork",
    "ggrepel",
    "scales",
    "colorspace",
    "gridExtra",
    "magick",
    "tidyr",
    "dplyr"        # zuletzt laden, damit dplyr::select nicht von AnnotationDbi überdeckt wird
  )
  
  # Bioconductor
  bioc_packages <- c(
    "AnnotationDbi",  # select() für die Ensembl -> Symbol Übersetzung
    "org.Hs.eg.db"    # humane Gen-Annotation
  )
  
  # GitHub
  github_packages <- c(
    "tauFisher" = "micnngo/tauFisher"
  )
  
  ist_installiert <- function(pkg) requireNamespace(pkg, quietly = TRUE)
  
  
  #### CRAN-Pakete installieren ####
  
  missing_cran <- cran_packages[!vapply(cran_packages, ist_installiert, logical(1))]
  
  if (length(missing_cran) > 0) {
    message("Fehlende CRAN-Pakete: ", paste(missing_cran, collapse = ", "))
    install.packages(missing_cran, dependencies = TRUE)
  } else {
    message("Alle CRAN-Pakete sind bereits installiert.")
  }
  
  
  #### Bioconductor-Pakete installieren ####
  
  missing_bioc <- bioc_packages[!vapply(bioc_packages, ist_installiert, logical(1))]
  
  if (length(missing_bioc) > 0) {
    if (!ist_installiert("BiocManager")) install.packages("BiocManager")
    message("Fehlende Bioconductor-Pakete: ", paste(missing_bioc, collapse = ", "))
    BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
  } else {
    message("Alle Bioconductor-Pakete sind bereits installiert.")
  }
  
  
  #### GitHub-Pakete installieren ####
  
  missing_github <- names(github_packages)[!vapply(names(github_packages), ist_installiert, logical(1))]
  
  if (length(missing_github) > 0) {
    
    if (!ist_installiert("remotes"))  install.packages("remotes")
    if (!ist_installiert("pkgbuild")) install.packages("pkgbuild")
    
    # GitHub-Pakete werden kompiliert und brauchen Build-Tools
    # (Windows: Rtools passend zur R-Version)
    if (!pkgbuild::has_build_tools(debug = FALSE)) {
      stop(
        "\n========================================\n",
        "FEHLER: Build-Tools fehlen!\n",
        "Fehlende GitHub-Pakete: ", paste(missing_github, collapse = ", "), "\n",
        "Bitte Rtools passend zur R-Version installieren\n",
        "und R/RStudio neu starten.\n",
        "========================================\n"
      )
    }
    
    for (pkg in missing_github) {
      message("Installiere ", pkg, " von GitHub (", github_packages[[pkg]], ") ...")
      remotes::install_github(github_packages[[pkg]], dependencies = TRUE,
                              build_vignettes = FALSE)
    }
  } else {
    message("Alle GitHub-Pakete sind bereits installiert.")
  }
  
  
  #### Pakete laden ####
  # Reihenfolge: Bioconductor zuerst, dann GitHub, dann CRAN (dplyr zuletzt)
  
  all_packages <- c(bioc_packages, names(github_packages), cran_packages)
  failed <- character(0)
  
  for (pkg in all_packages) {
    ok <- tryCatch({
      suppressPackageStartupMessages(library(pkg, character.only = TRUE))
      TRUE
    }, error = function(e) {
      message("✗ ", pkg, ": ", conditionMessage(e))
      FALSE
    })
    if (ok) message("✓ ", pkg) else failed <- c(failed, pkg)
  }
  
  
  #### Ergebnis ####
  
  message("\n========================================")
  if (length(failed) == 0) {
    message("✓ Alle Pakete sind installiert und geladen.")
  } else {
    message("⚠ Probleme mit: ", paste(failed, collapse = ", "))
    message("Bitte die Fehlermeldungen oben prüfen.")
  }
  message("========================================")
  
  invisible(length(failed) == 0)
}