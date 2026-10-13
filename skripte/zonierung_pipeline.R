# =============================================================
# zonierung_pipeline.R
# Speicherort: Bachelorarbeit Pipeline/skripte/
# Benötigte Pakete (in load_packages_pipeline.R):
#   Seurat, ggplot2, patchwork, FNN (get.knn), dbscan (frNN)
#
# Alle Plots sind auskommentiert und werden nicht gespeichert.
# Jeder Plot-Block beginnt mit "# ---- PLOT: <Beschreibung> ----".
# Zum Aktivieren den Block markieren und in RStudio Strg+Shift+C drücken;
# der Plot erscheint dann im Plots-Fenster.
# =============================================================


# -------------------------------------------------------------
# Zonierung für ein Dataset
#   dataset:        normalisiertes Seurat-Objekt (NULL, wenn
#                   _zoniert.rds schon existiert)
#   dateiname:      z.B. "M1_Visium"
#   path_to_visium: Ordner mit den <Probe>_Visium-Unterordnern
# Rückgabe: das zonierte Seurat-Objekt
# -------------------------------------------------------------
zonierung <- function(dataset, dateiname, path_to_visium, num_zone = 8) {
  
  rds_file <- file.path(path_to_visium, dateiname,
                        paste0(dateiname, "_zoniert.rds"))
  message(dateiname, " -> ", rds_file)
  
  # ---- Bereits zoniert: nur laden
  if (file.exists(rds_file)) {
    return(readRDS(rds_file))
  }
  
  if (is.null(dataset)) {
    stop("Kein Dataset übergeben und keine Datei gefunden: ", rds_file)
  }
  
  #### pericentral-categorized Zonation Score ####
  # Landmark-Gene aus dem Paper
  LmG_pericentral <- c("CYP3A4", "ADH1B", "CYP1A2", "CYP2E1", "APOA2", "APOC1",
                       "ADH4", "ADH1A", "APOH", "AMBP", "GSTA2", "ADH1C",
                       "SLCO1B3", "AOX1", "APOA5", "DCXR", "RBP4", "OAT",
                       "CYP2C19", "GC")
  LmG_periportal  <- c("SERPINA1", "APOA1", "ALB", "C7", "NNMT", "HAMP", "ALDOB",
                       "ASS1", "CYP2A7", "MGP", "A2M", "FXYD2", "CCL21", "HAL",
                       "IGFBP2", "SDS", "AQP1", "CYP2A6", "FBLN1", "PTGDS")
  
  expressionsmatrix <- LayerData(dataset, assay = "Spatial", layer = "data")
  
  # nicht alle Gene in der Expressionsmatrix
  LmG_pericentral <- intersect(LmG_pericentral, rownames(expressionsmatrix))
  LmG_periportal  <- intersect(LmG_periportal,  rownames(expressionsmatrix))
  
  sum_pericentral <- Matrix::colSums(expressionsmatrix[LmG_pericentral, , drop = FALSE])
  sum_periportal  <- Matrix::colSums(expressionsmatrix[LmG_periportal,  , drop = FALSE])
  
  # portal -> 1, central -> 0; NA wenn beide Summen 0
  dataset$zonation_score <- ifelse(
    (sum_periportal + sum_pericentral) == 0,
    NA,
    sum_periportal / (sum_periportal + sum_pericentral)
  )
  
  #### Zonen über Quantile zuordnen ####
  thresh <- quantile(dataset$zonation_score, (0:num_zone) / num_zone, na.rm = TRUE)
  
  dataset$liver_Zone <- NA
  scores <- dataset$zonation_score
  
  for (i in seq_len(num_zone)) {
    # finde Spots zwischen thresh[i] und thresh[i+1] (Obergrenze ausgeschlossen,
    # wie im Original: der Spot mit dem höchsten Score bekommt dadurch keine Zone)
    spots <- names(scores)[!is.na(scores) & (scores >= thresh[i]) & (scores < thresh[i + 1])]
    dataset$liver_Zone[spots] <- i
  }
  
  #### Salt-and-pepper-Rauschen entfernen (Median im Radius) ####
  coords <- as.matrix(GetTissueCoordinates(dataset)[, c("x", "y")])
  
  # Abstand zum nächsten Nachbarn = Spot-Abstand (Visium: 100 µm)
  nn <- get.knn(coords, k = 1)
  step_size <- 100 / min(nn$nn.dist[, 1])
  coor_um <- coords * step_size
  
  radius_um <- 150
  nn_r <- frNN(coor_um, eps = radius_um)
  
  zone_index_med <- dataset$liver_Zone
  ind_rel_spots  <- which(dataset$liver_Zone > 0)
  
  for (i in ind_rel_spots) {
    ind_for_med <- intersect(nn_r$id[[i]], ind_rel_spots)
    ind_for_med <- setdiff(ind_for_med, i)
    if (length(ind_for_med) == 0) next
    zone_index_med[i] <- round(median(dataset$liver_Zone[ind_for_med], na.rm = TRUE))
  }
  
  dataset$med_liver_zone <- zone_index_med
  
  #### Seurat-Objekt speichern ####
  saveRDS(dataset, file = rds_file)
  
  dataset
}


# -------------------------------------------------------------
# Hauptfunktion: alle Proben zonieren
# -------------------------------------------------------------
zonierung_pipeline <- function(path_to_visium,
                               save_plots_to,
                               proben = c("M1","M2","M3","M4","M5","M6","M7","M8",
                                          "P2","P3","P6","P7","P14","P17","P18","P21"),
                               num_zone = 8) {
  
  datasets     <- list()
  summary_list <- list()
  
  for (i in proben) {
    filename <- paste0(i, "_Visium")
    rds_file <- file.path(path_to_visium, filename, paste0(filename, "_zoniert.rds"))
    
    # Rohdaten nur laden, wenn noch nicht zoniert (spart viel Zeit)
    if (file.exists(rds_file)) {
      dataset <- NULL
    } else {
      dataset <- Load10X_Spatial(data.dir = file.path(path_to_visium, filename),
                                 filename = paste0(i, "_filtered_feature_bc_matrix.h5"))
      dataset <- NormalizeData(dataset)
    }
    
    dataset <- zonierung(dataset, filename, path_to_visium, num_zone)
    
    ## UMIs pro Spot
    dataset$UMIs_per_bin <- Matrix::colSums(GetAssayData(dataset, layer = "counts"))
    summary_list[[i]] <- summary(dataset$UMIs_per_bin)
    
    datasets[[i]] <- dataset
  }
  
  # UMI-Zusammenfassung pro Probe ausgeben
  for (i in names(summary_list)) {
    cat("\n====================\n", i, "\n")
    print(summary_list[[i]])
  }
  
  # ---- PLOT: Geglättete Leberzonen (med_liver_zone) aller 16 Proben räumlich, 4 x 4 Raster ----
  # plots <- lapply(names(datasets), function(i) {
  #   SpatialFeaturePlot(subset(datasets[[i]], !is.na(liver_Zone)), features = "med_liver_zone",
  #                      pt.size.factor = 2) + ggtitle(paste("Liver Zones -", i))
  # })
  # print(wrap_plots(plots, ncol = 4, nrow = 4, guides = "collect") & theme(legend.position = "bottom"))
  
  # ---- PLOT: Zonation Score jeder Probe räumlich (0 = zentral, 1 = portal) ----
  # for (i in names(datasets)) {
  #   print(SpatialFeaturePlot(datasets[[i]], features = "zonation_score", pt.size.factor = 2) +
  #           ggtitle(paste("Zonation Score -", i)))
  # }
  
  # ---- PLOT: Histogramm des Zonation Scores jeder Probe mit den Quantil-Grenzen der Zonen (rot gestrichelt) ----
  # for (i in names(datasets)) {
  #   score  <- datasets[[i]]$zonation_score
  #   thresh <- quantile(score, (0:num_zone) / num_zone, na.rm = TRUE)
  #   print(ggplot(data.frame(score = score), aes(x = score)) +
  #           geom_histogram(bins = 100) +
  #           geom_vline(xintercept = thresh, colour = "red", linetype = "dashed", linewidth = 1) +
  #           theme_classic() + ggtitle(i))
  # }
  
  # ---- PLOT: Violinplot der UMIs pro Spot für jede Probe ----
  # for (i in names(datasets)) {
  #   print(VlnPlot(datasets[[i]], features = "UMIs_per_bin", pt.size = 0.1) + ggtitle(i))
  # }
  
  # ---- PLOT: Räumliche Expression der Clock-Gene in einer Beispielprobe (M2) ----
  # print(SpatialFeaturePlot(datasets[["M2"]], features = c("CRY1", "CRY2", "PER1", "PER2", "RORA", "CLOCK")))
  
  invisible(datasets)
}