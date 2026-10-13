# =============================================================
# tauFisher_pipeline.R
# Speicherort: Bachelorarbeit Pipeline/skripte/
#
# Speicherorte (aus config.R):
#   path_to_save_model -> ergebnisse/modelle/modelle  (Modelle, Trainingsdaten, MetaCycle)
#   path_to_save_pred  -> ergebnisse/modelle/pred     (Vorhersagen + ZT-Lookups)
#
# Drei Modelle:
#   JTK      MetaCycle und Genauswahl auf den ersten n_train_metacycle (=170)
#            Zeitpunkten; Training auf allen vollständigen Spalten (wie im
#            Original) oder mit training_nur_170 = TRUE nur auf den ersten 170
#   LS       wie JTK
#   LS_full  MetaCycle (nur LS), Genauswahl und Training auf allen Zeitpunkten
#   Fehlende Zeitpunkte bekommt MetaCycle als NA-Spalten (JTK braucht ein
#   gleichmäßiges Zeitraster); das Training nutzt nur die vollständigen.
#   save_plot_file_to  -> ergebnisse/plots (derzeit ungenutzt: alle Plots sind
#                         auskommentiert und werden nicht gespeichert)
#
# Alle Plots sind auskommentiert und werden nicht gespeichert.
# Jeder Plot-Block beginnt mit "# ---- PLOT: <Beschreibung> ----".
# Zum Aktivieren den Block markieren und in RStudio Strg+Shift+C drücken;
# der Plot erscheint dann im Plots-Fenster.
# =============================================================

# function übersetzt EN IDs in GenSymbole
#          prüft auf doppelte Gene
#          bennent Spaltennamen um
vorbereitung_training_table<- function(training_data) {
  #### Pre Processing ####
  # wir wollen keine doubles
  # Take the average expression if there are non-unique genes
  training_data$ID <- rownames(training_data)
  
  # Ensembl-Versionsnummer entfernen
  # ENSG00000290825.2 -> ENSG00000290825
  training_data$ID <- sub(
    "\\..*$",
    "",
    training_data$ID
  )
  
  #### Ensembl umwandlung ####
  library(org.Hs.eg.db)
  
  mapping <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = training_data$ID,
    keytype = "ENSEMBL",
    columns = c("ENSEMBL", "SYMBOL")
  )
  
  # Spalten umbenennen, damit der Rest des Codes unverändert bleibt
  colnames(mapping) <- c("ensembl_gene_id", "hgnc_symbol")
  
  training_data$ID <- mapping$hgnc_symbol[
    match(training_data$ID, mapping$ensembl_gene_id)
  ]
  
  print("training_data after maaping")
  print(training_data$ID)
  
  # nicht gemappte Gene entfernen
  training_data <- training_data[
    !is.na(training_data$ID) & training_data$ID != "",
  ]
  
  #### training vorbereitung ####
  training<- training_data[stats::complete.cases(training_data[,-1]), ] %>%
    dplyr::mutate(ID = toupper(ID)) %>%
    dplyr::group_by(ID) %>%
    dplyr::summarise_all(mean)  %>%
    as.data.frame()
  
  # Set the rownames to be the genes and remove the ID column 
  rownames(training) <- training$ID
  training <- training[,-1]
  
  return(training)
}

# leichte anpassung von find_periodic_genes von tauFisher 
# es wird geprüft ob alle gefundenen gene auch im prediction datenset vorhanden sind
find_periodic_genes <- function(input_file,
                                test_genes,
                                per_method=c("LS", "JTK"),
                                thres=1.01){
  
  cc_genes = toupper(c("Bmal1", "Arntl", "Dbp", "Nr1d1", "Nr1d2","Per1", "Per2", "Per3", "Cry1",  "Cry2"))
  
  test_genes = toupper(test_genes)
  
  if(any(per_method %in% "LS")){
    
    per_LS <- utils::read.delim(input_file) %>%
      dplyr::filter(round(LS_period) == 24, LS_pvalue < thres) %>%
      dplyr::arrange(LS_pvalue) %>%
      dplyr::mutate(CycID = toupper(CycID))
    
    genes_LS = per_LS[per_LS$CycID %in% test_genes, ] %>%
      dplyr::top_n(-10, LS_pvalue) %>%
      dplyr::select(CycID) %>% c(.)
    genes_LS = unique(c(genes_LS$CycID, cc_genes[cc_genes %in% per_LS$CycID]) )
    
  } else{ genes_LS = NULL}
  
  if(any(per_method %in% "JTK")){
    
    per_JTK <- utils::read.delim(input_file) %>%
      dplyr::filter(round(JTK_period) == 24, JTK_pvalue < thres) %>%
      dplyr::arrange(JTK_pvalue) %>%
      dplyr::mutate(CycID = toupper(CycID))
    
    genes_JTK = per_JTK[per_JTK$CycID %in% test_genes, ] %>%
      dplyr::top_n(-10, JTK_pvalue) %>%
      dplyr::select(CycID) %>% c(.)
    genes_JTK = unique(c(genes_JTK$CycID, cc_genes[cc_genes %in% per_JTK$CycID]))
    
  } else{ genes_JTK = NULL}
  
  
  if(any(!(per_method %in% c("JTK", "LS")))){
    print("Not a valid method included. Must be at least one of JTK or LS.")
  }
  genes_JTK <- genes_JTK[genes_JTK %in% test_genes]
  genes_LS <- genes_LS[genes_LS %in% test_genes]
  
  return(list(LS = genes_LS, JTK = genes_JTK))
}

# Feste Sortierreihenfolge fuer Sample-IDs in ALLEN Plots:
# M aufsteigend (M1..M8), P in der Reihenfolge P7,P6,P14,P18,P17,P2,P21,P3.
sample_level_order <- c(
  "M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8",
  "P7", "P6", "P14", "P18", "P17", "P2", "P21", "P3"
)

# Gibt eine Tabelle vollständig aus (alle Zeilen und Spalten)
print_full <- function(x, label = NULL) {
  if (is.null(x) || nrow(x) == 0) {
    message(if (!is.null(label)) paste0(label, ": ") else "", "keine Zeilen zum Anzeigen.")
    return(invisible(NULL))
  }
  if (!is.null(label)) cat(label, "\n")
  print(dplyr::as_tibble(x), n = Inf, width = Inf)
}


# Hilfsfunktion: wandelt einen Vektor von Sample-IDs (z.B. "M3" oder "P14")
# in einen Factor mit der festen sample_level_order um, sodass ggplot
# (facet_wrap, Farbskalen, Achsen) immer in dieser Reihenfolge sortiert.
# IDs, die nicht in sample_level_order vorkommen, werden alphabetisch ans
# Ende angehaengt, damit nichts verloren geht.
order_samples <- function(x) {
  extra <- sort(setdiff(unique(as.character(x)), sample_level_order))
  factor(x, levels = c(sample_level_order, extra))
}


tauFisher_pipeline <- function(path_to_xlsx,
                               path_to_gct,
                               path_to_visium,
                               path_to_save_model,
                               path_to_save_pred,
                               save_plot_file_to,
                               bulk,
                               zonen,
                               bin,
                               n_train_metacycle = 170,
                               training_nur_170 = FALSE){
  
  #### Speicherorte ####
  path_to_model <- path_to_save_model   # Modelle, Trainingsdaten, MetaCycle
  path_to_pred  <- path_to_save_pred    # Vorhersagen + ZT-Lookups
  dir.create(path_to_model, recursive = TRUE, showWarnings = FALSE)
  dir.create(path_to_pred,  recursive = TRUE, showWarnings = FALSE)
  
  # MetaCycle-Ein-/Ausgabe (vorher im Arbeitsverzeichnis unter "tauFisher/")
  metacycle_dir <- file.path(path_to_model, "metacycle")
  dir.create(metacycle_dir, recursive = TRUE, showWarnings = FALSE)
  out_file    <- file.path(metacycle_dir, "Dataset_full_adj.txt")
  out_dir     <- file.path(metacycle_dir, "Dataset_full_adj_meta2d")
  meta2d_file <- file.path(out_dir, "meta2d_Dataset_full_adj.txt")
  
  # zusätzlicher LS-Lauf auf der GESAMTEN Trainingstabelle (eigener Ordner,
  # damit er die 170-Spalten-Ausgabe nicht überschreibt)
  metacycle_dir_full <- file.path(path_to_model, "metacycle_full")
  dir.create(metacycle_dir_full, recursive = TRUE, showWarnings = FALSE)
  out_file_full    <- file.path(metacycle_dir_full, "Dataset_full_adj.txt")
  out_dir_full     <- file.path(metacycle_dir_full, "Dataset_full_adj_meta2d")
  meta2d_file_full <- file.path(out_dir_full, "meta2d_Dataset_full_adj.txt")
  
  model_file_JTK <- file.path(path_to_model, "model_talamanca_JTK.rds")
  model_file_LS  <- file.path(path_to_model, "model_talamanca_LS.rds")
  model_file_LS_full <- file.path(path_to_model, "model_talamanca_LS_full.rds")
  
  #### Import training Data ####
  timeTabelle <- readxl::read_xlsx(path_to_xlsx)
  
  expressionTabelle <- read.gct(path_to_gct)
  
  cat("ANZAHL SAMPLES IN GCT:", ncol(expressionTabelle))
  
  ## expression Tabelle vorbereiten: 
  # GTEX.ZZPU.0426.SM.5GZYH -> ZZPU -> jeweiliger Uhrzeit die in timeTabelle steht gerundet
  colnames(expressionTabelle) <- sub("^[^.]+\\.([^.]+)\\..*$", "\\1", colnames(expressionTabelle))
  
  ids <- colnames(expressionTabelle)
  
  hours <- timeTabelle$hour[match(ids, timeTabelle$SUBJ.ID)]
  
  expressionTabelle <- expressionTabelle[, !is.na(hours), drop = FALSE]
  
  colnames(expressionTabelle) <- round(hours[!is.na(hours)])
  
  ## Spalten sortieren und Umbennen 
  # 8 -> ZT_8_REP1 ...
  ord <- order(as.numeric(colnames(expressionTabelle)))
  expressionTabelle <- expressionTabelle[, ord]
  
  hours <- as.numeric(colnames(expressionTabelle))
  
  rep <- ave(hours, hours, FUN = seq_along)
  
  colnames(expressionTabelle) <- paste0("ZT_", hours, "_REP_", rep)
  
  expressionTabelle <- as.data.frame(expressionTabelle)
  
  training <- vorbereitung_training_table(expressionTabelle)
  rm(expressionTabelle)
  
  
  #### TRAINING DATA PREPROCESSING #####
  
  # ----------------------------------------------------------
  # Maximale Anzahl Replikate festlegen
  
  max_rep <- 21
  
  # Nur die ersten n Replikate pro ZT verwenden
  replicate_all <- as.numeric(
    sub("^ZT_\\d+_REP_(\\d+)$", "\\1", colnames(training))
  )
  
  keep_rep <- replicate_all <= max_rep
  
  training <- training[, keep_rep, drop = FALSE]
  
  # ----------------------------------------------------------
  # Zeit und Replikat aus Spaltennamen auslesen
  
  time <- as.numeric(sub("^ZT_(\\d+)_REP_\\d+$", "\\1", colnames(training)))
  replicate <- as.numeric(sub("^ZT_\\d+_REP_(\\d+)$", "\\1", colnames(training)))
  
  # Kontrolle
  cat("\nZeitpunkte vorhanden:\n")
  print(sort(unique(time)))
  
  cat("\nReplikate vorhanden:\n")
  print(sort(unique(replicate)))
  
  cat("\nAnzahl pro ZT:\n")
  print(table(time))
  
  cat("\nAnzahl pro Replikat:\n")
  print(table(replicate))
  
  
  # ----------------------------------------------------------
  # ZT24 entfernen -> 0 und 24 überschneiden sich 
  # ZT24 entspricht wieder ZT0.
  # Für einen 24-Stunden-Zyklus verwenden wir ZT0-ZT23.
  
  keep <- time < 24
  
  training <- training[, keep, drop = FALSE]
  
  time <- time[keep]
  replicate <- replicate[keep]
  
  
  # ----------------------------------------------------------
  # testen auf fehlende kombinationen
  # Erwartet: ZT 0-23 x REP 1-max_rep
  
  expected_combinations <- expand.grid(
    time = 0:23,
    replicate = 1:max_rep
  )
  
  observed_combinations <- unique(data.frame(time = time, replicate = replicate))
  
  missing_combinations <- dplyr::anti_join(
    expected_combinations,
    observed_combinations,
    by = c("time", "replicate")
  )
  
  if (nrow(missing_combinations) == 0) {
    cat("Keine Kombination fehlt.\n")
  } else {
    print(missing_combinations)
  }
  
  # Fehlende Samples als NA-Spalten hinzufügen für metacycle
  if (nrow(missing_combinations) > 0) {
    for (i in seq_len(nrow(missing_combinations))) {
      t <- missing_combinations$time[i]
      r <- missing_combinations$replicate[i]
      new_name <- paste0("ZT_", t, "_REP_", r)
      training[[new_name]] <- NA_real_
    }
  }
  
  
  # ----------------------------------------------------------
  # Zeit und Replikat nach dem Ergänzen erneut auslesen
  
  time <- as.numeric(sub("^ZT_(\\d+)_REP_\\d+$", "\\1", colnames(training)))
  replicate <- as.numeric(sub("^ZT_\\d+_REP_(\\d+)$", "\\1", colnames(training)))
  
  
  # ----------------------------------------------------------
  # Nach ZT und Replikat sortieren aka ZT_0_REP_1, ZT_0_REP_2,...,ZT_1_REP_1,...
  
  ord <- order(time, replicate)
  
  training <- training[, ord, drop = FALSE]
  
  time <- time[ord]
  replicate <- replicate[ord]
  
  
  # ----------------------------------------------------------
  # Adjusted time für tauFisher
  #
  # REP1: ZT0 -> 0  ... ZT23 -> 23
  # REP2: ZT0 -> 24 ... ZT23 -> 47
  # REP3: ZT0 -> 48 ... ZT23 -> 71
  
  time_adj <- time + 24 * (replicate - 1)
  
  
  # ----------------------------------------------------------
  # Nach adjusted time sortieren
  
  ord_adj <- order(time_adj)
  
  training <- training[, ord_adj, drop = FALSE]
  
  time_adj <- time_adj[ord_adj]
  time <- time[ord_adj]
  replicate <- replicate[ord_adj]
  
  
  # ----------------------------------------------------------
  # Spaltennamen auf adjusted time setzen
  
  colnames(training) <- time_adj
  
  
  # MetaCycle-Input-Dateien werden direkt vor run_metacycle() geschrieben,
  # jeweils mit genau den Spalten, die auch fürs Training benutzt werden.
  
  
  #### LAST STEPS ####
  
  train_adj = training
  colnames(train_adj) <- time_adj
  train_adj = train_adj[, order(as.numeric(colnames(train_adj)))]
  time_adj <- time_adj[order(as.numeric(time_adj))]
  
  keep_clean <- colSums(is.na(train_adj)) == 0
  train_adj_clean <- train_adj[, keep_clean, drop = FALSE]
  time_adj_clean  <- time_adj[keep_clean]
  replicate_clean <- replicate[keep_clean]
  
  # ----------------------------------------------------------
  # Für den JTK- und LS-Lauf nur die ersten n_train_metacycle Spalten verwenden
  # (für MetaCycle UND für das Training mit train_tauFisher)
  
  idx_170 <- seq_len(n_train_metacycle)
  # für MetaCycle: dieselben Zeitpunkte, fehlende als NA-Spalten, damit das
  # Zeitraster gleichmäßig bleibt (sonst rechnet meta2d kein JTK, nur LS)
  train_adj_170 <- train_adj[, idx_170, drop = FALSE]
  time_adj_170  <- time_adj[idx_170]
  keep_clean_170      <- keep_clean[idx_170]
  train_adj_clean_170 <- train_adj[, idx_170, drop = FALSE][, keep_clean_170, drop = FALSE]
  time_adj_clean_170  <- time_adj[idx_170][keep_clean_170]
  replicate_clean_170 <- replicate[idx_170][keep_clean_170]
  
  cat("Training LS_full:", ncol(train_adj_clean), "vollständige von", ncol(train_adj), "Spalten\n")
  
  print(table(time_adj_clean_170 %% 24))
  
  # wie im Original: die ersten 170 Spalten inkl. NA-Spalten speichern
  # (plots_pipeline.R entfernt die NA-Spalten selbst)
  saveRDS(train_adj_170, file = file.path(path_to_model, "train_adj_JTK.rds"))
  saveRDS(time_adj_170,  file = file.path(path_to_model, "time_adj_JTK.rds"))
  # Trainingsdaten des LS_full-Modells (alle vollständigen Spalten), für plots_pipeline.R
  saveRDS(train_adj, file = file.path(path_to_model, "train_adj_LS_full.rds"))
  saveRDS(time_adj,  file = file.path(path_to_model, "time_adj_LS_full.rds"))
  
  
  #### Import Prediction Data ####
  dateipfad <- path_to_visium
  
  daten <- sample_level_order
  
  prediction_data <- list()
  
  ## gesamtes Datenset
  
  if(bulk == TRUE){
    pb_matrix_bulk <- NULL
    
    for (i in daten){
      dataset <- paste0("/",i,"_Visium/",i,"_Visium")
      rds_file <- paste0(dateipfad,dataset, "_zoniert.rds")
      data <- readRDS(rds_file)
      
      data$sample <- i 
      
      pseudobulk <- AggregateExpression(
        data,
        group.by = "sample",
        assays = "Spatial",
        return.seurat = FALSE
      )
      
      # an Gesamtmatrix anhängen
      pb <- pseudobulk$Spatial
      colnames(pb) <- i
      
      if (is.null(pb_matrix_bulk)) {
        pb_matrix_bulk <- pb
      } else {
        pb_matrix_bulk <- cbind(pb_matrix_bulk, pb)
      }
    }
    
    prediction_data$bulk <- pb_matrix_bulk
    rm(pb_matrix_bulk)
  }
  
  
  ## einzelne Zonen der Datensets
  
  if(zonen == TRUE){
    pb_matrix_zone <- NULL
    
    for (i in daten){
      
      dataset <- paste0("/", i, "_Visium/", i, "_Visium")
      rds_file <- paste0(dateipfad, dataset, "_zoniert.rds")
      
      data <- readRDS(rds_file)
      
      # NA-Zonen entfernen
      data <- subset(data, subset = !is.na(med_liver_zone))
      
      for (zone in 1:8){
        
        # nur diese einzelne Zone subseten
        zone_data <- subset(data, subset = med_liver_zone == zone)
        
        zone_data$sample <- paste0(i, "_Zone_", zone)
        
        pseudobulk <- AggregateExpression(
          object = zone_data,
          group.by = "sample",
          assays = "Spatial",
          slot = "counts",
          return.seurat = FALSE
        )
        
        pb <- pseudobulk$Spatial
        colnames(pb) <- paste0(i, "_Zone_", zone)
        
        # Matrix anhängen
        if (is.null(pb_matrix_zone)) {
          pb_matrix_zone <- pb
        } else {
          pb_matrix_zone <- cbind(pb_matrix_zone, pb)
        }
      }
    }
    
    prediction_data$zonen <- pb_matrix_zone
    rm(pb_matrix_zone)
  }
  
  
  ## einzelne Bins der Datensets
  
  if(bin == TRUE){
    pb_matrix_bin <- NULL
    
    for (i in daten){
      
      dataset <- paste0("/", i, "_Visium/", i, "_Visium")
      rds_file <- paste0(dateipfad, dataset, "_zoniert.rds")
      
      data <- readRDS(rds_file)
      
      # NA-Zonen entfernen
      data <- subset(data, subset = !is.na(med_liver_zone))
      
      # Count-Matrix aller Bins dieses Samples
      bin_counts <- GetAssayData(
        object = data,
        assay = "Spatial",
        layer = "counts"
      )
      
      # Samplenamen vor die originalen Visium-Barcodes setzen
      colnames(bin_counts) <- paste0(i, "_", colnames(bin_counts))
      
      # Matrix zusammenfügen
      if (is.null(pb_matrix_bin)) {
        pb_matrix_bin <- bin_counts
      } else {
        pb_matrix_bin <- cbind(pb_matrix_bin, bin_counts)
      }
    }
    
    prediction_data$bin <- pb_matrix_bin
    rm(pb_matrix_bin)
  }
  
  #### tauFisher Prediction ####
  
  ## Farben (werden nur von den auskommentierten Plots benutzt)
  zt_anchor_colours <- c(
    "0"  = "#E68A2E", "3"  = "#D9B54C", "6"  = "#8FAE7A", "9"  = "#4C8B8B",
    "12" = "#3D6FA8", "15" = "#6B4E8E", "18" = "#A6528E", "21" = "#D9636B",
    "24" = "#E68A2E"
  )
  zt_ramp <- colorRampPalette(zt_anchor_colours)
  zt_colours_24 <- zt_ramp(25)[1:24]
  names(zt_colours_24) <- 0:23
  zt_colours_continuous <- zt_ramp(100)
  
  zt_colour_scale <- function() {
    scale_fill_gradientn(
      colours = zt_colours_continuous,
      values = scales::rescale(0:23),
      limits = c(0, 23),
      breaks = seq(0, 23, 2),
      oob = scales::squish,
      name = "ZT"
    )
  }
  
  # Krankheitsstatus pro Sample
  status_of <- function(s) dplyr::case_when(
    grepl("^M", s)                     ~ "Healthy",
    s %in% c("P6", "P7", "P14", "P18") ~ "Steatosis",
    s == "P17"                         ~ "Severe steatosis",
    s %in% c("P2", "P21")              ~ "Fibrosis",
    s == "P3"                          ~ "Portal flow deprivation",
    TRUE                               ~ "Other"
  )
  status_cols <- c("Healthy" = "#4A90E2", "Steatosis" = "#2ca02c",
                   "Severe steatosis" = "#145214", "Fibrosis" = "#ff7f0e",
                   "Portal flow deprivation" = "#9467bd", "Other" = "grey50")
  
  group_colours <- c("M" = "#1f77b4", "P" = "#d62728")
  
  hour_labels <- function(x) {
    ifelse(x == 0, "0", paste0(ifelse(x > 0, "+", ""), x, "h"))
  }
  
  # Signifikanzsterne aus adjustierten p-Werten
  sig_stars <- function(p) case_when(
    is.na(p)  ~ "",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "ns"
  )
  
  # Projektion von Testdaten in den tauFisher-PCA-Raum
  # (exakt die Vorverarbeitung aus test_tauFisher(); nur für den PCA-Plot)
  tf_project <- function(tf, test_df) {
    f_diff  <- if (exists("create_DiffMatrix")) create_DiffMatrix else tauFisher:::create_DiffMatrix
    f_scale <- if (exists("scale_DiffMatrix"))  scale_DiffMatrix  else tauFisher:::scale_DiffMatrix
    rownames(test_df) <- toupper(rownames(test_df))
    genes <- tf$Genes[order(tf$Genes)]
    df <- t(data.frame(test_df[tf$Genes, , drop = FALSE]))
    df <- df[, order(colnames(df)), drop = FALSE]
    d  <- f_scale(f_diff(genes, log2(df + 1)))
    stats::predict(tf$PCA_embedding, newdata = d)
  }
  
  # PCA-Raum-Plot: Trainingsdaten, Entscheidungsregionen des multinom-Modells
  # und die projizierten Testproben. Wird unten nur aufgerufen, wenn der
  # Aufruf einkommentiert wird.
  plot_tf_pca_space <- function(tf, pred_data, predicted, type, method) {
    pca <- tf$PCA_embedding
    mod <- tf$Model
    
    proj <- tf_project(tf, pred_data)
    
    classify <- function(m) {
      d <- as.data.frame(m[, 1:2, drop = FALSE])
      names(d) <- c("PC1", "PC2")
      as.integer(as.character(predict(mod, newdata = d, type = "class"))) %% 24
    }
    zt_pred <- round(predicted) %% 24
    cat("PCA plot (", type, method, "): Nachklassifikation = test_tauFisher() bei",
        round(100 * mean(classify(proj) == zt_pred, na.rm = TRUE), 1), "% der Proben\n")
    
    train_t <- as.numeric(sub("_.*$", "", rownames(pca$x)))
    
    rng <- function(v) range(v, na.rm = TRUE) + c(-0.1, 0.1) * diff(range(v, na.rm = TRUE))
    x_rng <- rng(c(pca$x[, 1], proj[, 1]))
    y_rng <- rng(c(pca$x[, 2], proj[, 2]))
    grid <- expand.grid(PC1 = seq(x_rng[1], x_rng[2], length.out = 250),
                        PC2 = seq(y_rng[1], y_rng[2], length.out = 250))
    grid_hex <- zt_colours_24[as.character(classify(as.matrix(grid)))]
    
    ve <- round(100 * summary(pca)$importance[2, 1:2], 1)
    train_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                           ZT = factor(floor(train_t) %% 24, levels = 0:23))
    hour_lab <- aggregate(cbind(PC1, PC2) ~ ZT, data = train_df, FUN = mean)
    hour_lab[, c("PC1", "PC2")] <- hour_lab[, c("PC1", "PC2")] * 1.15
    
    genes_da  <- intersect(toupper(tf$Genes), toupper(rownames(pred_data)))
    expr_mat  <- as.matrix(pred_data)
    rownames(expr_mat) <- toupper(rownames(expr_mat))
    sample_id <- sub("^(M[0-9]+|P[0-9]+).*", "\\1", colnames(pred_data))
    test_df <- data.frame(
      ID = colnames(pred_data), PC1 = proj[, 1], PC2 = proj[, 2],
      Status  = factor(status_of(sample_id), levels = names(status_cols)),
      ZT_pred = factor(zt_pred, levels = 0:23),
      n_genes = colSums(expr_mat[genes_da, , drop = FALSE] > 0)
    )
    
    base <- ggplot() +
      geom_raster(data = grid, aes(PC1, PC2), fill = grid_hex, alpha = 0.25) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
      geom_point(data = train_df, aes(PC1, PC2, colour = ZT), size = 1.5, alpha = 0.8) +
      geom_text(data = hour_lab, aes(PC1, PC2, label = ZT, colour = ZT),
                size = 5, fontface = "bold", show.legend = FALSE) +
      scale_colour_manual(values = zt_colours_24, name = "ZT (training)", drop = FALSE) +
      coord_equal(expand = FALSE) +
      labs(x = paste0("PC1 (", ve[1], " %)"), y = paste0("PC2 (", ve[2], " %)"),
           caption = "Background: decision regions of the multinomial model") +
      theme_classic()
    
    print(base + labs(title = paste0("tauFisher PCA space (", method, ") - training data")))
    
    if (type == "bin") {
      set.seed(1)
      td <- test_df[sample.int(nrow(test_df), min(20000, nrow(test_df))), ]
      print(base +
              geom_point(data = td, aes(PC1, PC2, fill = ZT_pred), shape = 21, stroke = 0.1,
                         size = 1.2, alpha = 0.4) +
              scale_fill_manual(values = zt_colours_24, name = "Predicted ZT", drop = FALSE) +
              labs(title = paste0("tauFisher PCA space (", method, ") - bins by predicted ZT")))
      print(base +
              geom_point(data = td[order(td$n_genes, decreasing = TRUE), ],
                         aes(PC1, PC2, fill = n_genes), shape = 21, stroke = 0.1, size = 1.2, alpha = 0.5) +
              scale_fill_gradient(low = "red", high = "grey20", name = "Detected\npredictor genes") +
              labs(title = paste0("tauFisher PCA space (", method, ") - bins by detected predictor genes")))
    } else {
      p <- base +
        geom_point(data = test_df, aes(PC1, PC2, fill = Status), shape = 21,
                   size = if (type == "bulk") 3.5 else 2.5, alpha = 0.85) +
        scale_fill_manual(values = status_cols, name = "Disease status", drop = TRUE) +
        labs(title = paste0("tauFisher PCA space (", method, ") - projected ", type, " samples"))
      if (type == "bulk") p <- p + ggrepel::geom_text_repel(data = test_df, aes(PC1, PC2, label = ID),
                                                            size = 3, max.overlaps = Inf)
      print(p)
    }
    invisible(test_df)
  }
  
  
  # schreibt die MetaCycle-Input-Datei mit genau den übergebenen Spalten
  write_metacycle_input <- function(df, file) {
    utils::write.table(cbind(Genes = rownames(df), df), file = file,
                       sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  }
  
  for(type in names(prediction_data)){
    
    #### Prediction Data Preprocessing ####
    pred_data <- prediction_data[[type]]
    pred_data[is.na(pred_data)] <- 0
    
    # 11 ist gesuchte Uhrzeit (Platzhalter, tauFisher braucht eine "truth")
    pred_time <- rep(11, length(colnames(pred_data)))
    
    #### Modelle trainieren oder laden ####
    if(!(file.exists(model_file_JTK)) | !(file.exists(model_file_LS))){
      
      # MetaCycle + find_periodic_genes auf denselben Zeitpunkten wie das Training;
      # fehlende Zeitpunkte bleiben als NA drin (gleichmäßiges Raster für JTK)
      write_metacycle_input(train_adj_170, out_file)
      run_metacycle(df=train_adj_170,
                    timepoints=time_adj_170,
                    out_file=out_file,
                    out_dir=out_dir,
                    method=c("JTK","LS"),
                    parallelize = TRUE, nCores = 1 )
      
      # Don't set a threshold for p-value
      genes = find_periodic_genes(input_file=meta2d_file,
                                  test_genes=rownames(pred_data),
                                  per_method=c("JTK","LS"),
                                  thres=1)
      print(genes)
      
      # training_nur_170 = FALSE: wie im Original auf ALLEN vollständigen Spalten
      # training_nur_170 = TRUE:  nur auf den vollständigen unter den ersten 170
      if (training_nur_170) {
        tr_df <- train_adj_clean_170; tr_time <- time_adj_clean_170; tr_rep <- replicate_clean_170
      } else {
        tr_df <- train_adj_clean;     tr_time <- time_adj_clean;     tr_rep <- replicate_clean
      }
      max_rep_train <- max(tr_rep, na.rm = TRUE)
      cat("Training JTK/LS auf", ncol(tr_df), "Spalten, nrep =", max_rep_train, "\n")
      
      trained_JTK <- train_tauFisher(periodic_genes=genes$JTK,
                                     train_df=tr_df,
                                     train_time=tr_time,
                                     nrep=max_rep_train,
                                     numbasis=5)
      
      trained_LS <- train_tauFisher(periodic_genes=genes$LS,
                                    train_df=tr_df,
                                    train_time=tr_time,
                                    nrep=max_rep_train,
                                    numbasis=5)
      
      print("MODELLE WERDEN GESPEICHERT")
      saveRDS(trained_JTK, file = model_file_JTK)
      saveRDS(trained_LS, file = model_file_LS)
      
      # ---- PLOT: FDA-Kurven der tauFisher-Prädiktorgene (JTK) im Training, einmal nur Kurven, einmal mit Messpunkten ----
      # genes_used <- trained_JTK$Genes[trained_JTK$Genes %in% rownames(train_adj_clean)]
      # df_fda <- t(train_adj_clean[genes_used, ])
      # fda_curves <- get_FDAcurves(dat = df_fda, time = time_adj_clean, numbasis = 5)
      # fda_df <- fda_curves %>%
      #   pivot_longer(cols = -fda_time, names_to = "Gene", values_to = "Expression")
      # original_df <- train_adj_clean[genes_used, ] %>%
      #   as.data.frame() %>%
      #   mutate(Gene = rownames(.)) %>%
      #   pivot_longer(cols = -Gene, names_to = "Time", values_to = "Expression") %>%
      #   mutate(Time = as.numeric(Time))
      # print(ggplot(fda_df, aes(fda_time, Expression)) +
      #         geom_line(color = "steelblue", linewidth = 1) +
      #         facet_wrap(~ Gene, scales = "free_y", ncol = 4) +
      #         scale_x_continuous(breaks = seq(0, max(time_adj_clean), by = 4)) +
      #         labs(title = "FDA-Kurven der von tauFisher verwendeten Gene", x = "Zeit (h)", y = "Expression") +
      #         theme_classic())
      # print(ggplot() +
      #         geom_point(data = original_df, aes(Time, Expression), color = "black", size = 2) +
      #         geom_line(data = fda_df, aes(fda_time, Expression), color = "steelblue", linewidth = 1.2) +
      #         facet_wrap(~Gene, scales = "free_y", ncol = 4) +
      #         labs(title = "Genexpression und tauFisher FDA-Kurve", x = "Zeit (h)", y = "Expression") +
      #         theme_classic())
      
    } else {
      trained_JTK <- readRDS(model_file_JTK)
      trained_LS  <- readRDS(model_file_LS)
    }
    
    #### Zusätzliches LS-Modell aus der GESAMTEN Tabelle trainieren oder laden ####
    if(!file.exists(model_file_LS_full)){
      
      # MetaCycle + find_periodic_genes auf denselben Zeitpunkten wie das Training;
      # fehlende Zeitpunkte bleiben als NA drin (gleichmäßiges Raster)
      write_metacycle_input(train_adj, out_file_full)
      run_metacycle(df=train_adj,
                    timepoints=time_adj,
                    out_file=out_file_full,
                    out_dir=out_dir_full,
                    method="LS",
                    parallelize = TRUE, nCores = 1 )
      
      genes_full = find_periodic_genes(input_file=meta2d_file_full,
                                       test_genes=rownames(pred_data),
                                       per_method="LS",
                                       thres=1)
      print(genes_full)
      
      # Training auf ALLEN vollständigen Spalten
      trained_LS_full <- train_tauFisher(periodic_genes=genes_full$LS,
                                         train_df=train_adj_clean,
                                         train_time=time_adj_clean,
                                         nrep=max(replicate_clean, na.rm = TRUE),
                                         numbasis=5)
      
      print("MODELL LS (gesamte Tabelle) WIRD GESPEICHERT")
      saveRDS(trained_LS_full, file = model_file_LS_full)
      
    } else {
      trained_LS_full <- readRDS(model_file_LS_full)
    }
    
    #### predict circadian time JTK ####
    
    missing <- setdiff(trained_JTK$Genes, rownames(pred_data))
    if (length(missing) > 0) message("Fehlende Gene (JTK): ", paste(missing, collapse = ", "))
    cat("GENES_JTK:", trained_JTK$Genes, "\n")
    
    result_JTK <- test_tauFisher(chosen_genes=trained_JTK$Genes,
                                 test_df=pred_data,
                                 test_time=pred_time,
                                 pca_embedding=trained_JTK$PCA_embedding,
                                 trained_model=trained_JTK$Model)
    
    calc_error(result_JTK$truth, result_JTK$predicted)
    
    #### Predictions speichern (JTK) -> pred/ ####
    pred_result_df_JTK <- data.frame(
      ID = colnames(pred_data),
      ZT_predicted = result_JTK$predicted,
      ZT_truth = result_JTK$truth,
      method = "JTK",
      type = type,
      stringsAsFactors = FALSE
    )
    saveRDS(pred_result_df_JTK, file = file.path(path_to_pred, paste0(type, "_prediction_JTK.rds")))
    
    # Bulk-ZT Lookup (wird im "zonen"- und "bin"-Durchlauf gebraucht)
    if(type == "bulk"){
      bulk_zt_lookup_JTK <- setNames(result_JTK$predicted, colnames(pred_data))
      saveRDS(bulk_zt_lookup_JTK, file = file.path(path_to_pred, "bulk_zt_lookup_JTK.rds"))
    }
    
    # ---- PLOT: PCA-Raum von tauFisher (JTK) mit Entscheidungsregionen und projizierten Proben dieses Typs ----
    # plot_tf_pca_space(trained_JTK, pred_data, result_JTK$predicted, type, "JTK")
    
    # ---- PLOT: Oszillation der Prädiktorgene (JTK) im Training (Punkte + FDA-Kurve) und Heatmap nach Peak-Zeit; bei "bulk" zusätzlich die Visium-Bulk-Proben an ihrer vorhergesagten ZT ----
    # f_fda <- if (exists("get_FDAcurves")) get_FDAcurves else tauFisher:::get_FDAcurves
    # tr <- as.matrix(train_adj_clean)
    # rownames(tr) <- toupper(rownames(tr))
    # osz_genes <- sort(toupper(trained_JTK$Genes))
    # osz_genes <- osz_genes[osz_genes %in% rownames(tr)]
    # df_log <- log2(t(tr[osz_genes, , drop = FALSE]) + 1)
    # fda_long <- tidyr::pivot_longer(f_fda(df_log, time = time_adj_clean, numbasis = 5),
    #                                 -fda_time, names_to = "Gene", values_to = "Expr")
    # fda_long$Gene  <- toupper(fda_long$Gene)
    # fda_long$ZT    <- fda_long$fda_time %% 24
    # fda_long$Cycle <- floor(fda_long$fda_time / 24)
    # fda_long$z     <- ave(fda_long$Expr, fda_long$Gene, FUN = function(v) (v - mean(v)) / sd(v))
    # pts <- data.frame(ZT   = rep(time_adj_clean, times = length(osz_genes)) %% 24,
    #                   Gene = rep(osz_genes, each = nrow(df_log)),
    #                   Expr = as.vector(df_log))
    # hm <- aggregate(z ~ Gene + ZT_h, data = transform(fda_long, ZT_h = round(ZT) %% 24), FUN = mean)
    # gene_order <- names(sort(sapply(split(hm, hm$Gene), function(d) d$ZT_h[which.max(d$z)])))
    # lab <- setNames(osz_genes, osz_genes)
    # if (file.exists(meta2d_file)) {
    #   m2d <- utils::read.delim(meta2d_file)
    #   m2d$Gene <- toupper(m2d$CycID)
    #   jtk <- m2d[m2d$Gene %in% osz_genes & !duplicated(m2d$Gene), ]
    #   lab[jtk$Gene] <- sprintf("%s\npeak time ZT %.1f, p = %.1e",
    #                            jtk$Gene, jtk$JTK_adjphase %% 24, jtk$JTK_pvalue)
    # }
    # fda_long$Gene <- factor(fda_long$Gene, levels = gene_order)
    # pts$Gene      <- factor(pts$Gene,      levels = gene_order)
    # hm$Gene       <- factor(hm$Gene,       levels = rev(gene_order))
    # print(ggplot() +
    #         geom_point(data = pts, aes(ZT, Expr), colour = "black", size = 1.2, alpha = 0.5) +
    #         geom_line(data = fda_long, aes(ZT, Expr, group = Cycle), colour = "black", linewidth = 1) +
    #         facet_wrap(~ Gene, scales = "free_y", ncol = 4, labeller = as_labeller(lab)) +
    #         scale_x_continuous(limits = c(0, 24), breaks = seq(0, 24, 4)) +
    #         labs(title = "Oscillation of the tauFisher predictor genes (JTK) in the training data",
    #              x = "ZT (h)", y = "log2(expression + 1)") +
    #         theme_bw())
    # print(ggplot(hm, aes(ZT_h, Gene, fill = z)) +
    #         geom_tile(colour = "white") +
    #         scale_fill_gradient2(low = "#3D6FA8", mid = "white", high = "#D9636B", midpoint = 0, name = "z-score") +
    #         scale_x_continuous(breaks = seq(0, 23, 2), expand = c(0, 0)) +
    #         labs(title = "Predictor genes (JTK): FDA curves, z-scored, ordered by peak time", x = "ZT (h)", y = NULL) +
    #         theme_minimal() + theme(panel.grid = element_blank()))
    # if (type == "bulk") {
    #   te <- as.matrix(pred_data)
    #   rownames(te) <- toupper(rownames(te))
    #   g_te <- intersect(osz_genes, rownames(te))
    #   te_z <- t(scale(t(log2(te[g_te, , drop = FALSE] + 1))))
    #   te_df <- data.frame(Gene   = rep(g_te, times = ncol(te_z)),
    #                       Sample = rep(colnames(te_z), each = nrow(te_z)),
    #                       z      = as.vector(te_z))
    #   te_df$ZT     <- (round(result_JTK$predicted) %% 24)[match(te_df$Sample, colnames(pred_data))]
    #   te_df$Status <- factor(status_of(te_df$Sample), levels = names(status_cols))
    #   te_df$Gene   <- factor(te_df$Gene, levels = gene_order)
    #   print(ggplot() +
    #           geom_line(data = fda_long, aes(ZT, z, group = Cycle), colour = "black", linewidth = 1) +
    #           geom_point(data = te_df, aes(ZT, z, fill = Status), shape = 21, size = 2.8, alpha = 0.85,
    #                      position = position_jitter(width = 0.2, height = 0, seed = 1)) +
    #           facet_wrap(~ Gene, ncol = 4, labeller = as_labeller(lab)) +
    #           scale_fill_manual(values = status_cols, name = "Disease status", drop = TRUE) +
    #           scale_x_continuous(limits = c(-0.5, 24), breaks = seq(0, 24, 4)) +
    #           labs(title = "Visium bulk samples at their predicted ZT (JTK) vs. training curve",
    #                x = "ZT (h)", y = "z-score") +
    #           theme_bw() + theme(legend.position = "bottom"))
    # }
    
    
    ##### Bulk (JTK) #####
    if(bulk == TRUE & type == "bulk"){
      
      # ---- PLOT: Häufigkeit der vorhergesagten Bulk-ZT (JTK) als Balkendiagramm, getrennt für M, P und M + P ----
      # freq_plot <- data.frame(Sample = colnames(pred_data), ZT = round(result_JTK$predicted) %% 24) %>%
      #   mutate(Gruppe = ifelse(grepl("^M", Sample), "M", "P"))
      # freq_plot <- bind_rows(
      #   freq_plot %>% filter(Gruppe == "M") %>% count(ZT) %>% mutate(Panel = "M"),
      #   freq_plot %>% filter(Gruppe == "P") %>% count(ZT) %>% mutate(Panel = "P"),
      #   freq_plot %>% count(ZT) %>% mutate(Panel = "M + P")
      # ) %>%
      #   complete(Panel = c("M", "P", "M + P"), ZT = 0:23, fill = list(n = 0)) %>%
      #   mutate(Panel = factor(Panel, levels = c("M", "P", "M + P")), ZT = factor(ZT, levels = 0:23))
      # print(ggplot(freq_plot, aes(ZT, n, fill = ZT)) +
      #         geom_col(colour = "black") +
      #         facet_wrap(~Panel, ncol = 1) +
      #         scale_fill_manual(values = zt_colours_24, name = "ZT") +
      #         scale_x_discrete(breaks = as.character(seq(0, 23, by = 2))) +
      #         labs(title = "Häufigkeit der vorhergesagten ZT (JTK)", x = "Zeitgeberzeit (ZT)", y = "Häufigkeit") +
      #         theme_bw() + theme(legend.position = "bottom"))
      
      # ---- PLOT: Kreis-Plot (Polar) aller Bulk-Samples an ihrer vorhergesagten ZT (JTK), farbig nach M/P ----
      # plot_df <- data.frame(Sample = order_samples(colnames(pred_data)),
      #                       Zeit = result_JTK$predicted %% 24) %>%
      #   mutate(Typ = ifelse(grepl("^M", Sample), "M", "P")) %>%
      #   dplyr::arrange(Sample)
      # print(ggplot(plot_df, aes(Zeit, Sample, color = Typ)) +
      #         geom_point(size = 4) +
      #         ggrepel::geom_text_repel(aes(label = Sample), size = 3, max.overlaps = Inf, show.legend = FALSE) +
      #         scale_colour_manual(values = group_colours, name = "Gruppe") +
      #         coord_polar(theta = "x", start = 0) +
      #         scale_x_continuous(limits = c(0, 24), breaks = 0:23, labels = paste0(0:23, "h")) +
      #         scale_y_discrete(expand = expansion(add = c(3, 1)), limits = rev(levels(plot_df$Sample))) +
      #         labs(title = "Predicted ZT (JTK) pro Sample", x = "", y = "") +
      #         theme_minimal() + theme(axis.text.y = element_blank()))
      
      # ---- PLOT: Trainings-Expression der Prädiktorgene (JTK, Replikat 1) über die Zeit, alle zusammen und einzeln ----
      # plot_data_long <- as.data.frame(train_adj_clean[trained_JTK$Genes, replicate_clean == 1, drop = FALSE]) %>%
      #   mutate(Gene = rownames(.)) %>%
      #   pivot_longer(cols = -Gene, names_to = "Time", values_to = "Expression") %>%
      #   mutate(Time = as.numeric(Time))
      # p_all <- ggplot(plot_data_long, aes(Time, Expression, color = Gene, group = Gene)) +
      #   geom_line(linewidth = 1) + geom_point(size = 2) +
      #   labs(title = "Alle Gene", x = "Zeit", y = "Expression") + theme_classic()
      # p_facet <- ggplot(plot_data_long, aes(Time, Expression, group = Gene)) +
      #   geom_line(linewidth = 1, color = "steelblue") + geom_point(size = 2, color = "steelblue") +
      #   facet_wrap(~Gene, scales = "free_y") +
      #   labs(x = "Zeit", y = "Expression") + theme_classic()
      # print(p_all / p_facet)
    }
    
    ##### Zonen (JTK) #####
    if(zonen == TRUE & type == "zonen"){
      
      if(!exists("bulk_zt_lookup_JTK", inherits = FALSE)){
        bulk_zt_lookup_JTK <- readRDS(file.path(path_to_pred, "bulk_zt_lookup_JTK.rds"))
      }
      
      # Zonen-ZT Lookup (wird im "bin"-Durchlauf gebraucht)
      zonen_zt_lookup_JTK <- setNames(result_JTK$predicted, colnames(pred_data))
      saveRDS(zonen_zt_lookup_JTK, file = file.path(path_to_pred, "zonen_zt_lookup_JTK.rds"))
      
      # Zonen-Vorhersagen als Tabelle: Zone, Sample_ID, Gruppe (M/P), Delta zu Zone 1
      zone_df <- data.frame(
        Sample = colnames(pred_data),
        ZT = result_JTK$predicted
      )
      zone_df$Zone <- factor(sub(".*_(Zone_[0-9]+)$", "\\1", zone_df$Sample), levels = paste0("Zone_", 1:8))
      zone_df$Sample_ID <- order_samples(sub("_Zone_[0-9]+$", "", zone_df$Sample))
      zone_df <- zone_df %>%
        dplyr::arrange(Sample_ID) %>%
        mutate(Group = case_when(
          grepl("^M", Sample_ID) ~ "M",
          grepl("^P", Sample_ID) ~ "P",
          TRUE ~ NA_character_
        )) %>%
        filter(!is.na(Group)) %>%
        group_by(Sample_ID) %>%
        mutate(ZT_diff = ZT - ZT[Zone == "Zone_1"]) %>%
        ungroup()
      
      bulk_lookup_df <- data.frame(
        Sample_ID = order_samples(names(bulk_zt_lookup_JTK)),
        Bulk_ZT   = as.numeric(bulk_zt_lookup_JTK)
      )
      
      # ---- PLOT: Vorhergesagte ZT (JTK) über die Zonen, alle Samples zusammen und je Sample mit gestrichelter Bulk-ZT ----
      # p_all <- ggplot(zone_df, aes(Zone, ZT, group = Sample_ID, colour = Sample_ID)) +
      #   geom_line(linewidth = 1) + geom_point(size = 2) +
      #   scale_y_continuous(limits = c(0,24), breaks = 0:24) +
      #   labs(title = "Alle Samples", x = "Zone", y = "ZT") + theme_bw()
      # p_individual <- ggplot(zone_df, aes(Zone, ZT, group = Sample_ID)) +
      #   geom_line(colour = "steelblue", linewidth = 0.8) + geom_point(size = 2) +
      #   geom_hline(data = bulk_lookup_df, aes(yintercept = Bulk_ZT),
      #              linetype = "dashed", colour = "grey40", linewidth = 0.4) +
      #   facet_wrap(~Sample_ID, ncol = 4) +
      #   scale_x_discrete(labels = 1:8) +
      #   scale_y_continuous(limits = c(0,24), breaks = seq(0,24,2)) +
      #   labs(x = "Zone", y = "ZT") + theme_bw()
      # print(p_all / p_individual + plot_layout(heights = c(1, 3)))
      
      # ---- PLOT: Heatmap Sample x Zone mit der vorhergesagten ZT (JTK) ----
      # print(ggplot(zone_df, aes(Zone, factor(Sample_ID, levels = rev(levels(droplevels(Sample_ID)))), fill = ZT)) +
      #         geom_tile(colour = "white", linewidth = 0.5) +
      #         geom_text(aes(label = round(ZT, 1)), size = 3) +
      #         scale_fill_gradientn(colours = zt_colours_continuous, values = scales::rescale(0:24),
      #                              limits = c(0, 24), breaks = seq(0, 24, 2), oob = scales::squish, name = "ZT") +
      #         labs(title = "Vorhergesagte ZT (JTK) pro Sample und Zone", x = "Zone", y = "Sample") +
      #         theme_minimal() + theme(panel.grid = element_blank()))
      
      # ---- PLOT: Heatmap Sample x Zone mit Delta ZT relativ zu Zone 1 (JTK) ----
      # print(ggplot(zone_df, aes(Zone, factor(Sample_ID, levels = rev(levels(droplevels(Sample_ID)))), fill = ZT_diff)) +
      #         geom_tile(colour = "white", linewidth = 0.5) +
      #         geom_text(aes(label = round(ZT_diff, 1)), size = 3) +
      #         scale_fill_gradient2(low = "#3D6FA8", mid = "white", high = "#D9636B", midpoint = 0,
      #                              name = "Delta ZT\nvs Zone 1") +
      #         labs(title = "Delta ZT (JTK) relativ zu Zone 1, pro Sample", x = "Zone", y = "Sample") +
      #         theme_minimal() + theme(panel.grid = element_blank()))
      
      # ---- PLOT: Räumliche Plots aller 16 Proben, jede Zone in der Farbe ihrer vorhergesagten ZT (JTK); danach Liver Zones und Predicted ZT nebeneinander, verteilt auf 2 Seiten ----
      # plots <- list()
      # plots_pair <- list()
      # for(i in daten){
      #   obj <- readRDS(file.path(path_to_visium, paste0(i, "_Visium"), paste0(i, "_Visium_zoniert.rds")))
      #   zone_time <- round(zone_df$ZT[zone_df$Sample_ID == i])
      #   obj$PredictedTime <- as.numeric(zone_time[obj$med_liver_zone])
      #   p_time <- SpatialFeaturePlot(obj, features = "PredictedTime", pt.size.factor = 2) +
      #     zt_colour_scale() + ggtitle(i)
      #   plots[[i]] <- p_time
      #   p_zones <- SpatialFeaturePlot(subset(obj, !is.na(liver_Zone)), features = "med_liver_zone",
      #                                 pt.size.factor = 2) + ggtitle(paste("Liver Zones -", i))
      #   plots_pair[[i]] <- p_zones | p_time
      # }
      # print(wrap_plots(plots, ncol = 4, nrow = 4, guides = "collect") & theme(legend.position = "bottom"))
      # half <- ceiling(length(plots_pair) / 2)
      # print(wrap_plots(plots_pair[1:half], ncol = 2, guides = "collect"))
      # print(wrap_plots(plots_pair[(half + 1):length(plots_pair)], ncol = 2, guides = "collect"))
      
      # ---- PLOT: Delta ZT relativ zu Zone 1 (JTK) als Linien, alle Samples / nur M / nur P, einmal farbig nach Gruppe und einmal nach Sample ----
      # y_range <- range(zone_df$ZT_diff, na.rm = TRUE)
      # make_plot <- function(df, title, colour_by) {
      #   ggplot(df, aes(Zone, ZT_diff, group = Sample_ID, colour = .data[[colour_by]])) +
      #     geom_line(linewidth = 1) + geom_point(size = 2) +
      #     geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
      #     scale_y_continuous(limits = y_range, breaks = pretty(y_range), labels = hour_labels) +
      #     labs(title = title, x = "Zone", y = "\u0394 ZT relativ zu Zone 1") + theme_bw()
      # }
      # for (colour_by in c("Group", "Sample_ID")) {
      #   print(make_plot(zone_df, "Alle Samples (M & P)", colour_by) /
      #           make_plot(filter(zone_df, Group == "M"), "Nur M-Samples", colour_by) /
      #           make_plot(filter(zone_df, Group == "P"), "Nur P-Samples", colour_by))
      # }
      
      # ---- PLOT: Boxplots der Bin-ZT (JTK) pro Zone, M vs P gepoolt und pro Sample (braucht bin_zone_lookup_JTK.rds aus einem früheren bin-Durchlauf) ----
      # bin_plot_df <- readRDS(file.path(path_to_pred, "bin_zone_lookup_JTK.rds")) %>%
      #   mutate(Sample_ID = order_samples(Sample), Group = ifelse(grepl("^M", Sample), "M", "P"))
      # print(ggplot(bin_plot_df, aes(Zone, ZT, colour = Group, fill = Group)) +
      #         geom_boxplot(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.4) +
      #         scale_colour_manual(values = group_colours) + scale_fill_manual(values = group_colours) +
      #         labs(title = "Bin-ZT (JTK) pro Zone, M vs P (alle Bins gepoolt)", x = "Zone", y = "ZT (Bin-Ebene)") +
      #         theme_bw())
      # print(ggplot(bin_plot_df, aes(Zone, ZT)) +
      #         geom_boxplot(fill = "steelblue", alpha = 0.5) +
      #         facet_wrap(~Sample_ID, ncol = 4) +
      #         labs(title = "Bin-ZT (JTK) pro Zone, pro Sample", x = "Zone", y = "ZT (Bin-Ebene)") +
      #         theme_bw())
      
      
      ##### Statistik: M vs P pro Zone (t-Test & Wilcoxon, BH-korrigiert) #####
      
      summary_mp <- zone_df %>%
        group_by(Zone, Group) %>%
        summarise(mean_ZT = mean(ZT, na.rm = TRUE), sd_ZT = sd(ZT, na.rm = TRUE),
                  n = n(), .groups = "drop")
      
      test_mp <- zone_df %>%
        group_by(Zone) %>%
        summarise(
          p_ttest  = tryCatch(t.test(ZT ~ Group, data = cur_data())$p.value, error = function(e) NA_real_),
          p_wilcox = tryCatch(wilcox.test(ZT ~ Group, data = cur_data(), exact = FALSE)$p.value, error = function(e) NA_real_),
          .groups = "drop"
        ) %>%
        mutate(
          p_ttest_adj  = p.adjust(p_ttest, method = "BH"),
          p_wilcox_adj = p.adjust(p_wilcox, method = "BH"),
          sig_ttest    = sig_stars(p_ttest_adj),
          sig_wilcox   = sig_stars(p_wilcox_adj)
        )
      print_full(test_mp, "M vs P pro Zone (JTK)")
      
      # ---- PLOT: Mittelwert ± SD der ZT pro Zone, M vs P, mit Signifikanzsternen (t-Test und Wilcoxon) ----
      # label_df <- summary_mp %>%
      #   group_by(Zone) %>%
      #   summarise(y_top = max(mean_ZT + sd_ZT, na.rm = TRUE), .groups = "drop") %>%
      #   left_join(test_mp, by = "Zone")
      # print(ggplot(summary_mp, aes(Zone, mean_ZT, colour = Group, group = Group)) +
      #         geom_line(linewidth = 1, position = position_dodge(width = 0.3)) +
      #         geom_point(size = 2.5, position = position_dodge(width = 0.3)) +
      #         geom_errorbar(aes(ymin = mean_ZT - sd_ZT, ymax = mean_ZT + sd_ZT), width = 0.2,
      #                       position = position_dodge(width = 0.3)) +
      #         geom_text(data = label_df, aes(Zone, y_top + 0.8, label = paste0("t: ", sig_ttest)), inherit.aes = FALSE) +
      #         geom_text(data = label_df, aes(Zone, y_top + 1.6, label = paste0("W: ", sig_wilcox)), inherit.aes = FALSE) +
      #         scale_colour_manual(values = group_colours) +
      #         labs(title = "ZT (JTK): Mittelwert ± SD, M vs P pro Zone", x = "Zone", y = "ZT") +
      #         theme_bw())
      
      # ---- PLOT: Boxplot der ZT pro Zone, M vs P, mit Signifikanzsternen (t-Test und Wilcoxon) ----
      # label_df <- zone_df %>%
      #   group_by(Zone) %>%
      #   summarise(y_top = max(ZT, na.rm = TRUE), .groups = "drop") %>%
      #   left_join(test_mp, by = "Zone")
      # print(ggplot(zone_df, aes(Zone, ZT, colour = Group, fill = Group)) +
      #         geom_boxplot(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.4) +
      #         geom_text(data = label_df, aes(Zone, y_top + 0.8, label = paste0("t: ", sig_ttest)), inherit.aes = FALSE) +
      #         geom_text(data = label_df, aes(Zone, y_top + 1.6, label = paste0("W: ", sig_wilcox)), inherit.aes = FALSE) +
      #         scale_colour_manual(values = group_colours) + scale_fill_manual(values = group_colours) +
      #         labs(title = "ZT (JTK): Boxplot, M vs P pro Zone", x = "Zone", y = "ZT") +
      #         theme_bw())
      
      
      ##### Statistik: Delta ZT zu Zone 1, getestet gegen 0 (innerhalb M bzw. P) #####
      
      summary_within <- zone_df %>%
        group_by(Zone, Group) %>%
        summarise(mean_diff = mean(ZT_diff, na.rm = TRUE), sd_diff = sd(ZT_diff, na.rm = TRUE),
                  n = n(), .groups = "drop")
      
      test_within <- zone_df %>%
        filter(Zone != "Zone_1") %>%
        group_by(Zone, Group) %>%
        summarise(
          p_ttest  = tryCatch(t.test(ZT_diff, mu = 0)$p.value, error = function(e) NA_real_),
          p_wilcox = tryCatch(wilcox.test(ZT_diff, mu = 0, exact = FALSE)$p.value, error = function(e) NA_real_),
          .groups = "drop"
        ) %>%
        group_by(Group) %>%
        mutate(
          p_ttest_adj  = p.adjust(p_ttest, method = "BH"),
          p_wilcox_adj = p.adjust(p_wilcox, method = "BH"),
          sig_ttest    = sig_stars(p_ttest_adj),
          sig_wilcox   = sig_stars(p_wilcox_adj)
        ) %>%
        ungroup()
      print_full(test_within, "Delta zu Zone 1 vs 0 (JTK)")
      
      # ---- PLOT: Mittelwert ± SD des Delta ZT zu Zone 1 pro Zone, M und P getrennt, mit Signifikanzsternen ----
      # label_df <- summary_within %>%
      #   mutate(y_top = mean_diff + sd_diff) %>%
      #   left_join(test_within, by = c("Zone", "Group"))
      # print(ggplot(summary_within, aes(Zone, mean_diff, colour = Group, group = Group)) +
      #         geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
      #         geom_line(linewidth = 1, position = position_dodge(width = 0.3)) +
      #         geom_point(size = 2.5, position = position_dodge(width = 0.3)) +
      #         geom_errorbar(aes(ymin = mean_diff - sd_diff, ymax = mean_diff + sd_diff), width = 0.2,
      #                       position = position_dodge(width = 0.3)) +
      #         geom_text(data = label_df, aes(Zone, y_top + 0.6, label = paste0("t: ", sig_ttest), colour = Group),
      #                   position = position_dodge(width = 0.3), size = 3, show.legend = FALSE) +
      #         geom_text(data = label_df, aes(Zone, y_top + 1.2, label = paste0("W: ", sig_wilcox), colour = Group),
      #                   position = position_dodge(width = 0.3), size = 3, show.legend = FALSE) +
      #         scale_colour_manual(values = group_colours) +
      #         scale_y_continuous(labels = hour_labels) +
      #         labs(title = "Delta ZT (JTK) relativ zu Zone 1, getestet gegen 0", x = "Zone", y = "Delta ZT relativ zu Zone 1") +
      #         theme_bw())
      
      # ---- PLOT: Boxplot des Delta ZT zu Zone 1 pro Zone, M und P getrennt, mit Signifikanzsternen ----
      # label_df <- zone_df %>%
      #   group_by(Zone, Group) %>%
      #   summarise(y_top = max(ZT_diff, na.rm = TRUE), .groups = "drop") %>%
      #   left_join(test_within, by = c("Zone", "Group"))
      # print(ggplot(zone_df, aes(Zone, ZT_diff, colour = Group, fill = Group)) +
      #         geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
      #         geom_boxplot(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.4) +
      #         geom_text(data = label_df, aes(Zone, y_top + 0.6, label = paste0("t: ", sig_ttest), colour = Group),
      #                   position = position_dodge(width = 0.7), size = 3, show.legend = FALSE) +
      #         geom_text(data = label_df, aes(Zone, y_top + 1.4, label = paste0("W: ", sig_wilcox), colour = Group),
      #                   position = position_dodge(width = 0.7), size = 3, show.legend = FALSE) +
      #         scale_colour_manual(values = group_colours) + scale_fill_manual(values = group_colours) +
      #         scale_y_continuous(labels = hour_labels) +
      #         labs(title = "Delta ZT (JTK, Boxplot) relativ zu Zone 1, getestet gegen 0", x = "Zone", y = "Delta ZT relativ zu Zone 1") +
      #         theme_bw())
      
      
      ##### Statistik: paarweise Zonen-Vergleiche, gepaart über Sample_ID #####
      
      pairwise_paired_results <- lapply(c("M", "P"), function(g) {
        wide_df <- zone_df %>%
          filter(Group == g) %>%
          dplyr::select(Sample_ID, Zone, ZT) %>%
          pivot_wider(names_from = Zone, values_from = ZT)
        combs <- combn(levels(droplevels(zone_df$Zone[zone_df$Group == g])), 2, simplify = FALSE)
        bind_rows(lapply(combs, function(pair) {
          z1 <- wide_df[[pair[1]]]
          z2 <- wide_df[[pair[2]]]
          valid <- complete.cases(z1, z2)
          data.frame(
            Group = g, Zone1 = pair[1], Zone2 = pair[2], n_paired = sum(valid),
            p_ttest_paired  = tryCatch(t.test(z1[valid], z2[valid], paired = TRUE)$p.value, error = function(e) NA_real_),
            p_wilcox_paired = tryCatch(wilcox.test(z1[valid], z2[valid], paired = TRUE, exact = FALSE)$p.value, error = function(e) NA_real_)
          )
        }))
      }) %>%
        bind_rows() %>%
        group_by(Group) %>%
        mutate(
          p_ttest_adj  = p.adjust(p_ttest_paired, method = "BH"),
          p_wilcox_adj = p.adjust(p_wilcox_paired, method = "BH"),
          sig_ttest    = sig_stars(p_ttest_adj),
          sig_wilcox   = sig_stars(p_wilcox_adj)
        ) %>%
        ungroup()
      print_full(pairwise_paired_results, "Paarweise Zonen-Vergleiche (JTK)")
      
      ##### Statistik: Friedman-Test (unterscheiden sich die Zonen insgesamt?) #####
      
      friedman_results <- lapply(c("M", "P"), function(g) {
        mat <- zone_df %>%
          filter(Group == g) %>%
          dplyr::select(Sample_ID, Zone, ZT) %>%
          pivot_wider(names_from = Zone, values_from = ZT) %>%
          filter(complete.cases(.))
        mat <- as.matrix(mat[, -1])
        fried <- friedman.test(mat)
        data.frame(Group = g, statistic = fried$statistic, p_value = fried$p.value, n_complete = nrow(mat))
      }) %>%
        bind_rows()
      print_full(friedman_results, "Friedman-Test (JTK)")
      
      ##### Statistik: Bootstrap paarweiser Zonen-Differenzen (gepaart über Sample_ID) #####
      
      library(boot)
      set.seed(123)
      
      bootstrap_pairwise <- lapply(c("M", "P"), function(g) {
        wide_df <- zone_df %>%
          filter(Group == g) %>%
          dplyr::select(Sample_ID, Zone, ZT) %>%
          pivot_wider(names_from = Zone, values_from = ZT) %>%
          filter(complete.cases(.))
        combs <- combn(levels(droplevels(zone_df$Zone)), 2, simplify = FALSE)
        bind_rows(lapply(combs, function(pair) {
          diff_stat <- function(data, indices) {
            d <- data[indices, ]
            mean(d[[pair[1]]] - d[[pair[2]]], na.rm = TRUE)
          }
          boot_out <- boot(data = wide_df, statistic = diff_stat, R = 10000)
          ci <- boot.ci(boot_out, type = "perc")
          boot_diffs <- boot_out$t
          observed_diff <- boot_out$t0
          data.frame(
            Group = g, Zone1 = pair[1], Zone2 = pair[2],
            observed_diff = observed_diff,
            ci_lower = ci$percent[4],
            ci_upper = ci$percent[5],
            p_boot = mean(abs(boot_diffs - mean(boot_diffs)) >= abs(observed_diff - mean(boot_diffs)))
          )
        }))
      }) %>%
        bind_rows() %>%
        group_by(Group) %>%
        mutate(p_boot_adj = p.adjust(p_boot, method = "BH"), sig_boot = sig_stars(p_boot_adj)) %>%
        ungroup()
      print_full(bootstrap_pairwise, "Bootstrap paarweise Zonen-Differenzen (JTK)")
      
      # ---- PLOT: Bootstrap-Konfidenzintervalle aller paarweisen Zonen-Differenzen, M und P getrennt ----
      # print(bootstrap_pairwise %>%
      #         mutate(Comparison = paste(Zone1, "vs", Zone2)) %>%
      #         ggplot(aes(observed_diff, Comparison, colour = Group)) +
      #         geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
      #         geom_pointrange(aes(xmin = ci_lower, xmax = ci_upper)) +
      #         facet_wrap(~ Group) +
      #         scale_colour_manual(values = group_colours) +
      #         labs(title = "Bootstrap-Konfidenzintervalle (JTK): paarweise Zonen-Differenzen",
      #              x = "Differenz (Zone1 - Zone2)", y = NULL) +
      #         theme_bw() + theme(legend.position = "none"))
      
      # ---- PLOT: Heatmap Zone x Zone der Bootstrap-Differenzen mit Signifikanzsternen, M und P getrennt ----
      # heatmap_df <- bind_rows(
      #   bootstrap_pairwise,
      #   bootstrap_pairwise %>% rename(Zone1 = Zone2, Zone2 = Zone1) %>% mutate(observed_diff = -observed_diff),
      #   expand.grid(Group = c("M", "P"), Zone1 = paste0("Zone_", 1:8), stringsAsFactors = FALSE) %>%
      #     mutate(Zone2 = Zone1, observed_diff = 0, sig_boot = "")
      # ) %>%
      #   mutate(Zone1 = factor(Zone1, levels = paste0("Zone_", 1:8)),
      #          Zone2 = factor(Zone2, levels = paste0("Zone_", 1:8)))
      # print(ggplot(heatmap_df, aes(Zone2, Zone1, fill = observed_diff)) +
      #         geom_tile(colour = "white") +
      #         geom_text(aes(label = sig_boot), size = 3.5) +
      #         facet_wrap(~ Group) +
      #         scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0,
      #                              name = "Differenz\n(Zone1 - Zone2)") +
      #         labs(title = "Bootstrap-Differenzen (JTK) zwischen allen Zonen-Paaren", x = "Zone 2", y = "Zone 1") +
      #         theme_bw())
      
      ##### Statistik: Bootstrap-Mittelwerte pro Zone mit 95%-CI #####
      
      bootstrap_zone_means <- lapply(c("M", "P"), function(g) {
        sub_df <- zone_df %>% filter(Group == g)
        bind_rows(lapply(levels(droplevels(sub_df$Zone)), function(z) {
          boot_out <- boot(data = sub_df$ZT[sub_df$Zone == z],
                           statistic = function(data, indices) mean(data[indices], na.rm = TRUE),
                           R = 10000)
          ci <- boot.ci(boot_out, type = "perc")
          data.frame(Group = g, Zone = z, mean_ZT = boot_out$t0,
                     ci_lower = ci$percent[4], ci_upper = ci$percent[5])
        }))
      }) %>%
        bind_rows() %>%
        mutate(Zone = factor(Zone, levels = paste0("Zone_", 1:8)))
      print_full(bootstrap_zone_means, "Bootstrap-Mittelwerte pro Zone (JTK)")
      
      # ---- PLOT: ZT pro Zone als Bootstrap-Mittelwert mit 95%-CI, M vs P ----
      # print(ggplot(bootstrap_zone_means, aes(Zone, mean_ZT, colour = Group, group = Group)) +
      #         geom_line(linewidth = 1, position = position_dodge(width = 0.3)) +
      #         geom_pointrange(aes(ymin = ci_lower, ymax = ci_upper), position = position_dodge(width = 0.3)) +
      #         scale_colour_manual(values = group_colours) +
      #         labs(title = "ZT (JTK) pro Zone mit Bootstrap-95%-Konfidenzintervallen", x = "Zone",
      #              y = "ZT (Mittelwert ± Bootstrap-CI)") +
      #         theme_bw())
      
      ##### Statistik: Abweichung Zone - Bulk (zirkulär), getestet gegen 0 #####
      
      bulk_df <- zone_df
      bulk_df$Bulk_ZT <- bulk_zt_lookup_JTK[as.character(bulk_df$Sample_ID)]
      missing_bulk <- unique(as.character(bulk_df$Sample_ID[is.na(bulk_df$Bulk_ZT)]))
      if(length(missing_bulk) > 0){
        message("Kein Bulk-Wert (JTK) gefunden für: ", paste(missing_bulk, collapse = ", "))
      }
      bulk_df <- bulk_df %>%
        filter(!is.na(Bulk_ZT)) %>%
        mutate(ZT_diff_bulk = ((ZT - Bulk_ZT + 12) %% 24) - 12)  # kürzester Weg auf der 24h-Uhr
      
      test_vs_bulk <- bulk_df %>%
        group_by(Zone, Group) %>%
        summarise(
          p_ttest  = tryCatch(t.test(ZT_diff_bulk, mu = 0)$p.value, error = function(e) NA_real_),
          p_wilcox = tryCatch(wilcox.test(ZT_diff_bulk, mu = 0, exact = FALSE)$p.value, error = function(e) NA_real_),
          .groups = "drop"
        ) %>%
        group_by(Group) %>%
        mutate(
          p_ttest_adj  = p.adjust(p_ttest, method = "BH"),
          p_wilcox_adj = p.adjust(p_wilcox, method = "BH"),
          sig_ttest    = sig_stars(p_ttest_adj),
          sig_wilcox   = sig_stars(p_wilcox_adj)
        ) %>%
        ungroup()
      print_full(test_vs_bulk, "Zone - Bulk vs 0 (JTK)")
      
      overall_vs_bulk <- bulk_df %>%
        group_by(Zone) %>%
        summarise(
          mean_diff = mean(ZT_diff_bulk, na.rm = TRUE),
          sd_diff   = sd(ZT_diff_bulk, na.rm = TRUE),
          n         = n(),
          p_ttest   = tryCatch(t.test(ZT_diff_bulk, mu = 0)$p.value, error = function(e) NA_real_),
          .groups = "drop"
        ) %>%
        mutate(p_ttest_adj = p.adjust(p_ttest, method = "BH"))
      print_full(overall_vs_bulk, "Mittlere Abweichung Zone - Bulk, alle Samples (JTK)")
      
      # ---- PLOT: Heatmap Sample x Zone der Abweichung Zone - Bulk (JTK) ----
      # print(ggplot(bulk_df, aes(Zone, factor(Sample_ID, levels = rev(levels(droplevels(Sample_ID)))), fill = ZT_diff_bulk)) +
      #         geom_tile(colour = "white", linewidth = 0.5) +
      #         geom_text(aes(label = round(ZT_diff_bulk, 1)), size = 3) +
      #         scale_fill_gradient2(low = "#3D6FA8", mid = "white", high = "#D9636B", midpoint = 0,
      #                              name = "Delta ZT\n(Zone - Bulk)") +
      #         labs(title = "Delta ZT (JTK, Zone - Bulk) pro Sample und Zone", x = "Zone", y = "Sample") +
      #         theme_minimal() + theme(panel.grid = element_blank()))
      
      # ---- PLOT: Boxplot der Abweichung Zone - Bulk pro Zone, M und P getrennt, mit Signifikanzsternen ----
      # label_df <- bulk_df %>%
      #   group_by(Zone, Group) %>%
      #   summarise(y_top = max(ZT_diff_bulk, na.rm = TRUE), .groups = "drop") %>%
      #   left_join(test_vs_bulk, by = c("Zone", "Group"))
      # print(ggplot(bulk_df, aes(Zone, ZT_diff_bulk, colour = Group, fill = Group)) +
      #         geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
      #         geom_boxplot(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.4) +
      #         geom_text(data = label_df, aes(Zone, y_top + 0.6, label = paste0("t: ", sig_ttest), colour = Group),
      #                   position = position_dodge(width = 0.7), size = 3, show.legend = FALSE) +
      #         geom_text(data = label_df, aes(Zone, y_top + 1.4, label = paste0("W: ", sig_wilcox), colour = Group),
      #                   position = position_dodge(width = 0.7), size = 3, show.legend = FALSE) +
      #         scale_colour_manual(values = group_colours) + scale_fill_manual(values = group_colours) +
      #         scale_y_continuous(labels = hour_labels) +
      #         labs(title = "Delta ZT (JTK, Zone - Bulk), getestet gegen 0", x = "Zone", y = "Delta ZT (Zone - Bulk)") +
      #         theme_bw())
      
      # ---- PLOT: Balken der mittleren Abweichung jeder Zone vom Bulk (alle Samples gepoolt, ± SD) ----
      # print(ggplot(overall_vs_bulk, aes(Zone, mean_diff)) +
      #         geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
      #         geom_col(fill = "steelblue", alpha = 0.6) +
      #         geom_errorbar(aes(ymin = mean_diff - sd_diff, ymax = mean_diff + sd_diff), width = 0.2) +
      #         scale_y_continuous(labels = hour_labels) +
      #         labs(title = "Mittlere Abweichung (JTK) jeder Zone vom Bulk-Ergebnis", x = "Zone",
      #              y = "Mittleres Delta ZT (Zone - Bulk)") +
      #         theme_bw())
      
      # Alle Statistik-Ergebnisse speichern -> pred/
      saveRDS(
        list(summary_mp = summary_mp, test_mp = test_mp,
             summary_within = summary_within, test_within = test_within,
             pairwise_paired = pairwise_paired_results, friedman = friedman_results,
             bootstrap_pairwise = bootstrap_pairwise, bootstrap_zone_means = bootstrap_zone_means,
             test_vs_bulk = test_vs_bulk, overall_vs_bulk = overall_vs_bulk),
        file = file.path(path_to_pred, "zonen_statistik_JTK.rds")
      )
    }
    
    ##### Bins (JTK) #####
    if(bin == TRUE & type == "bin"){
      
      # Bin -> Zone Mapping aus den zonierten Seurat-Objekten
      zone_map <- bind_rows(lapply(daten, function(i) {
        obj <- readRDS(file.path(path_to_visium, paste0(i, "_Visium"), paste0(i, "_Visium_zoniert.rds")))
        data.frame(Bin = paste0(i, "_", colnames(obj)), Zone = obj$med_liver_zone,
                   stringsAsFactors = FALSE)
      })) %>%
        filter(!is.na(Zone), Zone > 0, Zone <= 8)   # unzonierte Spots raus
      
      bin_zone_df <- data.frame(
        Bin = colnames(pred_data),
        ZT = result_JTK$predicted,
        stringsAsFactors = FALSE
      ) %>%
        mutate(Sample = sub("^(M[0-9]+|P[0-9]+).*", "\\1", Bin)) %>%
        left_join(zone_map, by = "Bin") %>%
        filter(!is.na(Zone)) %>%
        mutate(Zone = factor(paste0("Zone_", Zone), levels = paste0("Zone_", 1:8)))
      
      # Bin-Zone-Lookup (ZT gerundet, 0-23) speichern -> pred/
      bin_zone_lookup_JTK <- bin_zone_df %>% mutate(ZT = round(ZT) %% 24)
      saveRDS(bin_zone_lookup_JTK, file = file.path(path_to_pred, "bin_zone_lookup_JTK.rds"))
      
      # ---- PLOT: Pro Sample ein Balkendiagramm der Bin-ZT (JTK) mit 8 Panels (eine pro Zone), gestrichelt die Bulk-ZT, gepunktet die Zonen-ZT ----
      # if(!exists("bulk_zt_lookup_JTK", inherits = FALSE)){
      #   bulk_zt_lookup_JTK <- readRDS(file.path(path_to_pred, "bulk_zt_lookup_JTK.rds"))
      # }
      # if(!exists("zonen_zt_lookup_JTK", inherits = FALSE)){
      #   zonen_zt_lookup_JTK <- readRDS(file.path(path_to_pred, "zonen_zt_lookup_JTK.rds"))
      # }
      # zonen_lookup_df <- data.frame(ID = names(zonen_zt_lookup_JTK), Zonen_ZT = zonen_zt_lookup_JTK) %>%
      #   mutate(Zone = factor(sub(".*_(Zone_[0-9]+)$", "\\1", ID), levels = paste0("Zone_", 1:8)),
      #          Sample = sub("_Zone_[0-9]+$", "", ID),
      #          x_pos = round(Zonen_ZT) %% 24 + 1)
      # for (s in intersect(sample_level_order, bin_zone_lookup_JTK$Sample)) {
      #   p <- ggplot(filter(bin_zone_lookup_JTK, Sample == s), aes(factor(ZT), fill = factor(ZT))) +
      #     geom_bar(colour = "black") +
      #     facet_wrap(~Zone, ncol = 4, drop = FALSE) +
      #     geom_vline(xintercept = round(bulk_zt_lookup_JTK[s]) %% 24 + 1, linetype = "dashed", colour = "grey40") +
      #     geom_vline(data = filter(zonen_lookup_df, Sample == s), aes(xintercept = x_pos), linetype = "dotted") +
      #     scale_fill_manual(values = zt_colours_24, name = "ZT") +
      #     scale_x_discrete(limits = as.character(0:23), breaks = seq(0, 23, by = 2)) +
      #     labs(title = paste("Sample (JTK):", s), x = "Zeitgeberzeit (ZT)", y = "Anzahl Bins") +
      #     theme_bw() + theme(legend.position = "bottom")
      #   print(p)
      # }
      
      # ---- PLOT: Balkendiagramm der Bin-ZT (JTK) pro Sample (alle Zonen zusammen) ----
      # print(ggplot(mutate(bin_zone_lookup_JTK, Sample = order_samples(Sample)), aes(factor(ZT), fill = factor(ZT))) +
      #         geom_bar(colour = "black") +
      #         facet_wrap(~Sample, ncol = 4) +
      #         scale_fill_manual(values = zt_colours_24, name = "ZT") +
      #         scale_x_discrete(limits = as.character(0:23), breaks = seq(0, 23, by = 2)) +
      #         labs(x = "Zeitgeberzeit (ZT)", y = "Anzahl Bins") +
      #         theme_bw() + theme(legend.position = "bottom"))
      
      # ---- PLOT: Räumliche Plots aller 16 Proben mit der vorhergesagten ZT (JTK) pro Bin; danach Liver Zones und Bin-ZT nebeneinander, verteilt auf 2 Seiten ----
      # plots <- list()
      # plots_pair <- list()
      # for(i in daten){
      #   obj <- readRDS(file.path(path_to_visium, paste0(i, "_Visium"), paste0(i, "_Visium_zoniert.rds")))
      #   colnames(obj) <- paste0(i, "_", colnames(obj))
      #   obj$PredictedTime <- result_JTK$predicted[match(colnames(obj), colnames(pred_data))]
      #   p_time <- SpatialFeaturePlot(obj, features = "PredictedTime", pt.size.factor = 2) +
      #     zt_colour_scale() + ggtitle(i)
      #   plots[[i]] <- p_time
      #   p_zones <- SpatialFeaturePlot(subset(obj, !is.na(liver_Zone)), features = "med_liver_zone",
      #                                 pt.size.factor = 2) + ggtitle(paste("Liver Zones -", i))
      #   plots_pair[[i]] <- p_zones | p_time
      # }
      # print(wrap_plots(plots, ncol = 4, nrow = 4, guides = "collect") & theme(legend.position = "bottom"))
      # half <- ceiling(length(plots_pair) / 2)
      # print(wrap_plots(plots_pair[1:half], ncol = 2, guides = "collect"))
      # print(wrap_plots(plots_pair[(half + 1):length(plots_pair)], ncol = 2, guides = "collect"))
      
      
      ##### Statistik: Bootstrap auf Bin-Ebene, pro Sample #####
      # Innerhalb eines Samples wird auf den einzelnen BINS resampled. Bins
      # verschiedener Zonen sind nicht gepaart, daher UNGEPAARTES Bootstrap:
      # Zone und Zone_1 werden unabhaengig aus ihren eigenen Bins gezogen.
      
      set.seed(123)
      R_boot <- 2000          # ACHTUNG: rechenintensiv, bei Bedarf reduzieren
      min_bins_per_zone <- 5  # Zonen mit weniger Bins werden ausgeschlossen
      
      samples <- sample_level_order[sample_level_order %in% unique(bin_zone_df$Sample)]
      
      # Delta ZT relativ zu Zone 1 pro Sample und Zone
      bin_delta_zone1 <- lapply(samples, function(s) {
        sub_df <- bin_zone_df %>% filter(Sample == s)
        zone_counts <- table(droplevels(sub_df$Zone))
        zone_levels <- names(zone_counts)[zone_counts >= min_bins_per_zone]
        if (!("Zone_1" %in% zone_levels) || length(zone_levels) < 2) return(NULL)
        
        vals_zone1 <- sub_df$ZT[sub_df$Zone == "Zone_1"]
        
        bind_rows(lapply(setdiff(zone_levels, "Zone_1"), function(z) {
          vals_z <- sub_df$ZT[sub_df$Zone == z]
          boot_diffs <- replicate(R_boot, {
            mean(sample(vals_z, length(vals_z), replace = TRUE)) -
              mean(sample(vals_zone1, length(vals_zone1), replace = TRUE))
          })
          ci <- quantile(boot_diffs, probs = c(0.025, 0.975))
          data.frame(
            Sample = s, Zone = z,
            n_zone = length(vals_z), n_zone1 = length(vals_zone1),
            ZT_diff = mean(vals_z) - mean(vals_zone1),
            ci_lower = as.numeric(ci[1]), ci_upper = as.numeric(ci[2])
          )
        }))
      }) %>%
        bind_rows() %>%
        mutate(
          Sample = order_samples(Sample),
          Zone = factor(Zone, levels = paste0("Zone_", 2:8)),
          Group = ifelse(grepl("^M", Sample), "M", "P")
        )
      print_full(bin_delta_zone1, "Bin-Ebene: Delta zu Zone 1 (JTK)")
      
      # Delta-Werte gepoolt, pro Zone innerhalb M bzw. P gegen 0 getestet
      summary_bin <- bin_delta_zone1 %>%
        group_by(Zone, Group) %>%
        summarise(mean_diff = mean(ZT_diff, na.rm = TRUE), sd_diff = sd(ZT_diff, na.rm = TRUE),
                  n = n(), .groups = "drop")
      
      test_bin <- bin_delta_zone1 %>%
        group_by(Zone, Group) %>%
        summarise(
          p_ttest  = tryCatch(t.test(ZT_diff, mu = 0)$p.value, error = function(e) NA_real_),
          p_wilcox = tryCatch(wilcox.test(ZT_diff, mu = 0, exact = FALSE)$p.value, error = function(e) NA_real_),
          .groups = "drop"
        ) %>%
        group_by(Group) %>%
        mutate(
          p_ttest_adj  = p.adjust(p_ttest, method = "BH"),
          p_wilcox_adj = p.adjust(p_wilcox, method = "BH"),
          sig_ttest    = sig_stars(p_ttest_adj),
          sig_wilcox   = sig_stars(p_wilcox_adj)
        ) %>%
        ungroup()
      print_full(test_bin, "Bin-Ebene: Delta zu Zone 1 vs 0 (JTK)")
      
      # ---- PLOT: Mittelwert ± SD des Bin-Delta zu Zone 1 pro Zone, M und P getrennt, mit Signifikanzsternen ----
      # label_df <- summary_bin %>%
      #   mutate(y_top = mean_diff + sd_diff) %>%
      #   left_join(test_bin, by = c("Zone", "Group"))
      # print(ggplot(summary_bin, aes(Zone, mean_diff, colour = Group, group = Group)) +
      #         geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
      #         geom_line(linewidth = 1, position = position_dodge(width = 0.3)) +
      #         geom_point(size = 2.5, position = position_dodge(width = 0.3)) +
      #         geom_errorbar(aes(ymin = mean_diff - sd_diff, ymax = mean_diff + sd_diff), width = 0.2,
      #                       position = position_dodge(width = 0.3)) +
      #         geom_text(data = label_df, aes(Zone, y_top + 0.6, label = paste0("t: ", sig_ttest), colour = Group),
      #                   position = position_dodge(width = 0.3), size = 3, show.legend = FALSE) +
      #         geom_text(data = label_df, aes(Zone, y_top + 1.2, label = paste0("W: ", sig_wilcox), colour = Group),
      #                   position = position_dodge(width = 0.3), size = 3, show.legend = FALSE) +
      #         scale_colour_manual(values = group_colours) +
      #         scale_y_continuous(labels = hour_labels) +
      #         labs(title = "Delta ZT (JTK, Bin-Ebene) relativ zu Zone 1, getestet gegen 0",
      #              x = "Zone", y = "Delta ZT relativ zu Zone 1") +
      #         theme_bw())
      
      # Bootstrap-Mittelwert der Bin-ZT pro Zone mit 95%-CI, pro Sample
      bin_zone_means_boot <- lapply(samples, function(s) {
        sub_df <- bin_zone_df %>% filter(Sample == s)
        zone_counts <- table(droplevels(sub_df$Zone))
        zone_levels <- names(zone_counts)[zone_counts >= min_bins_per_zone]
        bind_rows(lapply(zone_levels, function(z) {
          vals <- sub_df$ZT[sub_df$Zone == z]
          boot_means <- replicate(R_boot, mean(sample(vals, length(vals), replace = TRUE)))
          ci <- quantile(boot_means, probs = c(0.025, 0.975))
          data.frame(Sample = s, Zone = z, mean_ZT = mean(vals),
                     ci_lower = as.numeric(ci[1]), ci_upper = as.numeric(ci[2]), n = length(vals))
        }))
      }) %>%
        bind_rows() %>%
        mutate(
          Sample = order_samples(Sample),
          Zone = factor(Zone, levels = paste0("Zone_", 1:8)),
          Group = ifelse(grepl("^M", Sample), "M", "P")
        )
      print_full(bin_zone_means_boot, "Bin-Ebene: Bootstrap-Mittelwerte pro Zone (JTK)")
      
      # ---- PLOT: Bin-ZT pro Zone als Bootstrap-Mittelwert mit 95%-CI, ein Panel pro Sample ----
      # print(ggplot(bin_zone_means_boot, aes(Zone, mean_ZT, colour = Group, group = Sample)) +
      #         geom_line(linewidth = 0.6, alpha = 0.7) +
      #         geom_pointrange(aes(ymin = ci_lower, ymax = ci_upper), size = 0.35) +
      #         facet_wrap(~ Sample, ncol = 4) +
      #         scale_colour_manual(values = group_colours) +
      #         labs(title = "Bin-ZT (JTK): Mittelwert pro Zone mit Bootstrap-95%-CI, pro Sample",
      #              x = "Zone", y = "ZT (Bin-Ebene)") +
      #         theme_bw() + theme(legend.position = "none"))
      
      # Statistik-Ergebnisse speichern -> pred/
      saveRDS(
        list(bin_delta_zone1 = bin_delta_zone1, summary_bin = summary_bin,
             test_bin = test_bin, bin_zone_means_boot = bin_zone_means_boot),
        file = file.path(path_to_pred, "bin_statistik_JTK.rds")
      )
    }
    
    
    #### predict circadian time LS ####
    
    missing <- setdiff(trained_LS$Genes, rownames(pred_data))
    if (length(missing) > 0) message("Fehlende Gene (LS): ", paste(missing, collapse = ", "))
    cat("GENES_LS:", trained_LS$Genes, "\n")
    
    result_LS <- test_tauFisher(chosen_genes=trained_LS$Genes,
                                test_df=pred_data,
                                test_time=pred_time,
                                pca_embedding=trained_LS$PCA_embedding,
                                trained_model=trained_LS$Model)
    
    calc_error(result_LS$truth, result_LS$predicted)
    
    #### Predictions speichern (LS) -> pred/ ####
    pred_result_df_LS <- data.frame(
      ID = colnames(pred_data),
      ZT_predicted = result_LS$predicted,
      ZT_truth = result_LS$truth,
      method = "LS",
      type = type,
      stringsAsFactors = FALSE
    )
    saveRDS(pred_result_df_LS, file = file.path(path_to_pred, paste0(type, "_prediction_LS.rds")))
    
    
    #### predict circadian time LS (gesamte Tabelle) ####
    
    missing <- setdiff(trained_LS_full$Genes, rownames(pred_data))
    if (length(missing) > 0) message("Fehlende Gene (LS_full): ", paste(missing, collapse = ", "))
    cat("GENES_LS_full:", trained_LS_full$Genes, "\n")
    
    result_LS_full <- test_tauFisher(chosen_genes=trained_LS_full$Genes,
                                     test_df=pred_data,
                                     test_time=pred_time,
                                     pca_embedding=trained_LS_full$PCA_embedding,
                                     trained_model=trained_LS_full$Model)
    
    calc_error(result_LS_full$truth, result_LS_full$predicted)
    
    #### Predictions speichern (LS, gesamte Tabelle) -> pred/ ####
    pred_result_df_LS_full <- data.frame(
      ID = colnames(pred_data),
      ZT_predicted = result_LS_full$predicted,
      ZT_truth = result_LS_full$truth,
      method = "LS_full",
      type = type,
      stringsAsFactors = FALSE
    )
    saveRDS(pred_result_df_LS_full, file = file.path(path_to_pred, paste0(type, "_prediction_LS_full.rds")))
    
    ##### Plots LS #####
    if(bulk == TRUE & type == "bulk"){
      
      # ---- PLOT: Häufigkeit der vorhergesagten Bulk-ZT (LS) als Balkendiagramm, getrennt für M, P und M + P ----
      # freq_plot <- data.frame(Sample = colnames(pred_data), ZT = round(result_LS$predicted) %% 24) %>%
      #   mutate(Gruppe = ifelse(grepl("^M", Sample), "M", "P"))
      # freq_plot <- bind_rows(
      #   freq_plot %>% filter(Gruppe == "M") %>% count(ZT) %>% mutate(Panel = "M"),
      #   freq_plot %>% filter(Gruppe == "P") %>% count(ZT) %>% mutate(Panel = "P"),
      #   freq_plot %>% count(ZT) %>% mutate(Panel = "M + P")
      # ) %>%
      #   complete(Panel = c("M", "P", "M + P"), ZT = 0:23, fill = list(n = 0)) %>%
      #   mutate(Panel = factor(Panel, levels = c("M", "P", "M + P")), ZT = factor(ZT, levels = 0:23))
      # print(ggplot(freq_plot, aes(ZT, n, fill = ZT)) +
      #         geom_col(colour = "black") +
      #         facet_wrap(~Panel, ncol = 1) +
      #         scale_fill_manual(values = zt_colours_24, name = "ZT") +
      #         scale_x_discrete(breaks = as.character(seq(0, 23, by = 2))) +
      #         labs(title = "Häufigkeit der vorhergesagten ZT (LS)", x = "Zeitgeberzeit (ZT)", y = "Häufigkeit") +
      #         theme_bw() + theme(legend.position = "bottom"))
      
      # ---- PLOT: Pro Bulk-Sample eine Uhr mit der vorhergesagten ZT (LS) als Punkt ----
      # plot_df <- data.frame(Sample = order_samples(colnames(pred_data)), ZT = result_LS$predicted) %>%
      #   mutate(x = sin(2*pi*ZT/24), y = cos(2*pi*ZT/24))
      # circle <- data.frame(t = seq(0, 2*pi, length.out = 500)) %>% mutate(x = cos(t), y = sin(t))
      # hours <- data.frame(hour = seq(0, 21, 3)) %>%
      #   mutate(x = 1.15*sin(2*pi*hour/24), y = 1.15*cos(2*pi*hour/24),
      #          label = ifelse(hour == 0, "0/24", as.character(hour)))
      # print(ggplot() +
      #         geom_path(data = circle, aes(x, y), linewidth = 0.8) +
      #         geom_point(data = plot_df, aes(x, y), colour = "steelblue", size = 3) +
      #         geom_text(data = hours, aes(x, y, label = label), size = 3.5) +
      #         facet_wrap(~Sample, ncol = 4) +
      #         coord_equal(xlim = c(-1.3, 1.3), ylim = c(-1.3, 1.3)) +
      #         theme_void())
      
      # ---- PLOT: Trainings-Expression der Prädiktorgene (LS, Replikat 1) über die Zeit, alle zusammen und einzeln ----
      # plot_data_long <- as.data.frame(train_adj_clean[trained_LS$Genes, replicate_clean == 1, drop = FALSE]) %>%
      #   mutate(Gene = rownames(.)) %>%
      #   pivot_longer(cols = -Gene, names_to = "Time", values_to = "Expression") %>%
      #   mutate(Time = as.numeric(Time))
      # p_all <- ggplot(plot_data_long, aes(Time, Expression, color = Gene, group = Gene)) +
      #   geom_line(linewidth = 1) + geom_point(size = 2) +
      #   labs(title = "Alle Gene", x = "Zeit", y = "Expression") + theme_classic()
      # p_facet <- ggplot(plot_data_long, aes(Time, Expression, group = Gene)) +
      #   geom_line(linewidth = 1, color = "steelblue") + geom_point(size = 2, color = "steelblue") +
      #   facet_wrap(~Gene, scales = "free_y") +
      #   labs(x = "Zeit", y = "Expression") + theme_classic()
      # print(p_all / p_facet)
    }
    
    if(zonen == TRUE & type == "zonen"){
      
      # ---- PLOT: Vorhergesagte ZT (LS) über die Zonen, alle Samples zusammen und je Sample; danach Delta ZT zu Zone 1 (alle / M / P, farbig nach Gruppe und nach Sample) ----
      # zone_df_LS <- data.frame(Sample = colnames(pred_data), ZT = result_LS$predicted) %>%
      #   mutate(Zone = factor(sub(".*_(Zone_[0-9]+)$", "\\1", Sample), levels = paste0("Zone_", 1:8)),
      #          Sample_ID = order_samples(sub("_Zone_[0-9]+$", "", Sample)),
      #          Group = ifelse(grepl("^M", Sample_ID), "M", "P")) %>%
      #   dplyr::arrange(Sample_ID) %>%
      #   group_by(Sample_ID) %>%
      #   mutate(ZT_diff = ZT - ZT[Zone == "Zone_1"]) %>%
      #   ungroup()
      # p_all <- ggplot(zone_df_LS, aes(Zone, ZT, group = Sample_ID, colour = Sample_ID)) +
      #   geom_line(linewidth = 1) + geom_point(size = 2) +
      #   scale_y_continuous(limits = c(0, 24.5), breaks = 0:24) +
      #   labs(title = "Alle Samples (LS)", x = "Zone", y = "ZT") + theme_bw()
      # p_individual <- ggplot(zone_df_LS, aes(Zone, ZT, group = Sample_ID)) +
      #   geom_line(colour = "steelblue", linewidth = 0.8) + geom_point(size = 2) +
      #   facet_wrap(~Sample_ID, ncol = 4) +
      #   scale_x_discrete(labels = 1:8) +
      #   scale_y_continuous(limits = c(0, 24), breaks = seq(0, 24, 2)) +
      #   labs(x = "Zone", y = "ZT") + theme_bw()
      # print(p_all / p_individual + plot_layout(heights = c(1, 3)))
      # y_range <- range(zone_df_LS$ZT_diff, na.rm = TRUE)
      # make_plot <- function(df, title, colour_by) {
      #   ggplot(df, aes(Zone, ZT_diff, group = Sample_ID, colour = .data[[colour_by]])) +
      #     geom_line(linewidth = 1) + geom_point(size = 2) +
      #     geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
      #     scale_y_continuous(limits = y_range, breaks = pretty(y_range), labels = hour_labels) +
      #     labs(title = title, x = "Zone", y = "\u0394 ZT relativ zu Zone 1") + theme_bw()
      # }
      # for (colour_by in c("Group", "Sample_ID")) {
      #   print(make_plot(zone_df_LS, "Alle Samples (M & P)", colour_by) /
      #           make_plot(filter(zone_df_LS, Group == "M"), "Nur M-Samples", colour_by) /
      #           make_plot(filter(zone_df_LS, Group == "P"), "Nur P-Samples", colour_by))
      # }
      
      # ---- PLOT: Räumliche Plots aller 16 Proben, jede Zone in der Farbe ihrer vorhergesagten ZT (LS) ----
      # plots <- list()
      # for(i in daten){
      #   obj <- readRDS(file.path(path_to_visium, paste0(i, "_Visium"), paste0(i, "_Visium_zoniert.rds")))
      #   zone_time <- round(result_LS$predicted[grep(paste0("^", i, "_Zone_"), colnames(pred_data))])
      #   obj$PredictedTime <- as.numeric(zone_time[obj$med_liver_zone])
      #   plots[[i]] <- SpatialFeaturePlot(obj, features = "PredictedTime", pt.size.factor = 2) +
      #     zt_colour_scale() + ggtitle(i)
      # }
      # print(wrap_plots(plots, ncol = 4, nrow = 4, guides = "collect") & theme(legend.position = "bottom"))
    }
    
    if(bin == TRUE & type == "bin"){
      
      # ---- PLOT: Balkendiagramm der Bin-ZT (LS) pro Sample ----
      # plot_df <- data.frame(Bin = colnames(pred_data), ZT = round(result_LS$predicted) %% 24) %>%
      #   mutate(Sample = order_samples(sub("^(M[0-9]+|P[0-9]+).*", "\\1", Bin)))
      # print(ggplot(plot_df, aes(factor(ZT), fill = factor(ZT))) +
      #         geom_bar(colour = "black") +
      #         facet_wrap(~Sample, ncol = 4) +
      #         scale_fill_manual(values = zt_colours_24, name = "ZT") +
      #         scale_x_discrete(limits = as.character(0:23), breaks = seq(0, 23, by = 2)) +
      #         labs(x = "Zeitgeberzeit (ZT)", y = "Anzahl Bins") +
      #         theme_bw() + theme(legend.position = "bottom"))
      
      # ---- PLOT: Räumliche Plots aller 16 Proben mit der vorhergesagten ZT (LS) pro Bin ----
      # plots <- list()
      # for(i in daten){
      #   obj <- readRDS(file.path(path_to_visium, paste0(i, "_Visium"), paste0(i, "_Visium_zoniert.rds")))
      #   colnames(obj) <- paste0(i, "_", colnames(obj))
      #   obj$PredictedTime <- result_LS$predicted[match(colnames(obj), colnames(pred_data))]
      #   plots[[i]] <- SpatialFeaturePlot(obj, features = "PredictedTime", pt.size.factor = 2) +
      #     zt_colour_scale() + ggtitle(i)
      # }
      # print(wrap_plots(plots, ncol = 4, nrow = 4, guides = "collect") & theme(legend.position = "bottom"))
    }
    
  }
  
}