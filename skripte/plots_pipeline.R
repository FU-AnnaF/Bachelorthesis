## ================================================================
## plots_pipeline.R
## Speicherort: Bachelorarbeit Pipeline/skripte/
##
## erstelle_thesis_plots(): alle Abbildungen fuer die Bachelorarbeit
## als PNG (optional zusaetzlich PDF) in LaTeX-gerechter Groesse.
##
## - liest NUR gespeicherte Ergebnisse (kein Training, kein tauFisher-Lauf)
## - methode waehlt das Modell: "JTK", "LS" (je 170 Spalten) oder
##   "LS_full" (gesamte Trainingstabelle)
## - jede Abbildung wird direkt in Zielgroesse gespeichert
##   (Breite = textbreite_cm). In LaTeX mit width=\textwidth einbinden,
##   dann erscheint die Schrift in der echten Groesse base_size (pt).
## - Labels Englisch, Farben nach Krankheitsstatus wie im SEM&SD-Plot
##
## Erzeugte Dateien (Nummer = Reihenfolge in der Funktion):
##   01_zonierung_alle_samples
##   02a_oszillation_training, 02b_oszillation_bulk_samples
##   03a_pca_space_training, 03_pca_space_bulk
##   04_kreisplot_predicted_zt_bulk
##   05_zonen_linien_pro_sample
##   06_sem_sd_pro_sample
##   07_flatness_R_netchange
##   08a_spots_pro_sample, 08b_spots_raeumlich
##   09_spots_pro_zone/<Sample>.png   (ein Bild pro Sample)
##   09_spots_pro_zone_M1_P6          (M1 oben, P6 unten in einer Abbildung)
##   10a_counts_praediktorgene_summe, 10b_counts_praediktorgene_je_gen
##   11_bootstrap_pseudobulk_zonen    (Pseudobulk-Bootstrap mit Zuruecklegen)
##   12_zonen_vs_prediction_seite1/2  (Zonen neben Zonen-CT, je 8 Samples)
##   13_violin_<GEN>                  (Default: CRY2, PER3)
##   14_delta_zone1_tests             (t-Test & Wilcoxon, M | P)
## ================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(ggrepel)
  library(Seurat)
})


## ---- Gemeinsame Farben / Labels ---------------------------------

sample_level_order <- c(
  "M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8",
  "P7", "P6", "P14", "P18", "P17", "P2", "P21", "P3"
)

sample_condition <- c(
  M1 = "Healthy", M2 = "Healthy", M3 = "Healthy", M4 = "Healthy",
  M5 = "Healthy", M6 = "Healthy", M7 = "Healthy", M8 = "Healthy",
  P6  = "Steatosis", P7  = "Steatosis", P14 = "Steatosis", P18 = "Steatosis",
  P17 = "Severe Steatosis",
  P2  = "Fibrosis", P21 = "Fibrosis",
  P3  = "Portal flow deprivation"
)

condition_colours <- c(
  "Healthy"                 = "#4A90E2",   # helleres Blau (vorher "blue")
  "Steatosis"               = "#2ca02c",
  "Severe Steatosis"        = "#145214",
  "Fibrosis"                = "#ff7f0e",
  "Portal flow deprivation" = "#9467bd"
)
condition_levels <- names(condition_colours)

group_colours   <- c("M" = unname(condition_colours["Healthy"]), "P" = "red")
group_labels    <- c("M" = "Healthy (M)", "P" = "Diseased (P)")
zt_group_labels <- c("vor_12" = "CT 0-11", "nach_12" = "CT 12-23")

zt_anchor_colours <- c(
  "0"  = "#E68A2E", "3"  = "#D9B54C", "6"  = "#8FAE7A", "9"  = "#4C8B8B",
  "12" = "#3D6FA8", "15" = "#6B4E8E", "18" = "#A6528E", "21" = "#D9636B",
  "24" = "#E68A2E"
)
zt_ramp <- colorRampPalette(zt_anchor_colours)
zt_colours_24 <- zt_ramp(25)[1:24]
names(zt_colours_24) <- 0:23
zt_colours_continuous <- zt_ramp(100)

# Zone 1 (central, rot) -> Zone 8 (portal, blau)
zone_colours <- colorRampPalette(c("#B2182B", "#EF8A62", "#FDDBC7",
                                   "#D1E5F0", "#67A9CF", "#2166AC"))(8)
zone_strip_labels <- setNames(
  c("Zone 1 (central)", paste("Zone", 2:7), "Zone 8 (portal)"),
  paste0("Zone_", 1:8)
)

# "M1" -> "M1 (Healthy)"
label_with_condition <- function(x) {
  x_chr <- as.character(x)
  cond  <- sample_condition[x_chr]
  ifelse(is.na(cond), x_chr, paste0(x_chr, " (", cond, ")"))
}
# zweizeilig fuer schmale Facetten: "P3\nPortal flow deprivation"
label_with_condition_2z <- function(x) {
  x_chr <- as.character(x)
  cond  <- sample_condition[x_chr]
  ifelse(is.na(cond), x_chr, paste0(x_chr, "\n", cond))
}

order_samples <- function(x) {
  extra <- sort(setdiff(unique(as.character(x)), sample_level_order))
  factor(x, levels = c(sample_level_order, extra))
}

add_condition <- function(df, col = "Sample_ID") {
  df$Condition <- factor(unname(sample_condition[as.character(df[[col]])]),
                         levels = condition_levels)
  df
}

# tauFisher-Funktion holen (exportiert, selbst gesourct oder intern)
.tf_fun <- function(name) {
  if (exists(name, mode = "function")) get(name, mode = "function")
  else utils::getFromNamespace(name, "tauFisher")
}

# zirkulaere Statistik (24 h)
.circ_mean <- function(x, period = 24) {
  a <- 2 * pi * x / period
  (atan2(mean(sin(a)), mean(cos(a))) * period / (2 * pi)) %% period
}
.circ_R <- function(x, period = 24) {
  a <- 2 * pi * x / period
  sqrt(mean(cos(a))^2 + mean(sin(a))^2)
}
# Net change: zirk. Mittel der letzten 2 Zonen minus der ersten 2 Zonen
.circ_net <- function(x, period = 24) {
  n <- length(x)
  start <- .circ_mean(x[1:2], period)
  ende  <- .circ_mean(x[(n - 1):n], period)
  (ende - start + period / 2) %% period - period / 2
}

.sig_stars <- function(p) {
  dplyr::case_when(is.na(p) ~ "", p < 0.001 ~ "***", p < 0.01 ~ "**",
                   p < 0.05 ~ "*", TRUE ~ "ns")
}
.hour_labels <- function(x) {
  ifelse(x == 0, "0", paste0(ifelse(x > 0, "+", ""), x, " h"))
}
zt_group_alpha <- c("vor_12" = 0.4, "nach_12" = 1)


## ================================================================
## Hauptfunktion
## ================================================================

erstelle_thesis_plots <- function(
    out_dir,
    pfad_pred,                    # Ordner mit <typ>_prediction_<methode>.rds   (ergebnisse/modelle/pred)
    pfad_modell_dir,              # Ordner mit Modellen, train_adj, metacycle   (ergebnisse/modelle/modelle)
    visium_dir,                   # Ordner mit <Sample>_Visium/<Sample>_Visium_zoniert.rds
    methode        = "JTK",       # "JTK" (Standard), "LS" oder "LS_full"
    pfad_modell    = file.path(pfad_modell_dir, paste0("model_talamanca_", methode, ".rds")),
    pfad_train_adj = file.path(pfad_modell_dir, paste0("train_adj_", train_suffix, ".rds")),
    pfad_time_adj  = file.path(pfad_modell_dir, paste0("time_adj_",  train_suffix, ".rds")),
    pfad_meta2d    = file.path(pfad_modell_dir, metacycle_ordner, "Dataset_full_adj_meta2d",
                               "meta2d_Dataset_full_adj.txt"),   # Peak-Phase in den Facetten (optional)
    plots = c("zonierung", "oszillation", "pca_bulk", "kreis_bulk",
              "zonen_linien", "sem_sd", "flatness",
              "spots_sample", "spots_raeumlich", "spots_zonen",
              "counts_praediktorgene", "bootstrap", "zonen_vs_prediction",
              "violin", "test_zone1"),
    formats        = "png",       # z.B. c("png", "pdf")
    dpi            = 300,
    base_size      = 9,           # Schriftgroesse in pt (passend zu 10-11pt Fliesstext)
    textbreite_cm  = 16,          # \the\textwidth deines Dokuments in cm
    titel          = TRUE,        # FALSE = keine Plot-Titel (Titel steht dann in \caption)
    pt_size_factor = 2,
    counts_quantil = 0.995,       # x-Achse der Count-Histogramme bis zu diesem Quantil
    violin_gene    = c("CRY2", "PER3"),
    violin_nur_positive = FALSE,  # TRUE = nur Spots mit Expression > 0 in den Violins
    R_boot         = 500,         # Bootstrap-Wiederholungen pro Sample x Zone (Pseudobulk)
    spots_zonen_kombi = c("M1", "P6")   # Samples untereinander in einer Abbildung (09); NULL = aus
) {
  
  # Methode festlegen (muss vor den Pfad-Defaults ausgewertet werden)
  methode <- match.arg(methode, c("JTK", "LS", "LS_full"))
  # JTK und LS wurden auf denselben 170 Spalten trainiert, LS_full auf allen
  train_suffix     <- if (methode == "LS_full") "LS_full" else "JTK"
  metacycle_ordner <- if (methode == "LS_full") "metacycle_full" else "metacycle"
  # Spalte der Peak-Phase in der meta2d-Ausgabe
  phase_methode    <- if (methode == "JTK") "JTK" else "LS"
  M <- methode   # kurz fuer Titel
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  W   <- textbreite_cm
  ttl <- function(x) if (titel) x else NULL
  th  <- function() {
    theme_bw(base_size = base_size) +
      theme(legend.position = "bottom",
            plot.title = element_text(face = "bold"),
            strip.background = element_rect(fill = "grey95"),
            strip.text = element_text(size = base_size - 1))
  }
  sp_th <- theme(
    plot.title   = element_text(size = base_size - 1, hjust = 0.5, face = "bold"),
    legend.title = element_text(size = base_size - 1),
    legend.text  = element_text(size = base_size - 2)
  )
  txt_size <- (base_size - 1) / ggplot2::.pt
  
  versuche <- function(name, expr) {
    message("\n=== ", name, " (", M, ") ===")
    tryCatch(expr, error = function(e) message("FEHLER bei '", name, "': ", conditionMessage(e)))
  }
  
  ## ---- Speichern -------------------------------------------------
  abbildungen <- list()
  png_dev <- if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png
  } else {
    function(filename, width, height, units = "in", res = 300, ...) {
      grDevices::png(filename, width = width, height = height, units = units,
                     res = res, type = "cairo", ...)
    }
  }
  
  speichern <- function(p, name, height, width = W, ordner = out_dir) {
    abbildungen[[name]] <<- p
    for (fmt in formats) {
      f <- file.path(ordner, paste0(name, ".", fmt))
      tryCatch({
        if (fmt == "png") {
          ggsave(f, p, width = width, height = height, units = "cm",
                 dpi = dpi, device = png_dev, bg = "white")
        } else if (fmt == "pdf") {
          ggsave(f, p, width = width, height = height, units = "cm",
                 device = grDevices::cairo_pdf, bg = "white")
        } else {
          stop("Format nicht unterstuetzt: ", fmt)
        }
        message("  gespeichert: ", f)
      }, error = function(e) message("  FEHLER beim Speichern von ", f, ": ", conditionMessage(e)))
    }
  }
  
  ## ---- Predictions einlesen ---------------------------------------
  lies_pred <- function(typ) {
    f <- file.path(pfad_pred, paste0(typ, "_prediction_", methode, ".rds"))
    if (!file.exists(f)) { message("Hinweis: fehlt ", f); return(NULL) }
    d <- readRDS(f)
    if ("type" %in% names(d)) d <- d[d$type == typ, , drop = FALSE]
    d
  }
  bulk_pred  <- lies_pred("bulk")
  zonen_pred <- lies_pred("zonen")
  bin_pred   <- lies_pred("bin")
  
  bulk_df <- NULL; bulk_lookup <- NULL
  if (!is.null(bulk_pred)) {
    bulk_df <- data.frame(Sample_ID = as.character(bulk_pred$ID),
                          ZT = bulk_pred$ZT_predicted %% 24,
                          stringsAsFactors = FALSE) %>% add_condition()
    bulk_lookup <- setNames(bulk_df$ZT, bulk_df$Sample_ID)
  }
  
  zonen_df <- NULL
  if (!is.null(zonen_pred)) {
    zonen_df <- data.frame(
      Sample_ID = sub("_Zone_[0-9]+$", "", zonen_pred$ID),
      Zone      = factor(sub("^.*_(Zone_[0-9]+)$", "\\1", zonen_pred$ID),
                         levels = paste0("Zone_", 1:8)),
      ZT        = zonen_pred$ZT_predicted,
      stringsAsFactors = FALSE
    ) %>%
      add_condition() %>%
      dplyr::arrange(order_samples(Sample_ID), Zone)
  }
  
  bin_df <- NULL
  if (!is.null(bin_pred)) {
    bin_df <- data.frame(
      ID        = as.character(bin_pred$ID),
      Sample_ID = sub("^(M[0-9]+|P[0-9]+).*", "\\1", bin_pred$ID),
      ZT        = round(bin_pred$ZT_predicted) %% 24,
      stringsAsFactors = FALSE
    )
  }
  
  ## ---- Modell -------------------------------------------------------
  modell <- NULL
  if (any(c("oszillation", "pca_bulk", "counts_praediktorgene", "bootstrap") %in% plots)) {
    if (file.exists(pfad_modell)) {
      modell <- readRDS(pfad_modell)
    } else {
      message("Hinweis: Modell fehlt (", pfad_modell, ")")
    }
  }
  pred_genes <- if (!is.null(modell$Genes)) toupper(modell$Genes) else
    c("ABHD17C", "CRY1", "CRY2", "FANCE", "NPAS2", "NR1D2",
      "NR6A1", "PER3", "PPARD", "RIC8B", "SPACA4")
  
  ## ---- Seurat-Objekte EINMAL pro Sample laden -----------------------
  brauche_seurat <- any(c("zonierung", "oszillation", "pca_bulk", "spots_raeumlich",
                          "spots_zonen", "counts_praediktorgene", "bootstrap",
                          "zonen_vs_prediction", "violin") %in% plots)
  
  zon_plots <- list(); spot_plots <- list(); pb_list <- list()
  spot_meta <- list(); gene_counts <- list()
  pair_plots <- list(); expr_list <- list(); boot_dat <- list()
  
  if (brauche_seurat) {
    for (s in sample_level_order) {
      f <- file.path(visium_dir, paste0(s, "_Visium"), paste0(s, "_Visium_zoniert.rds"))
      if (!file.exists(f)) { message("Hinweis: Visium-Objekt fehlt: ", f); next }
      message("Lade ", s, " ...")
      obj <- readRDS(f)
      cnt <- tryCatch(Seurat::GetAssayData(obj, assay = "Spatial", layer = "counts"),
                      error = function(e) Seurat::GetAssayData(obj, assay = "Spatial", slot = "counts"))
      ids <- paste0(s, "_", colnames(obj))
      
      # Pseudobulk wie in der Pipeline: Summe der Counts ueber alle Spots
      pb_list[[s]] <- Matrix::rowSums(cnt)
      
      spot_meta[[s]] <- data.frame(ID = ids, Zone_num = as.numeric(obj$med_liver_zone),
                                   stringsAsFactors = FALSE)
      
      if ("counts_praediktorgene" %in% plots && !is.null(bin_df)) {
        idx  <- match(pred_genes, toupper(rownames(cnt)))
        g_ok <- pred_genes[!is.na(idx)]
        keep <- which(ids %in% bin_df$ID)
        if (length(g_ok) > 0 && length(keep) > 0) {
          m <- as.matrix(cnt[idx[!is.na(idx)], keep, drop = FALSE])
          gene_counts[[s]] <- data.frame(
            ID        = rep(ids[keep], each = length(g_ok)),
            Sample_ID = s,
            Gene      = rep(g_ok, times = length(keep)),
            Counts    = as.vector(m),
            stringsAsFactors = FALSE
          )
        }
      }
      
      if ("zonierung" %in% plots) {
        obj_z <- subset(obj, subset = !is.na(liver_Zone))
        zon_plots[[s]] <- suppressMessages(
          SpatialFeaturePlot(obj_z, features = "med_liver_zone", pt.size.factor = pt_size_factor) +
            scale_fill_gradientn(colours = zone_colours, limits = c(1, 8), breaks = 1:8,
                                 name = "Zone (1 = central, 8 = portal)") +
            ggtitle(label_with_condition(s)) + sp_th
        )
        rm(obj_z)
      }
      
      if ("spots_raeumlich" %in% plots && !is.null(bin_df)) {
        obj$PredictedTime <- bin_df$ZT[match(ids, bin_df$ID)]
        spot_plots[[s]] <- suppressMessages(
          SpatialFeaturePlot(obj, features = "PredictedTime", pt.size.factor = pt_size_factor) +
            scale_fill_gradientn(colours = zt_colours_continuous, limits = c(0, 24),
                                 breaks = seq(0, 24, 4), oob = scales::squish,
                                 na.value = "grey85", name = "Predicted CT") +
            ggtitle(label_with_condition(s)) + sp_th
        )
      }
      
      # Praediktorgen-Counts der zonierten Spots fuer den Pseudobulk-Bootstrap
      if ("bootstrap" %in% plots && !is.null(modell$Genes)) {
        zone_v <- suppressWarnings(as.integer(as.character(obj$med_liver_zone)))
        keep_z <- which(!is.na(zone_v) & zone_v >= 1 & zone_v <= 8)
        idx_b  <- match(toupper(modell$Genes), toupper(rownames(cnt)))
        Xb <- cnt[ifelse(is.na(idx_b), 1L, idx_b), keep_z, drop = FALSE]
        if (any(is.na(idx_b))) {
          Xb[is.na(idx_b), ] <- 0
          message("  ", s, ": fehlende Praediktorgene als 0: ",
                  paste(modell$Genes[is.na(idx_b)], collapse = ", "))
        }
        Xb <- Matrix::drop0(methods::as(Xb, "CsparseMatrix"))
        dimnames(Xb) <- list(toupper(modell$Genes), colnames(obj)[keep_z])
        boot_dat[[s]] <- list(X = Xb, zone = zone_v[keep_z])
      }
      
      # Zonen neben Zonen-Prediction
      if ("zonen_vs_prediction" %in% plots && !is.null(zonen_df)) {
        zt_s <- zonen_df[zonen_df$Sample_ID == s, ]
        zone_zt <- setNames(round(zt_s$ZT) %% 24, sub("Zone_", "", as.character(zt_s$Zone)))
        obj$ZoneTime <- unname(zone_zt[as.character(obj$med_liver_zone)])
        obj_z <- subset(obj, subset = !is.na(liver_Zone))
        p_zone <- suppressMessages(
          SpatialFeaturePlot(obj_z, features = "med_liver_zone", pt.size.factor = pt_size_factor) +
            scale_fill_gradientn(colours = zone_colours, limits = c(1, 8), breaks = 1:8,
                                 na.value = "grey85", name = "Zone (1 = central, 8 = portal)") +
            ggtitle(paste0(label_with_condition(s), ": zones")) + sp_th
        )
        p_time <- suppressMessages(
          SpatialFeaturePlot(obj_z, features = "ZoneTime", pt.size.factor = pt_size_factor,
                             image.alpha = 0) +   # Gewebebild im Hintergrund aus
            scale_fill_gradientn(colours = zt_colours_continuous, limits = c(0, 24),
                                 breaks = seq(0, 24, 4), oob = scales::squish,
                                 na.value = "grey85", name = "Predicted CT (zone)") +
            ggtitle(paste0(s, ": predicted CT")) + sp_th
        )
        pair_plots[[s]] <- p_zone | p_time
        rm(obj_z)
      }
      
      # Expression (log-normalisiert wie NormalizeData) fuer die Violin-Gene
      if ("violin" %in% plots && !is.null(bin_df)) {
        idx  <- match(toupper(violin_gene), toupper(rownames(cnt)))
        g_ok <- toupper(violin_gene)[!is.na(idx)]
        keep <- which(ids %in% bin_df$ID)
        if (length(g_ok) > 0 && length(keep) > 0) {
          lib <- Matrix::colSums(cnt)[keep]
          lib[lib == 0] <- NA
          m <- as.matrix(cnt[idx[!is.na(idx)], keep, drop = FALSE])
          m <- log1p(sweep(m, 2, lib, "/") * 1e4)
          expr_list[[s]] <- data.frame(
            ID         = rep(ids[keep], each = length(g_ok)),
            Sample_ID  = s,
            Gene       = rep(g_ok, times = length(keep)),
            Expression = as.vector(m),
            stringsAsFactors = FALSE
          )
        }
      }
      
      rm(obj, cnt); gc(verbose = FALSE)
    }
  }
  
  # Pseudobulk-Matrix Gene x Samples
  pb_mat <- NULL
  if (length(pb_list) > 0) {
    genes_all <- Reduce(union, lapply(pb_list, names))
    pb_mat <- sapply(pb_list, function(v) { x <- v[genes_all]; x[is.na(x)] <- 0; x })
    rownames(pb_mat) <- genes_all
  }
  
  # Bin -> Zone
  bin_zone_df <- NULL
  if (!is.null(bin_df) && length(spot_meta) > 0) {
    bin_zone_df <- bin_df %>%
      dplyr::inner_join(dplyr::bind_rows(spot_meta), by = "ID") %>%
      dplyr::filter(!is.na(Zone_num), Zone_num >= 1, Zone_num <= 8) %>%
      dplyr::mutate(Zone = factor(paste0("Zone_", Zone_num), levels = paste0("Zone_", 1:8)))
  }
  
  
  ## ================================================================
  ## 01 Zonierung aller Samples
  ## ================================================================
  if ("zonierung" %in% plots) versuche("zonierung", {
    if (length(zon_plots) == 0) stop("keine Zonierungsplots erzeugt")
    p <- wrap_plots(zon_plots, ncol = 4, guides = "collect") &
      theme(legend.position = "bottom",
            legend.key.width = unit(1.5, "cm"), legend.key.height = unit(0.3, "cm"))
    if (titel) p <- p + plot_annotation(
      title = "Liver zonation of all samples",
      theme = theme(plot.title = element_text(face = "bold", size = base_size + 1)))
    speichern(p, "01_zonierung_alle_samples", height = W * 1.15)
  })
  
  
  ## ================================================================
  ## 02 Oszillation der Praediktorgene (Training + Bulk-Samples)
  ## ================================================================
  if ("oszillation" %in% plots) versuche("oszillation", {
    if (!file.exists(pfad_train_adj) || !file.exists(pfad_time_adj))
      stop(basename(pfad_train_adj), " / ", basename(pfad_time_adj), " fehlen in ", dirname(pfad_train_adj))
    
    tr <- as.matrix(readRDS(pfad_train_adj))
    tt <- as.numeric(readRDS(pfad_time_adj))
    if (length(tt) != ncol(tr)) stop("time_adj passt nicht zu train_adj")
    
    ok <- colSums(is.na(tr)) == 0
    tr <- tr[, ok, drop = FALSE]; tt <- tt[ok]
    rownames(tr) <- toupper(rownames(tr))
    
    genes  <- sort(intersect(pred_genes, rownames(tr)))
    df_log <- log2(t(tr[genes, , drop = FALSE]) + 1)
    
    fda <- .tf_fun("get_FDAcurves")(df_log, time = tt, numbasis = 5)
    fda_long <- fda %>%
      tidyr::pivot_longer(-fda_time, names_to = "Gene", values_to = "Expr") %>%
      dplyr::mutate(Gene = toupper(Gene), ZT = fda_time %% 24, Cycle = floor(fda_time / 24)) %>%
      dplyr::group_by(Gene) %>%
      dplyr::mutate(z = (Expr - mean(Expr)) / sd(Expr)) %>%
      dplyr::ungroup()
    
    # Gene nach Peak der FDA-Kurve sortieren
    gene_order <- fda_long %>%
      dplyr::mutate(ZT_h = round(ZT) %% 24) %>%
      dplyr::group_by(Gene, ZT_h) %>%
      dplyr::summarise(z = mean(z), .groups = "drop_last") %>%
      dplyr::slice_max(z, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(ZT_h) %>%
      dplyr::pull(Gene)
    gene_order <- c(gene_order, setdiff(genes, gene_order))
    
    # Peak time je Gen aus der FDA-Kurve
    # Kurve ueber alle Zyklen auf ein 0.1-h-Raster gemittelt, Maximum = Peak
    peak_df <- fda_long %>%
      dplyr::mutate(ZT_r = round(ZT, 1) %% 24) %>%
      dplyr::group_by(Gene, ZT_r) %>%
      dplyr::summarise(z = mean(z), .groups = "drop_last") %>%
      dplyr::slice_max(z, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::transmute(Gene, Peak = ZT_r)
    
    lab <- setNames(sprintf("%s\npeak CT %.1f", peak_df$Gene, peak_df$Peak), peak_df$Gene)
    lab <- c(lab, setNames(setdiff(genes, names(lab)), setdiff(genes, names(lab))))
    
    # optional: zusaetzlich die Phase aus meta2d (JTK_adjphase bzw. LS_adjphase)
    phase_col <- paste0(phase_methode, "_adjphase")
    if (!is.null(pfad_meta2d) && file.exists(pfad_meta2d)) {
      m2d <- utils::read.delim(pfad_meta2d)
      m2d$Gene <- toupper(m2d$CycID)
      m2d <- m2d[m2d$Gene %in% names(lab) & !duplicated(m2d$Gene), ]
      if (phase_col %in% names(m2d)) {
        lab[m2d$Gene] <- sprintf("%s (%s phase %.1f)", lab[m2d$Gene], phase_methode,
                                 m2d[[phase_col]] %% 24)
      }
    }
    peak_df$Gene <- factor(peak_df$Gene, levels = gene_order)
    
    fda_long$Gene <- factor(fda_long$Gene, levels = gene_order)
    pts <- data.frame(
      ZT   = rep(tt %% 24, times = length(genes)),
      Gene = factor(rep(genes, each = nrow(df_log)), levels = gene_order),
      Expr = as.vector(df_log)
    )
    n_rows <- ceiling(length(genes) / 4)
    
    p1 <- ggplot() +
      geom_point(data = pts, aes(ZT, Expr), size = 0.6, alpha = 0.4) +
      geom_vline(data = peak_df, aes(xintercept = Peak), linetype = "dashed",
                 colour = "grey50", linewidth = 0.4) +
      geom_line(data = fda_long, aes(ZT, Expr, group = Cycle), linewidth = 0.6) +
      facet_wrap(~ Gene, ncol = 4, scales = "free_y", labeller = as_labeller(lab)) +
      scale_x_continuous(limits = c(0, 24), breaks = seq(0, 24, 6)) +
      labs(title = ttl(paste0("Oscillation of the tauFisher predictor genes (", M, ") in the training data")),
           x = "CT (h)", y = "log2(expression + 1)") +
      th()
    speichern(p1, "02a_oszillation_training", height = 3.6 * n_rows + 1.5)
    
    if (!is.null(pb_mat) && !is.null(bulk_lookup)) {
      te <- pb_mat; rownames(te) <- toupper(rownames(te))
      g_te <- intersect(genes, rownames(te))
      te_z <- t(scale(t(log2(te[g_te, , drop = FALSE] + 1))))
      te_df <- data.frame(
        Gene      = factor(rep(g_te, times = ncol(te_z)), levels = gene_order),
        Sample_ID = rep(colnames(te_z), each = nrow(te_z)),
        z         = as.vector(te_z),
        stringsAsFactors = FALSE
      ) %>% add_condition()
      te_df$ZT <- unname(bulk_lookup[te_df$Sample_ID])
      
      p2 <- ggplot() +
        geom_vline(data = peak_df, aes(xintercept = Peak), linetype = "dashed",
                   colour = "grey50", linewidth = 0.4) +
        geom_line(data = fda_long, aes(ZT, z, group = Cycle), linewidth = 0.6) +
        geom_point(data = te_df, aes(ZT, z, fill = Condition), shape = 21, size = 1.8, alpha = 0.9,
                   position = position_jitter(width = 0.2, height = 0, seed = 1)) +
        facet_wrap(~ Gene, ncol = 4, labeller = as_labeller(lab)) +
        scale_fill_manual(values = condition_colours, name = "Condition", drop = TRUE) +
        scale_x_continuous(limits = c(-0.5, 24.5), breaks = seq(0, 24, 6)) +
        labs(title = ttl(paste0("Visium bulk samples at their predicted CT vs. training curve (", M, ")")),
             x = "CT (h)", y = "z-score") +
        th() + guides(fill = guide_legend(nrow = 2))
      speichern(p2, "02b_oszillation_bulk_samples", height = 3.6 * n_rows + 2.5)
    } else {
      message("Hinweis: Pseudobulk oder bulk_prediction fehlt - 02b uebersprungen")
    }
  })
  
  
  ## ================================================================
  ## 03 tauFisher PCA space, Bulk-Ebene
  ## ================================================================
  if ("pca_bulk" %in% plots) versuche("pca_bulk", {
    if (is.null(modell$PCA_embedding) || is.null(modell$Model))
      stop("Modell unvollstaendig (PCA_embedding/Model fehlt): ", pfad_modell)
    if (is.null(pb_mat)) stop("Pseudobulk fehlt (Visium-Objekte nicht geladen)")
    # predict() fuer multinom ist eine S3-Methode aus nnet -> Namespace muss geladen sein
    if (!requireNamespace("nnet", quietly = TRUE)) stop("Paket 'nnet' ist nicht installiert")
    
    f_diff  <- .tf_fun("create_DiffMatrix")
    f_scale <- .tf_fun("scale_DiffMatrix")
    
    # exakt die Vorverarbeitung aus test_tauFisher()
    te <- pb_mat; rownames(te) <- toupper(rownames(te))
    genes <- modell$Genes[order(modell$Genes)]
    df <- t(data.frame(te[modell$Genes, , drop = FALSE], check.names = FALSE))
    df <- df[, order(colnames(df)), drop = FALSE]
    proj <- stats::predict(modell$PCA_embedding, newdata = f_scale(f_diff(genes, log2(df + 1))))
    
    pca <- modell$PCA_embedding
    classify <- function(m) {
      d <- as.data.frame(m[, 1:2, drop = FALSE]); names(d) <- c("PC1", "PC2")
      as.integer(as.character(predict(modell$Model, newdata = d, type = "class"))) %% 24
    }
    
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
    
    test_df <- data.frame(Sample_ID = rownames(proj), PC1 = proj[, 1], PC2 = proj[, 2],
                          stringsAsFactors = FALSE) %>% add_condition()
    
    # Grundplot: Entscheidungsregionen + Trainingsdaten (gleiche Achsen fuer beide Versionen)
    p_basis <- ggplot() +
      geom_raster(data = grid, aes(PC1, PC2), fill = grid_hex, alpha = 0.25) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
      geom_point(data = train_df, aes(PC1, PC2, colour = ZT), size = 0.9, alpha = 0.7) +
      geom_text(data = hour_lab, aes(PC1, PC2, label = ZT, colour = ZT),
                size = base_size / ggplot2::.pt, fontface = "bold") +
      scale_colour_manual(values = zt_colours_24, guide = "none", drop = FALSE) +
      coord_equal(expand = FALSE) +
      labs(x = paste0("PC1 (", ve[1], " %)"), y = paste0("PC2 (", ve[2], " %)")) +
      th()
    
    # Version OHNE projizierte Bulk-Punkte
    p_ohne <- p_basis +
      labs(title = ttl(paste0("tauFisher PCA space (", M, ") of the training data")))
    # coord_equal: Hoehe aus dem Seitenverhaeltnis der Daten ableiten,
    # sonst entsteht oben/unten leerer Rand
    h_panel <- (W - 1.5) * diff(y_rng) / diff(x_rng)
    h_ohne  <- min(max(h_panel + (if (titel) 2.5 else 1.7), 6), 22)
    h_bulk  <- h_ohne + 1.2                              # + Condition-Legende
    speichern(p_ohne, "03a_pca_space_training", height = h_ohne)
    
    p <- p_basis +
      geom_point(data = test_df, aes(PC1, PC2, fill = Condition),
                 shape = 21, size = 2.6, colour = "black") +
      geom_text_repel(data = test_df, aes(PC1, PC2, label = Sample_ID),
                      size = txt_size, max.overlaps = Inf, min.segment.length = 0) +
      scale_fill_manual(values = condition_colours, name = "Condition", drop = TRUE) +
      labs(title = ttl(paste0("tauFisher PCA space (", M, ") with projected bulk samples"))) +
      guides(fill = guide_legend(nrow = 2))
    speichern(p, "03_pca_space_bulk", height = h_bulk)
  })
  
  
  ## ================================================================
  ## 04 Kreisplot predicted CT pro Sample (Bulk)
  ## ================================================================
  if ("kreis_bulk" %in% plots) versuche("kreis_bulk", {
    if (is.null(bulk_df)) stop("bulk_prediction fehlt")
    kd <- bulk_df %>% dplyr::mutate(Sample_ID = droplevels(order_samples(Sample_ID)))
    
    ## Uhr in kartesischen Koordinaten: jedes Sample hat einen eigenen Ring
    ## (aussen M1 ... innen P3), CT 0 oben, im Uhrzeigersinn. Jeder Name steht
    ## neben seinem Punkt, entlang seines Rings gedreht -> keine Ueberlagerung.
    loch  <- 2                                    # leerer Innenkreis
    n_s   <- nlevels(kd$Sample_ID)
    r_max <- loch + n_s
    kd <- kd %>%
      dplyr::mutate(r = loch + (n_s + 1 - as.integer(Sample_ID)),   # M1 = aussen
                    theta = 2 * pi * ZT / 24,
                    x = r * sin(theta), y = r * cos(theta))
    
    # Label kurz hinter dem Punkt auf demselben Ring, tangential gedreht;
    # in der unteren Haelfte um 180 Grad gedreht, damit er nicht auf dem Kopf steht
    namen <- kd %>%
      dplyr::mutate(
        theta_l = theta + 0.55 / r,
        xl = r * sin(theta_l), yl = r * cos(theta_l),
        unten = cos(theta_l) < 0,
        winkel = -theta_l * 180 / pi + ifelse(unten, 180, 0),
        hj = ifelse(unten, 1, 0)
      )
    
    kreis <- function(r) {
      t <- seq(0, 2 * pi, length.out = 361)
      data.frame(x = r * sin(t), y = r * cos(t), r = r)
    }
    ringe   <- dplyr::bind_rows(lapply(loch + seq_len(n_s), kreis))
    aussen  <- kreis(r_max + 0.7)
    # jede Stunde ein Strich + Beschriftung; alle 3 h kraeftiger und mit "CT"
    stunden <- data.frame(h = 0:23) %>%
      dplyr::mutate(haupt = h %% 3 == 0,
                    theta = 2 * pi * h / 24,
                    xend = (r_max + 0.7) * sin(theta), yend = (r_max + 0.7) * cos(theta),
                    xl = (r_max + 2) * sin(theta), yl = (r_max + 2) * cos(theta),
                    lab = ifelse(haupt, paste0("CT ", h), as.character(h)))
    
    p <- ggplot() +
      geom_path(data = ringe, aes(x, y, group = r), colour = "grey88", linewidth = 0.3) +
      geom_path(data = aussen, aes(x, y), colour = "grey40", linewidth = 0.5) +
      geom_segment(data = dplyr::filter(stunden, !haupt),
                   aes(x = 0, y = 0, xend = xend, yend = yend),
                   colour = "grey92", linewidth = 0.25) +
      geom_segment(data = dplyr::filter(stunden, haupt),
                   aes(x = 0, y = 0, xend = xend, yend = yend),
                   colour = "grey75", linewidth = 0.35) +
      geom_text(data = dplyr::filter(stunden, haupt), aes(xl, yl, label = lab),
                size = txt_size, fontface = "bold") +
      geom_text(data = dplyr::filter(stunden, !haupt), aes(xl, yl, label = lab),
                size = (base_size - 2) / ggplot2::.pt, colour = "grey35") +
      geom_point(data = kd, aes(x, y, colour = Condition), size = 2) +
      geom_text(data = namen,
                aes(xl, yl, label = Sample_ID, colour = Condition, angle = winkel, hjust = hj),
                vjust = 0.5, size = (base_size - 2) / ggplot2::.pt, fontface = "bold",
                show.legend = FALSE) +
      scale_colour_manual(values = condition_colours, name = "Condition", drop = TRUE) +
      coord_equal(xlim = c(-r_max - 3, r_max + 3), ylim = c(-r_max - 3, r_max + 3)) +
      labs(title = ttl(paste0("Predicted circadian phase (", M, ") per sample"))) +
      theme_void(base_size = base_size) +
      theme(legend.position = "bottom",
            plot.title = element_text(face = "bold", hjust = 0.5),
            plot.caption = element_text(colour = "grey40", size = base_size - 2)) +
      guides(colour = guide_legend(nrow = 2))
    speichern(p, "04_kreisplot_predicted_zt_bulk", height = W * 1.08, width = W)
  })
  
  
  ## ================================================================
  ## 05 Linienplot pro Sample mit Zonen-Prediction (+ Bulk-CT)
  ## ================================================================
  if ("zonen_linien" %in% plots) versuche("zonen_linien", {
    if (is.null(zonen_df)) stop("zonen_prediction fehlt")
    zd <- zonen_df %>% dplyr::mutate(Sample_ID = order_samples(Sample_ID))
    
    p <- ggplot(zd, aes(Zone, ZT, group = Sample_ID, colour = Condition)) +
      geom_line(linewidth = 0.7) +
      geom_point(size = 1.3) +
      facet_wrap(~ Sample_ID, ncol = 4) +   # Krankheit steht nur in der Farblegende
      scale_x_discrete(labels = 1:8) +
      scale_y_continuous(limits = c(0, 24), breaks = seq(0, 24, 6)) +
      scale_colour_manual(values = condition_colours, name = "Condition", drop = TRUE) +
      labs(title = ttl(paste0("Predicted circadian phase (", M, ") per zone and sample")),
           x = "Zone (1 = central, 8 = portal)", y = "Predicted CT") +
      th() + guides(colour = guide_legend(nrow = 2))
    
    if (!is.null(bulk_df)) {
      bl <- bulk_df %>% dplyr::mutate(Sample_ID = order_samples(Sample_ID))
      p <- p +
        geom_hline(data = bl, aes(yintercept = ZT, linetype = "Sample CT"),
                   colour = "grey40", linewidth = 0.4, inherit.aes = FALSE) +
        scale_linetype_manual(values = c("Sample CT" = "dashed"), name = NULL)
    }
    speichern(p, "05_zonen_linien_pro_sample", height = W * 1.05)
  })
  
  
  ## ================================================================
  ## 06 SEM & SD der 8 Zonenwerte pro Sample
  ## ================================================================
  if ("sem_sd" %in% plots) versuche("sem_sd", {
    if (is.null(zonen_df)) stop("zonen_prediction fehlt")
    sd_df <- zonen_df %>%
      dplyr::group_by(Sample_ID, Condition) %>%
      dplyr::summarise(n = dplyr::n(), mean = mean(ZT), sd = sd(ZT), .groups = "drop") %>%
      dplyr::mutate(sem = sd / sqrt(n)) %>%
      dplyr::arrange(mean) %>%
      dplyr::mutate(Sample_ID = factor(Sample_ID, levels = Sample_ID))
    
    p <- ggplot(sd_df, aes(Sample_ID, mean, colour = Condition)) +
      geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.25, linewidth = 0.5, alpha = 0.5) +
      geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.15, linewidth = 0.9) +
      geom_point(size = 2) +
      scale_colour_manual(values = condition_colours, name = "Condition", drop = TRUE) +
      labs(title = ttl(paste0("SEM & SD of the 8 zone values per sample (", M, ")")),
           subtitle = ttl("Outer/light = \u00b1SD, inner/dark = \u00b1SEM (n = 8 zones)"),
           x = "Sample", y = "Predicted CT (mean across zones)") +
      th() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      guides(colour = guide_legend(nrow = 2))
    speichern(p, "06_sem_sd_pro_sample", height = W)
  })
  
  
  ## ================================================================
  ## 07 Flatness (R) und Net Change pro Sample
  ## ================================================================
  if ("flatness" %in% plots) versuche("flatness", {
    if (is.null(zonen_df)) stop("zonen_prediction fehlt")
    fl <- zonen_df %>%
      dplyr::arrange(Sample_ID, Zone) %>%
      dplyr::group_by(Sample_ID) %>%
      dplyr::summarise(R = .circ_R(ZT), net_change = .circ_net(ZT), .groups = "drop") %>%
      add_condition() %>%
      dplyr::arrange(dplyr::desc(R)) %>%
      dplyr::mutate(Sample_ID = factor(Sample_ID, levels = Sample_ID))
    
    nc_max <- max(abs(fl$net_change), 1, na.rm = TRUE)
    
    p_R <- ggplot(fl, aes(Sample_ID, R, fill = Condition)) +
      geom_col(width = 0.7) +
      geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
      scale_fill_manual(values = condition_colours, name = "Condition", drop = TRUE) +
      coord_cartesian(ylim = c(0, 1)) +
      labs(x = NULL, y = "R (flatness)") +
      th() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    
    p_nc <- ggplot(fl, aes(Sample_ID, net_change, fill = Condition)) +
      geom_col(width = 0.7) +
      geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.3) +
      scale_fill_manual(values = condition_colours, name = "Condition", drop = TRUE) +
      coord_cartesian(ylim = c(-nc_max, nc_max)) +
      labs(x = "Sample (sorted by R)", y = "Net change\n(zones 7-8 vs. 1-2, h)") +
      th() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    p <- (p_R / p_nc) + plot_layout(guides = "collect") &
      theme(legend.position = "bottom")
    if (titel) p <- p + plot_annotation(
      title = paste0("Flatness of the zonal CT profile (circular R) and net change (", M, ")"),
      theme = theme(plot.title = element_text(face = "bold", size = base_size + 1)))
    speichern(p, "07_flatness_R_netchange", height = W * 0.7)
  })
  
  ## ================================================================
  ## 08a Spot-Prediction pro Sample (Balken) / 08b raeumlich
  ## ================================================================
  if ("spots_sample" %in% plots) versuche("spots_sample", {
    if (is.null(bin_df)) stop("bin_prediction fehlt")
    bd <- bin_df %>%
      dplyr::mutate(Sample_ID = order_samples(Sample_ID),
                    ZT = factor(ZT, levels = 0:23))
    
    p <- ggplot(bd, aes(ZT, fill = ZT)) +
      geom_bar(colour = "black", linewidth = 0.1) +
      facet_wrap(~ Sample_ID, ncol = 4, scales = "free_y",
                 labeller = as_labeller(label_with_condition_2z)) +
      scale_fill_manual(values = zt_colours_24, guide = "none", drop = FALSE) +
      scale_x_discrete(drop = FALSE, breaks = as.character(seq(0, 23, 4))) +
      labs(title = ttl(paste0("Predicted CT (", M, ") of all spots, per sample")),
           x = "Predicted CT", y = "Number of spots") +
      th()
    
    if (!is.null(bulk_df)) {
      bl <- bulk_df %>% dplyr::mutate(Sample_ID = order_samples(Sample_ID),
                                      x = round(ZT) %% 24 + 1)
      p <- p +
        geom_vline(data = bl, aes(xintercept = x, linetype = "Sample CT"),
                   colour = "grey30", linewidth = 0.4, inherit.aes = FALSE) +
        scale_linetype_manual(values = c("Sample CT" = "dashed"), name = NULL)
    }
    speichern(p, "08a_spots_pro_sample", height = W * 1.05)
  })
  
  if ("spots_raeumlich" %in% plots) versuche("spots_raeumlich", {
    if (length(spot_plots) == 0) stop("keine raeumlichen Spot-Plots erzeugt")
    p <- wrap_plots(spot_plots, ncol = 4, guides = "collect") &
      theme(legend.position = "bottom",
            legend.key.width = unit(1.5, "cm"), legend.key.height = unit(0.3, "cm"))
    if (titel) p <- p + plot_annotation(
      title = paste0("Predicted CT (", M, ") per spot"),
      theme = theme(plot.title = element_text(face = "bold", size = base_size + 1)))
    speichern(p, "08b_spots_raeumlich", height = W * 1.15)
  })
  
  
  ## ================================================================
  ## 09 Spot-Prediction pro Sample, aufgeteilt in Zonen (1 Bild/Sample)
  ## ================================================================
  if ("spots_zonen" %in% plots) versuche("spots_zonen", {
    if (is.null(bin_zone_df)) stop("Bin-Zonen-Zuordnung fehlt")
    dir_z <- file.path(out_dir, "09_spots_pro_zone")
    dir.create(dir_z, showWarnings = FALSE)
    
    samples <- sample_level_order[sample_level_order %in% unique(bin_zone_df$Sample_ID)]
    zonen_plots <- list()
    for (s in samples) {
      df_s <- bin_zone_df %>%
        dplyr::filter(Sample_ID == s) %>%
        dplyr::mutate(ZT = factor(ZT, levels = 0:23))
      
      p <- ggplot(df_s, aes(ZT, fill = ZT)) +
        geom_bar(colour = "black", linewidth = 0.1) +
        facet_wrap(~ Zone, ncol = 4, drop = FALSE, labeller = as_labeller(zone_strip_labels)) +
        scale_fill_manual(values = zt_colours_24, guide = "none", drop = FALSE) +
        scale_x_discrete(drop = FALSE, breaks = as.character(seq(0, 23, 4))) +
        labs(title = ttl(paste0(paste0("Predicted CT (", M, ") of spots per zone - "), label_with_condition(s))),
             x = "Predicted CT", y = "Number of spots") +
        th()
      
      if (!is.null(bulk_lookup) && s %in% names(bulk_lookup)) {
        p <- p + geom_vline(data = data.frame(x = round(bulk_lookup[[s]]) %% 24 + 1, Linie = "Sample CT"),
                            aes(xintercept = x, linetype = Linie),
                            colour = "grey30", linewidth = 0.5, inherit.aes = FALSE)
      }
      if (!is.null(zonen_df)) {
        zl <- zonen_df %>%
          dplyr::filter(Sample_ID == s) %>%
          dplyr::transmute(Zone, x = round(ZT) %% 24 + 1, Linie = "Zone CT")
        if (nrow(zl) > 0) {
          p <- p + geom_vline(data = zl, aes(xintercept = x, linetype = Linie),
                              colour = "black", linewidth = 0.6, inherit.aes = FALSE)
        }
      }
      p <- p + scale_linetype_manual(values = c("Sample CT" = "dashed", "Zone CT" = "dotted"), name = NULL)
      
      speichern(p, s, height = W * 0.6, ordner = dir_z)
      zonen_plots[[s]] <- p
    }
    
    # mehrere Samples untereinander in EINER Abbildung (Default M1 oben, P6 unten)
    kombi <- spots_zonen_kombi[spots_zonen_kombi %in% names(zonen_plots)]
    if (length(kombi) >= 2) {
      teile <- lapply(kombi, function(s) {
        zonen_plots[[s]] + labs(title = label_with_condition(s)) +
          theme(plot.title = element_text(size = base_size, face = "bold"))
      })
      p_k <- wrap_plots(teile, ncol = 1) +
        plot_layout(guides = "collect") &
        theme(legend.position = "bottom")
      if (titel) p_k <- p_k + plot_annotation(
        title = paste0("Predicted circadian phase (", M, ") of spots per zone"),
        theme = theme(plot.title = element_text(face = "bold", size = base_size + 1)))
      speichern(p_k, paste0("09_spots_pro_zone_", paste(kombi, collapse = "_")),
                height = W * 0.55 * length(kombi) + 1.5)
    } else if (length(spots_zonen_kombi) > 0) {
      message("Hinweis: fuer die Kombi-Abbildung nicht alle Samples vorhanden: ",
              paste(setdiff(spots_zonen_kombi, names(zonen_plots)), collapse = ", "))
    }
  })
  
  
  ## ================================================================
  ## 10 Counts der Praediktorgene pro Spot
  ## ================================================================
  if ("counts_praediktorgene" %in% plots) versuche("counts_praediktorgene", {
    if (length(gene_counts) == 0) stop("keine Praediktorgen-Counts gesammelt")
    
    cg <- dplyr::bind_rows(gene_counts) %>%
      dplyr::inner_join(dplyr::select(bin_df, ID, ZT), by = "ID") %>%
      dplyr::mutate(
        Group     = factor(substr(Sample_ID, 1, 1), levels = c("M", "P")),
        ZT_Gruppe = factor(ifelse(ZT < 12, "vor_12", "nach_12"), levels = c("vor_12", "nach_12")),
        Gene      = factor(Gene, levels = pred_genes[pred_genes %in% Gene])
      )
    
    ca <- cg %>%
      dplyr::group_by(ID, Group, ZT_Gruppe) %>%
      dplyr::summarise(Counts = sum(Counts), .groups = "drop")
    
    st_a <- ca %>% dplyr::group_by(Group, ZT_Gruppe) %>%
      dplyr::summarise(mean = mean(Counts), .groups = "drop")
    st_g <- cg %>% dplyr::group_by(Gene, Group, ZT_Gruppe) %>%
      dplyr::summarise(mean = mean(Counts), .groups = "drop")
    
    fl_lab <- labeller(Group = group_labels, ZT_Gruppe = zt_group_labels)
    xa_max <- max(1, stats::quantile(ca$Counts, counts_quantil))
    xg_max <- max(1, stats::quantile(cg$Counts, counts_quantil))
    cut_txt <- function(q) paste0("x axis shown up to the ", 100 * counts_quantil,
                                  "% quantile (", q, " counts); dashed line = mean")
    
    pa <- ggplot(ca, aes(Counts, fill = Group)) +
      geom_bar(width = 0.8) +
      geom_vline(data = st_a, aes(xintercept = mean, linetype = "Mean"), colour = "black") +
      facet_grid(~ Group + ZT_Gruppe, labeller = fl_lab) +
      scale_fill_manual(values = group_colours, labels = group_labels, name = "Group") +
      scale_linetype_manual(values = c("Mean" = "dashed"), name = NULL) +
      coord_cartesian(xlim = c(-0.5, xa_max + 0.5)) +
      labs(title = ttl(paste0("Counts per spot - sum over all predictor genes (", M, ")")),
           subtitle = ttl(cut_txt(xa_max)),
           x = "Counts per spot (sum of predictor genes)", y = "Number of spots") +
      th()
    speichern(pa, "10a_counts_praediktorgene_summe", height = W * 0.45)
    
    pg <- ggplot(cg, aes(Counts, fill = Group)) +
      geom_bar(width = 0.8) +
      geom_vline(data = st_g, aes(xintercept = mean, linetype = "Mean"), colour = "black") +
      geom_text(data = st_g, aes(x = Inf, y = Inf, label = sprintf("mean = %.3f", mean)),
                inherit.aes = FALSE, hjust = 1.1, vjust = 1.5, size = (base_size - 2) / ggplot2::.pt) +
      facet_grid(Gene ~ Group + ZT_Gruppe, scales = "free_y", labeller = fl_lab) +
      scale_fill_manual(values = group_colours, labels = group_labels, name = "Group") +
      scale_linetype_manual(values = c("Mean" = "dashed"), name = NULL) +
      coord_cartesian(xlim = c(-0.5, xg_max + 0.5)) +
      labs(title = ttl(paste0("Counts per spot - per predictor gene (", M, ")")),
           subtitle = ttl(cut_txt(xg_max)),
           x = "Counts per spot (single gene)", y = "Number of spots") +
      th() +
      theme(strip.text.y = element_text(angle = 0))
    speichern(pg, "10b_counts_praediktorgene_je_gen",
              height = min(23, 2 * nlevels(cg$Gene) + 3))
  })
  
  
  ## ================================================================
  ## 11 Bootstrap mit Zuruecklegen auf Pseudobulk-Ebene:
  ##    pro Sample und Zone n = Zonengroesse Spots MIT Zuruecklegen ziehen
  ##    -> Counts summieren (Pseudobulk) -> tauFisher-Vorhersage, R_boot mal.
  ##    Punkt = zirkulaeres Mittel, Balken = 95%-Intervall,
  ##    Kreuz = Vorhersage aus allen Spots der Zone.
  ## ================================================================
  if ("bootstrap" %in% plots) versuche("bootstrap", {
    if (length(boot_dat) == 0) stop("keine Spot-Daten gesammelt (Modell oder Visium-Objekte fehlen)")
    if (!requireNamespace("nnet", quietly = TRUE)) stop("Paket 'nnet' ist nicht installiert")
    f_test <- .tf_fun("test_tauFisher")
    
    # test_tauFisher() verliert bei nur EINER Spalte intern die Matrix-Dimension
    # -> Einzelspalte doppeln, 1. Ergebnis nehmen
    predict_pb <- function(Mx) {
      Mx <- as.matrix(Mx)
      einzeln <- ncol(Mx) == 1
      if (einzeln) Mx <- cbind(Mx, Mx)
      if (is.null(colnames(Mx))) colnames(Mx) <- paste0("pb_", seq_len(ncol(Mx)))
      res <- f_test(chosen_genes = modell$Genes, test_df = Mx, test_time = rep(0, ncol(Mx)),
                    pca_embedding = modell$PCA_embedding, trained_model = modell$Model)
      zt <- round(res$predicted) %% 24
      if (einzeln) zt[1] else zt
    }
    circ_diff <- function(a, b) ((a - b + 12) %% 24) - 12
    circ_summary <- function(zt) {
      a  <- 2 * pi * zt / 24
      C  <- mean(cos(a)); S <- mean(sin(a))
      mu <- (atan2(S, C) * 24 / (2 * pi)) %% 24
      q  <- stats::quantile(circ_diff(zt, mu), c(0.025, 0.975), names = FALSE)
      data.frame(circ_mean = mu, lo = mu + q[1], hi = mu + q[2])
    }
    
    set.seed(123)
    summ <- list()
    for (s in names(boot_dat)) {
      X <- boot_dat[[s]]$X; zone <- boot_dat[[s]]$zone
      zones <- sort(unique(zone))
      # Referenz: alle Spots je Zone summiert, alle Zonen in einem Aufruf
      full <- do.call(cbind, lapply(zones, function(z) Matrix::rowSums(X[, zone == z, drop = FALSE])))
      rownames(full) <- rownames(X)
      colnames(full) <- paste0("Zone_", zones)
      ref_zt <- setNames(predict_pb(full), colnames(full))
      for (z in zones) {
        Xz <- X[, zone == z, drop = FALSE]
        N  <- ncol(Xz)
        ref <- ref_zt[[paste0("Zone_", z)]]
        idx <- sample.int(N, N * R_boot, replace = TRUE)
        Wm  <- Matrix::sparseMatrix(i = idx, j = rep(seq_len(R_boot), each = N), x = 1,
                                    dims = c(N, R_boot))
        zt  <- predict_pb(Xz %*% Wm)
        summ[[length(summ) + 1]] <- cbind(
          data.frame(Sample_ID = s, Zone = paste0("Zone_", z), n_spots = N, ref_all_spots = ref),
          circ_summary(zt))
      }
      message("  Bootstrap ", s, " fertig")
    }
    
    summ <- dplyr::bind_rows(summ) %>%
      add_condition() %>%
      dplyr::mutate(ref_unwrapped = circ_mean + circ_diff(ref_all_spots, circ_mean),
                    Sample_ID = order_samples(Sample_ID),
                    Zone = factor(Zone, levels = paste0("Zone_", 1:8)))
    saveRDS(summ, file.path(out_dir, "11_bootstrap_pseudobulk_zonen_summary.rds"))
    
    # Bulk-CT je Sample als Linie; y-Achse ist "ausgerollt" (Labels %% 24) ->
    # Bulk-CT auf die Lage der Zonenwerte des Samples bringen
    bulk_bl <- NULL
    if (!is.null(bulk_lookup)) {
      bulk_bl <- summ %>%
        dplyr::group_by(Sample_ID) %>%
        dplyr::summarise(zentrum = .circ_mean(circ_mean), .groups = "drop") %>%
        dplyr::mutate(bulk = unname(bulk_lookup[as.character(Sample_ID)]),
                      y = zentrum + circ_diff(bulk, zentrum)) %>%
        dplyr::filter(!is.na(y))
    }
    
    p11 <- ggplot(summ, aes(Zone, circ_mean, colour = Condition, group = Sample_ID)) +
      geom_line(alpha = 0.7) +
      geom_pointrange(aes(ymin = lo, ymax = hi), size = 0.2) +
      geom_point(aes(y = ref_unwrapped, shape = "All spots of the zone"),
                 colour = "black", size = 1.6) +
      facet_wrap(~ Sample_ID, ncol = 4) +
      scale_colour_manual(values = condition_colours, name = "Condition", drop = TRUE) +
      scale_shape_manual(values = c("All spots of the zone" = 4), name = NULL) +
      scale_x_discrete(labels = 1:8) +
      scale_y_continuous(breaks = seq(-24, 48, 4), labels = function(x) x %% 24) +
      labs(title = ttl(paste0("Predicted circadian phase (", M, ") per zone: pseudobulk bootstrap")),
           subtitle = ttl(paste0("Spots of each zone resampled with replacement (", R_boot,
                                 " replicates)\nPoint = circular mean, bar = 95% interval")),
           x = "Zone (1 = central, 8 = portal)", y = "Predicted CT (pseudobulk)") +
      th() + theme(legend.box = "horizontal", legend.box.just = "top") +
      guides(colour = guide_legend(nrow = 3, order = 1),
             shape  = guide_legend(order = 2))
    if (!is.null(bulk_bl) && nrow(bulk_bl) > 0) {
      p11 <- p11 +
        geom_hline(data = bulk_bl, aes(yintercept = y, linetype = "Sample CT"),
                   colour = "grey40", linewidth = 0.4) +
        scale_linetype_manual(values = c("Sample CT" = "dashed"), name = NULL) +
        guides(linetype = guide_legend(order = 3))
    }
    speichern(p11, "11_bootstrap_pseudobulk_zonen", height = W * 1.1)
  })
  
  
  ## ================================================================
  ## 12 Visium-Zonen neben Zonen-Prediction (2 Seiten a 8 Samples)
  ## ================================================================
  if ("zonen_vs_prediction" %in% plots) versuche("zonen_vs_prediction", {
    if (length(pair_plots) == 0) stop("keine Zonen/Prediction-Paare erzeugt")
    n_s  <- length(pair_plots)
    half <- ceiling(n_s / 2)
    seiten <- list(names(pair_plots)[1:half])
    if (n_s > half) seiten[[2]] <- names(pair_plots)[(half + 1):n_s]
    
    for (i in seq_along(seiten)) {
      p <- wrap_plots(pair_plots[seiten[[i]]], ncol = 2, guides = "collect") &
        theme(legend.position = "bottom",
              legend.box = "vertical",
              legend.spacing.y = unit(0.1, "cm"),
              legend.key.width = unit(1.2, "cm"), legend.key.height = unit(0.25, "cm")) &
        guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5))
      if (titel) p <- p + plot_annotation(
        title = paste0("Liver zones and predicted CT per zone (", M, ") - part ", i, "/", length(seiten)),
        theme = theme(plot.title = element_text(face = "bold", size = base_size + 1)))
      speichern(p, paste0("12_zonen_vs_prediction_seite", i), height = W * 1.3)
    }
  })
  
  
  ## ================================================================
  ## 13 Violinplots (Default CRY2, PER3): vor vs. nach CT 12, M | P
  ## ================================================================
  if ("violin" %in% plots) versuche("violin", {
    if (length(expr_list) == 0) stop("keine Expressionsdaten gesammelt")
    
    vd <- dplyr::bind_rows(expr_list) %>%
      dplyr::filter(is.finite(Expression)) %>%
      dplyr::inner_join(dplyr::select(bin_df, ID, ZT), by = "ID") %>%
      dplyr::mutate(ZT_Gruppe = factor(ifelse(ZT < 12, "vor_12", "nach_12"),
                                       levels = c("vor_12", "nach_12"))) %>%
      add_condition() %>%
      dplyr::mutate(Sample_ID = order_samples(Sample_ID))
    if (violin_nur_positive) vd <- dplyr::filter(vd, Expression > 0)
    
    y_lab <- if (violin_nur_positive) "Expression (log-normalised, spots > 0)" else
      "Expression (log-normalised)"
    m_ids <- sample_level_order[startsWith(sample_level_order, "M")]
    p_ids <- sample_level_order[startsWith(sample_level_order, "P")]
    
    for (g in unique(vd$Gene)) {
      d_g <- dplyr::filter(vd, Gene == g)
      y_min <- min(d_g$Expression); y_top <- max(d_g$Expression)
      if (y_top == y_min) y_top <- y_min + 1
      y_up <- y_top + 0.15 * (y_top - y_min)
      
      n_g <- d_g %>%
        dplyr::count(Sample_ID, ZT_Gruppe) %>%
        dplyr::mutate(y = y_top)
      
      mache_violin <- function(ids, gruppe) {
        d <- dplyr::filter(d_g, Sample_ID %in% ids)
        ggplot(d, aes(ZT_Gruppe, Expression, fill = Condition, alpha = ZT_Gruppe)) +
          # show.legend = TRUE: sonst leere Legendenkaestchen fuer fehlende Levels (ggplot2 >= 3.5)
          geom_violin(scale = "width", trim = TRUE, colour = "grey30", linewidth = 0.3,
                      show.legend = TRUE) +
          stat_summary(fun = median, geom = "point", size = 0.9, colour = "black", alpha = 1,
                       show.legend = FALSE) +
          geom_text(data = dplyr::filter(n_g, Sample_ID %in% ids),
                    aes(ZT_Gruppe, y, label = paste0("n=", n)),
                    inherit.aes = FALSE, vjust = -0.4, size = (base_size - 2.5) / ggplot2::.pt) +
          facet_wrap(~ Sample_ID, ncol = 4, labeller = as_labeller(label_with_condition_2z)) +
          # gleiche Levels in M- und P-Plot -> patchwork fasst die Legenden zusammen
          scale_fill_manual(values = condition_colours, name = "Condition",
                            limits = condition_levels, drop = FALSE) +
          scale_alpha_manual(values = zt_group_alpha, labels = zt_group_labels, name = "CT group",
                             guide = guide_legend(override.aes = list(fill = "grey30"))) +
          scale_x_discrete(labels = zt_group_labels) +
          coord_cartesian(ylim = c(y_min, y_up)) +
          labs(title = group_labels[[gruppe]], x = NULL, y = y_lab) +
          th() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
      }
      
      # Condition-Legende nur EINMAL zeichnen (aus dem P-Plot)
      p <- ((mache_violin(m_ids, "M") + guides(fill = "none")) /
              (mache_violin(p_ids, "P") +
                 guides(fill = guide_legend(nrow = 3, override.aes = list(alpha = 1), order = 1)))) +
        plot_layout(guides = "collect") &
        theme(legend.position = "bottom",
              legend.box = "horizontal",
              legend.box.just = "top",
              legend.key.size = unit(0.35, "cm")) &
        guides(alpha = guide_legend(ncol = 1, override.aes = list(fill = "grey30"), order = 2))
      if (titel) p <- p + plot_annotation(
        title = paste0(g, " - expression per spot, before vs. after CT 12 (", M, ")"),
        theme = theme(plot.title = element_text(face = "bold", size = base_size + 1)))
      speichern(p, paste0("13_violin_", g), height = W * 1.25)
    }
  })
  
  
  ## ================================================================
  ## 14 Delta CT relativ zu Zone 1, t-Test & Wilcoxon gegen 0, M | P
  ## ================================================================
  if ("test_zone1" %in% plots) versuche("test_zone1", {
    if (is.null(zonen_df)) stop("zonen_prediction fehlt")
    
    dz <- zonen_df %>%
      dplyr::mutate(Group = factor(substr(Sample_ID, 1, 1), levels = c("M", "P"))) %>%
      dplyr::group_by(Sample_ID) %>%
      dplyr::mutate(ZT_diff = ZT - ZT[Zone == "Zone_1"]) %>%
      dplyr::ungroup()
    
    sm <- dz %>%
      dplyr::group_by(Zone, Group) %>%
      dplyr::summarise(mean_diff = mean(ZT_diff), sd_diff = sd(ZT_diff), .groups = "drop")
    
    tests <- dz %>%
      dplyr::filter(Zone != "Zone_1") %>%
      dplyr::group_by(Zone, Group) %>%
      dplyr::summarise(
        p_ttest  = tryCatch(t.test(ZT_diff, mu = 0)$p.value, error = function(e) NA_real_),
        p_wilcox = tryCatch(wilcox.test(ZT_diff, mu = 0, exact = FALSE)$p.value, error = function(e) NA_real_),
        .groups = "drop"
      ) %>%
      dplyr::group_by(Group) %>%
      dplyr::mutate(p_ttest_adj  = p.adjust(p_ttest, method = "BH"),
                    p_wilcox_adj = p.adjust(p_wilcox, method = "BH"),
                    sig_t = .sig_stars(p_ttest_adj),
                    sig_w = .sig_stars(p_wilcox_adj)) %>%
      dplyr::ungroup()
    print(as.data.frame(tests))
    
    grp_lab  <- as_labeller(group_labels)
    test_sub <- ttl("Each zone tested against 0 within group; BH-adjusted. t = one-sample t-test, W = Wilcoxon signed-rank")
    
    off_a <- 0.07 * diff(range(c(sm$mean_diff - sm$sd_diff, sm$mean_diff + sm$sd_diff), na.rm = TRUE))
    lab_a <- tests %>%
      dplyr::left_join(sm, by = c("Zone", "Group")) %>%
      dplyr::mutate(y_t = mean_diff + sd_diff + off_a, y_w = y_t + off_a)
    
    p14a <- ggplot(sm, aes(Zone, mean_diff, colour = Group, group = Group)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
      geom_line(linewidth = 0.8) +
      geom_point(size = 1.8) +
      geom_errorbar(aes(ymin = mean_diff - sd_diff, ymax = mean_diff + sd_diff), width = 0.2) +
      geom_text(data = lab_a, aes(Zone, y_t, label = paste0("t: ", sig_t)),
                inherit.aes = FALSE, size = (base_size - 2) / ggplot2::.pt) +
      geom_text(data = lab_a, aes(Zone, y_w, label = paste0("W: ", sig_w)),
                inherit.aes = FALSE, size = (base_size - 2) / ggplot2::.pt) +
      facet_wrap(~ Group, labeller = grp_lab) +
      scale_x_discrete(labels = 1:8) +
      scale_colour_manual(values = group_colours, guide = "none") +
      scale_y_continuous(labels = .hour_labels) +
      labs(title = ttl(paste0("\u0394 CT (", M, ") relative to zone 1: mean \u00b1 SD")),
           subtitle = test_sub,
           x = "Zone (1 = central, 8 = portal)", y = "\u0394 CT relative to zone 1") +
      th()
    speichern(p14a, "14_delta_zone1_tests", height = W * 0.55)
  })
  
  message("\nFertig (", M, "). Abbildungen liegen in: ", out_dir)
  invisible(abbildungen)   # ggplot-Objekte zum Nachjustieren
}


## ================================================================
## Beispielaufruf (Pfade aus config.R)
## ================================================================
# abb <- erstelle_thesis_plots(
#   out_dir         = file.path(abbildungen_dir, "JTK"),
#   pfad_pred       = path_to_save_pred,
#   pfad_modell_dir = path_to_save_tauFisher_model,
#   visium_dir      = path_to_visium_data,
#   methode         = "JTK",          # oder "LS", "LS_full"
#   formats         = c("png", "pdf"),
#   textbreite_cm   = 16
# )
#
# nur einzelne Plots neu erzeugen (schnell, ohne Seurat-Objekte):
# erstelle_thesis_plots(..., plots = c("sem_sd", "flatness", "kreis_bulk", "zonen_linien"))