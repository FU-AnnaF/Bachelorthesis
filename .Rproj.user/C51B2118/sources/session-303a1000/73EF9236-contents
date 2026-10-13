# =============================================================
# ALL Code was sortet and partly commented with Claude Opus 5.5 (Pro plan) on 12.10.2026
#
#
# run_all.R – führt die komplette Pipeline aus
# Öffne zuerst "Bachelorarbeit Pipeline.Rproj", dann dieses Skript ausführen.
# =============================================================

# here wird schon für config.R gebraucht, also vor load_packages_pipeline()
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
library(here)

# ---- Ordner, in dem ALLE R-Skripte liegen
skript_dir <- here("skripte")

# ---- Einstellungen laden
source(file.path(skript_dir, "config.R"), encoding = "UTF-8")

# ---- Funktionsskripte laden (alle in skripte/)
funktions_skripte <- c(
  "load_packages_pipeline.R",
  "visium_data_vorbereitung.R",
  "zonierung_pipeline.R",
  "tauFisher_pipeline.R",
  "plots_pipeline.R"
)
for (f in funktions_skripte) {
  source(file.path(skript_dir, f), encoding = "UTF-8")
}

load_packages_pipeline()

# ---- 1. Entpackte Visium-Daten sortieren (nur einmalig)
if (schritte$vorbereitung) {
  visium_data_vorbereitung(visium_dir = path_to_visium_data)
}

# ---- 2. Zonen berechnen, *_zoniert.rds speichern
if (schritte$zonierung) {
  zonierung_pipeline(path_to_visium = path_to_visium_data,
                     save_plots_to  = save_plots_to)
}

# ---- 3. tauFisher mit Talamanca-Trainingsdaten
if (schritte$taufisher) {
  tauFisher_pipeline(path_to_xlsx       = path_to_xlsx_file,
                     path_to_gct        = path_to_gct_file,
                     path_to_visium     = path_to_visium_data,
                     path_to_save_model = path_to_save_tauFisher_model,
                     path_to_save_pred  = path_to_save_pred,
                     save_plot_file_to  = save_plots_to,
                     bulk = bulk, zonen = zonen, bin = bin,
                     training_nur_170 = training_nur_170)
}

# ---- 4. Thesis-Abbildungen, je Methode ein eigener Ordner
if (schritte$plots) {
  abb <- list()
  for (m in plot_methoden) {
    abb[[m]] <- erstelle_thesis_plots(
      out_dir         = file.path(abbildungen_dir, m),
      pfad_pred       = path_to_save_pred,
      pfad_modell_dir = path_to_save_tauFisher_model,
      visium_dir      = path_to_visium_data,
      methode         = m,
      formats         = c("png", "pdf"),
      textbreite_cm   = 16
      # plots = c("zonen_linien", "spots_sample")  # nur einzelne Plots neu erzeugen
    )
  }
}

message("FINISHED :)  Ergebnisse liegen in: ", ergebnis_dir)