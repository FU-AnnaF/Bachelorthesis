# =============================================================
# config.R – ALLE Einstellungen an einer Stelle
# Speicherort: Bachelorarbeit Pipeline/skripte/config.R
# Pfade sind relativ zum Projektordner "Bachelorarbeit Pipeline"
# (dort, wo die .Rproj liegt) – nicht relativ zu skripte/
# =============================================================

library(here)

# ---- Welche Schritte sollen laufen?
schritte <- list(
  vorbereitung = FALSE,  # entpackte Visium-Daten sortieren (nur einmal nötig)
  zonierung    = TRUE,
  taufisher    = TRUE,
  plots        = TRUE
)

# ---- Wofür soll tauFisher Vorhersagen machen?
bulk   <- TRUE
zonen  <- TRUE
bin    <- TRUE

# ---- Eingabedaten (werden nie verändert)
path_to_xlsx_file  <- here("daten", "raw", "science.add0846_table_s2.xlsx")
path_to_gct_file   <- here("daten", "raw", "gene_reads_adult_gtex_v11_liver.gct")

# ---- Visium-Daten: hierhin die ZIPs entpacken; danach sortiert
#      die Vorbereitung sie in <Probe>_Visium/ (inkl. *_zoniert.rds)
path_to_visium_data <- here("daten", "Visium")

# ---- tauFisher-Training für JTK und LS
# FALSE = wie im Original: Genauswahl auf 170 Spalten, Training auf allen
# TRUE  = Genauswahl UND Training nur auf den ersten 170 Spalten
# Nach dem Umstellen die Modelle in ergebnisse/modelle/modelle/ löschen!
training_nur_170 <- FALSE

# ---- Ergebnisse
#   ergebnisse/
#   ├── modelle/
#   │   ├── modelle/   trainierte tauFisher-Modelle
#   │   └── pred/      Vorhersagen (z.B. zonen_prediction_JTK.rds)
#   └── plots/         alle Plots und Thesis-Abbildungen
ergebnis_dir                 <- here("ergebnisse")
path_to_save_tauFisher_model <- file.path(ergebnis_dir, "modelle", "modelle")
path_to_save_pred            <- file.path(ergebnis_dir, "modelle", "pred")
save_plots_to                <- file.path(ergebnis_dir, "plots")
abbildungen_dir              <- save_plots_to

# Für welche Modelle sollen Thesis-Abbildungen erzeugt werden?
# Standard: nur JTK. Weitere möglich: c("JTK", "LS", "LS_full")
# (je Methode ein Unterordner in ergebnisse/plots/)
plot_methoden <- "JTK"

# ---- Ordner anlegen, falls sie fehlen
for (d in c(path_to_visium_data, path_to_save_tauFisher_model,
            path_to_save_pred, save_plots_to)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}