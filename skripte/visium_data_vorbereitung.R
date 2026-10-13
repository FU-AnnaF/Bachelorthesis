# =============================================================
# visium_data_vorbereitung.R
# Speicherort: Bachelorarbeit Pipeline/skripte/
#
# Voraussetzung: ALLE ZIP-Dateien sind bereits entpackt und liegen
# (in beliebigen Unterordnern) in visium_dir, z.B. daten/Visium.
# Das Skript entpackt nichts, sondern räumt nur auf und sortiert:
#   1. Non_Human- und Loupe-Dateien entfernen
#   2. Dateien pro Probe in <Probe>_Visium/ sortieren
#   3. Dateien im Probenordner in die Seurat-Struktur (spatial/) bringen
#   4. leere Restordner entfernen
#   5. tissue_lowres_image aus dem hires-Bild erzeugen
# ZIP-Dateien selbst werden nicht angefasst.
# =============================================================

visium_samples <- c("M1","M2","M3","M4","M5","M6","M7","M8",
                    "P2","P3","P6","P7","P14","P17","P18","P21")


#### Dateien nach Probe sortieren ####
sort_visium_files <- function(visium_dir, samples = visium_samples) {
  
  # laengste Prefixe zuerst pruefen, damit z.B. "P21" vor "P2" zugeordnet wird
  samples <- samples[order(-nchar(samples))]
  
  all_files <- list.files(visium_dir, full.names = TRUE, recursive = TRUE)
  
  # ZIP-Dateien nicht verschieben
  all_files <- all_files[tolower(tools::file_ext(all_files)) != "zip"]
  
  assigned <- rep(FALSE, length(all_files))
  
  for (s in samples) {
    target_dir <- file.path(visium_dir, paste0(s, "_Visium"))
    
    # Nur Dateien, deren Name mit der Probe beginnt und die noch nicht zugeordnet sind
    hit <- !assigned & grepl(paste0("^", s), basename(all_files))
    matches <- all_files[hit]
    
    if (length(matches) == 0) next
    
    dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Dateien, die schon im richtigen Ordner (oder dessen spatial/) liegen, bleiben
    schon_da <- normalizePath(dirname(matches), winslash = "/") %in%
      normalizePath(c(target_dir, file.path(target_dir, "spatial")), winslash = "/", mustWork = FALSE)
    
    file.rename(matches[!schon_da], file.path(target_dir, basename(matches[!schon_da])))
    
    assigned[hit] <- TRUE
  }
  
  invisible(NULL)
}


#### Dateien im Probenordner in Seurat-Struktur bringen ####
organize_visium_sample_folder <- function(sample_dir) {
  
  sample_name <- sub("_Visium$", "", basename(sample_dir))
  spatial_dir <- file.path(sample_dir, "spatial")
  dir.create(spatial_dir, showWarnings = FALSE)
  
  files <- list.files(sample_dir, full.names = TRUE)
  files <- files[!dir.exists(files)]
  
  for (f in files) {
    fname  <- basename(f)
    ext    <- tools::file_ext(fname)
    base   <- tools::file_path_sans_ext(fname)
    suffix <- sub(paste0("^", sample_name, "_"), "", base)
    
    ext_suffix <- if (nzchar(ext)) paste0(".", ext) else ""
    
    if (tolower(ext) == "cloupe") {
      # .cloupe-Dateien (Loupe Browser) werden nicht gebraucht
      file.remove(f)
      
    } else if (suffix == "tissue_positions_list") {
      file.rename(f, file.path(spatial_dir, fname))
      
    } else if (suffix %in% c("counts_ALL", "counts_UTT",
                             "scalefactors_json", "tissue_hires_image")) {
      file.rename(f, file.path(spatial_dir, paste0(suffix, ext_suffix)))
    }
    # alles andere (z.B. *_filtered_feature_bc_matrix.h5, *_zoniert.rds) bleibt liegen
  }
  
  invisible(NULL)
}

organize_all_visium_samples <- function(visium_dir, samples = visium_samples) {
  for (s in samples) {
    sample_dir <- file.path(visium_dir, paste0(s, "_Visium"))
    if (dir.exists(sample_dir)) {
      message("Sortiere: ", s, "_Visium")
      organize_visium_sample_folder(sample_dir)
    } else {
      warning(s, "_Visium existiert nicht – wurden alle ZIPs entpackt?", call. = FALSE)
    }
  }
  invisible(NULL)
}


#### Leere Ordner entfernen (z.B. die entpackten Human_...-Ordner) ####
remove_empty_dirs <- function(visium_dir) {
  dirs <- list.dirs(visium_dir, full.names = TRUE, recursive = TRUE)
  dirs <- setdiff(dirs, visium_dir)
  # tiefste Ordner zuerst, damit verschachtelte leere Ordner mit verschwinden
  dirs <- dirs[order(-nchar(dirs))]
  for (d in dirs) {
    if (length(list.files(d, all.files = TRUE, no.. = TRUE)) == 0) unlink(d, recursive = TRUE)
  }
  invisible(NULL)
}


#### Lowres-Bilder aus Hires-Bildern erzeugen ####
create_lowres_images <- function(visium_dir, lowres_max_dim = 600) {
  
  hires_files <- list.files(
    visium_dir,
    pattern = "^tissue_hires_image\\.(png|jpg|jpeg)$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  for (f in hires_files) {
    lowres_path <- file.path(dirname(f), paste0("tissue_lowres_image.", tools::file_ext(f)))
    
    if (file.exists(lowres_path)) {
      message(basename(dirname(dirname(f))), ": tissue_lowres_image existiert bereits, ueberspringe.")
      next
    }
    
    img <- magick::image_read(f)
    img <- magick::image_resize(img, paste0(lowres_max_dim, "x", lowres_max_dim))
    magick::image_write(img, path = lowres_path)
    message(basename(dirname(dirname(f))), ": tissue_lowres_image erstellt.")
  }
  
  invisible(NULL)
}


#### Hauptfunktion ####
# visium_dir: Ordner mit den bereits entpackten Visium-Daten (aus config.R)
visium_data_vorbereitung <- function(visium_dir = path_to_visium_data) {
  
  if (!dir.exists(visium_dir) || length(list.files(visium_dir)) == 0) {
    stop("Ordner '", visium_dir, "' fehlt oder ist leer.\n",
         "Bitte die Visium-ZIP-Dateien zuerst dorthin entpacken.")
  }
  
  # Non_Human-Dateien und -Ordner löschen (Groß-/Kleinschreibung egal)
  non_human <- list.files(visium_dir, pattern = "^Non_Human", full.names = TRUE,
                          include.dirs = TRUE, ignore.case = TRUE)
  non_human <- non_human[tolower(tools::file_ext(non_human)) != "zip"]
  unlink(non_human, recursive = TRUE)
  
  # Human_Loupe_files-Ordner löschen
  unlink(file.path(visium_dir, "Human_Loupe_files"), recursive = TRUE)
  
  # Dateien nach Probe sortieren und pro Probenordner strukturieren
  sort_visium_files(visium_dir)
  organize_all_visium_samples(visium_dir)
  
  # leere Restordner der entpackten Archive entfernen
  remove_empty_dirs(visium_dir)
  
  # Lowres-Bilder aus Hires-Bildern erzeugen
  create_lowres_images(visium_dir)
  
  message("Vorbereitung abgeschlossen.")
  invisible(visium_dir)
}