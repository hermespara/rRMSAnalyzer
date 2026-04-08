#' Plot 3D rRNA structure with annotated sites
#'
#' @description
#' This function visualizes the 3D structure of rRNA or the ribosome using `r3dmol`.
#' It highlights specifically annotated sites of interest. Optimized for large structures
#' like the human ribosome.
#'
#' @details
#' The function uses `r3dmol` to fetch a PDB structure and highlight residues corresponding
#' to the annotated sites in the `SummarizedExperiment` object. Mapping between
#' RNA names in the object and chains in the PDB structure is handled via `chain_mapping`.
#'
#' Default mappings for human 80S ribosome (PDB 4UG0):
#' - '28S' mapped to chain 'A5' (LSU rRNA)
#' - '5S' mapped to chain 'A8' (LSU 5S)
#' - '5.8S' mapped to chain 'A6' (LSU 5.8S)
#' - '18S' mapped to chain 'AA' (SSU rRNA)
#'
#' Note: Chain IDs may vary depending on PDB structure. Use `view_pdb_chains()` to inspect.
#'
#' @param ribo A SummarizedExperiment object.
#' @param pdb_id A character string specifying the PDB ID (default: "4UG0").
#' @param chain_mapping A named list or vector mapping RNA names to PDB chain IDs.
#'   Example: `c("28S" = "A5", "18S" = "AA")`. If NULL, uses defaults for 4UG0.
#' @param site_color Color for the highlighted sites (default: "red").
#' @param site_style Style for the highlighted sites: "sphere", "stick", or "cartoon" (default: "sphere").
#' @param label_sites Logical, whether to add labels to highlighted sites (default: TRUE).
#' @param cartoon_quality Quality of the cartoon representation (1-10, default: 3).
#'   Lower values improve performance. Use 1-3 for ribosome.
#' @param show_proteins Logical, whether to show ribosomal proteins (default: FALSE).
#'   Setting to FALSE dramatically improves performance for ribosome visualization.
#' @param background_color Background color (default: "white").
#' @param rna_color Color scheme for rRNA: "spectrum", "chain", "secondary", or specific color (default: "spectrum").
#' @param rna_opacity Opacity of RNA cartoon (0-1, default: 0.7).
#' @param sphere_radius Radius of spheres for highlighted sites (default: 2.5).
#' @param use_mmcif Logical, use mmCIF format instead of PDB (default: TRUE).
#'   mmCIF handles large structures better.
#' @param viewer_width Width of the viewer in pixels (default: 800).
#' @param viewer_height Height of the viewer in pixels (default: 600).
#' @param ... Additional arguments passed to `r3dmol::r3dmol()`.
#'
#' @return An `r3dmol` htmlwidget object.
#' @export
#'
#' @examples
#' \dontrun{
#' data("ribo_toy")
#' # View with default settings (rRNA cartoon, sites as spheres)
#' plot_3d_structure(ribo_toy)
#'
#' # Different color scheme
#' plot_3d_structure(ribo_toy, rna_color = "lightblue")
#'
#' # Show proteins (slower)
#' plot_3d_structure(ribo_toy, show_proteins = TRUE, cartoon_quality = 2)
#' }
plot_3d_structure <- function(ribo,
                              pdb_id = "4UG0",
                              chain_mapping = NULL,
                              site_color = "red",
                              site_style = "sphere",
                              label_sites = FALSE,
                              cartoon_quality = 3,
                              show_proteins = FALSE,
                              background_color = "white",
                              rna_color = "#D1D5DB",
                              rna_opacity = 0.55,
                              sphere_radius = 1.1,
                              use_mmcif = TRUE,
                              viewer_width = 900,
                              viewer_height = 700,
                              ...) {
    if (!requireNamespace("r3dmol", quietly = TRUE)) {
        cli::cli_abort("Package {.pkg r3dmol} is required.")
    }
    if (!requireNamespace("curl", quietly = TRUE)) {
        cli::cli_abort("Package {.pkg curl} is required.")
    }

    check_is_se(ribo)
    check_type(pdb_id, "character", "pdb_id", length = 1)
    check_type(site_color, "character", "site_color", length = 1)
    check_type(site_style, "character", "site_style", length = 1)
    check_type(label_sites, "logical", "label_sites", length = 1)
    check_type(cartoon_quality, "numeric", "cartoon_quality", length = 1)
    check_type(show_proteins, "logical", "show_proteins", length = 1)
    check_type(background_color, "character", "background_color", length = 1)
    check_type(rna_color, "character", "rna_color", length = 1)
    check_type(rna_opacity, "numeric", "rna_opacity", length = 1)
    check_type(sphere_radius, "numeric", "sphere_radius", length = 1)
    check_type(use_mmcif, "logical", "use_mmcif", length = 1)
    check_type(viewer_width, "numeric", "viewer_width", length = 1)
    check_type(viewer_height, "numeric", "viewer_height", length = 1)

    pdb_id <- toupper(pdb_id)
    if (pdb_id != "4UG0") {
        cli::cli_abort(
            "{.fn plot_3d_structure} currently supports only {.val 4UG0} to guarantee robust rRNA mapping."
        )
    }
    if (!isTRUE(use_mmcif)) {
        cli::cli_abort("{.val 4UG0} is supported only through mmCIF. Set {.arg use_mmcif = TRUE}.")
    }

    check_in_set(tolower(site_style), c("sphere", "stick", "cartoon"), "site_style")
    if (rna_opacity < 0 || rna_opacity > 1) {
        cli::cli_abort("{.arg rna_opacity} must be in [0, 1].")
    }

    annots <- get_annotation(ribo)
    if (nrow(annots) == 0) {
        cli::cli_abort("No annotated sites found in {.arg ribo}. Use {.fn annotate_site} first.")
    }

    if (is.null(chain_mapping)) {
        chain_mapping <- c("28S" = "A5", "5.8S" = "A6", "5S" = "A8", "18S" = "AA")
    }
    check_type(chain_mapping, "character", "chain_mapping")
    if (is.null(names(chain_mapping)) || any(names(chain_mapping) == "")) {
        cli::cli_abort("{.arg chain_mapping} must be a named character vector, e.g. c(\"18S\" = \"AA\").")
    }

    cif_url <- paste0("https://files.rcsb.org/download/", pdb_id, ".cif")
    cif_text <- tryCatch(
        {
            raw <- curl::curl_fetch_memory(cif_url)$content
            rawToChar(raw)
        },
        error = function(e) {
            cli::cli_abort("Failed to download {.val {pdb_id}} mmCIF from RCSB: {e$message}")
        }
    )

    mol <- r3dmol::r3dmol(
        viewer_spec = r3dmol::m_viewer_spec(
            cartoonQuality = cartoon_quality,
            backgroundColor = background_color,
            antialias = FALSE
        ),
        width = viewer_width,
        height = viewer_height,
        ...
    ) %>%
        r3dmol::m_add_model(data = cif_text, format = "cif")

    rna_chains <- unique(unname(chain_mapping))

    # Hide everything by default, then draw only the requested rRNA chains.
    mol <- mol %>% r3dmol::m_set_style(style = r3dmol::m_style_line(hidden = TRUE))

    if (isTRUE(show_proteins)) {
        mol <- mol %>%
            r3dmol::m_set_style(
                style = r3dmol::m_style_line(color = "#E5E7EB", opacity = 0.15)
            )
    }

    chain_palette <- c("18S" = "#BFDBFE", "28S" = "#FDE68A", "5.8S" = "#DDD6FE", "5S" = "#BBF7D0")
    for (rna_name in names(chain_mapping)) {
        chain_id <- unname(chain_mapping[[rna_name]])
        if (tolower(rna_color) == "spectrum" || tolower(rna_color) == "chain") {
            chain_col <- chain_palette[[rna_name]]
            if (is.null(chain_col)) chain_col <- "#4B5563"
        } else {
            chain_col <- rna_color
        }
        mol <- mol %>%
            r3dmol::m_set_style(
                sel = r3dmol::m_sel(chain = chain_id),
                style = r3dmol::m_style_line(color = chain_col, opacity = rna_opacity)
            )
    }

    .map_rna_to_chain <- function(rna_value, mapping) {
        idx <- which(tolower(names(mapping)) == tolower(as.character(rna_value)))
        if (length(idx) > 0) return(unname(mapping[[idx[1]]]))

        # Fallback: remove common prefixes/suffixes and re-match.
        norm <- tolower(gsub("[^a-z0-9\\.]", "", as.character(rna_value)))
        mapping_norm <- tolower(gsub("[^a-z0-9\\.]", "", names(mapping)))
        idx2 <- which(mapping_norm == norm)
        if (length(idx2) > 0) return(unname(mapping[[idx2[1]]]))

        NA_character_
    }

    .site_style_object <- function(style_name, color, radius) {
        style_name <- tolower(style_name)
        if (style_name == "sphere") return(r3dmol::m_style_sphere(color = color, radius = radius))
        if (style_name == "stick") return(r3dmol::m_style_stick(color = color, radius = 0.22))
        r3dmol::m_style_cartoon(color = color, opacity = 1)
    }

    # Keep only sites we can map to rRNA chains.
    annots$chain <- vapply(annots$rna, .map_rna_to_chain, character(1), mapping = chain_mapping)
    annots <- annots[!is.na(annots$chain) & annots$chain != "", ]

    if (nrow(annots) == 0) {
        cli::cli_abort("None of the annotated RNA names in {.arg ribo} match {.arg chain_mapping}.")
    }
    if (isTRUE(label_sites) && nrow(annots) > 40) {
        cli::cli_warn("More than 40 sites selected; labels disabled to avoid overlap.")
        label_sites <- FALSE
    }

    site_style_obj <- .site_style_object(site_style, site_color, sphere_radius)
    for (i in seq_len(nrow(annots))) {
        sel <- r3dmol::m_sel(chain = annots$chain[i], resi = as.character(annots$rnapos[i]))

        mol <- mol %>%
            r3dmol::m_add_style(sel = sel, style = site_style_obj)

        if (isTRUE(label_sites)) {
            mol <- mol %>%
                r3dmol::m_add_label(
                    text = as.character(annots$site[i]),
                    sel = sel,
                    style = r3dmol::m_style_label(
                        backgroundColor = "white",
                        fontColor = "black",
                        backgroundOpacity = 0.75,
                        borderColor = site_color,
                        borderThickness = 1.5,
                        fontSize = 11
                    )
                )
        }
    }

    mol %>%
        r3dmol::m_rotate(axis = "x", angle = -90) %>%
        r3dmol::m_rotate(axis = "y", angle = 180) %>%
        r3dmol::m_zoom_to(sel = r3dmol::m_sel(chain = rna_chains))
}



#' Helper function to view available chains in a PDB structure
#'
#' @param pdb_id PDB ID to inspect
#' @param use_mmcif Logical, use mmCIF format (default: TRUE, required for large structures)
#' @return Data frame with chain information
#' @export
view_pdb_chains <- function(pdb_id, use_mmcif = TRUE) {
    cli::cli_alert_info("Fetching structure information for {.val {pdb_id}}...")

    if (use_mmcif) {
        # Use mmCIF format for large structures like ribosomes
        if (!requireNamespace("bio3d", quietly = TRUE)) {
            cli::cli_abort("Package {.pkg bio3d} is required. Install with: install.packages('bio3d')")
        }

        # Download CIF file
        cif_url <- paste0("https://files.rcsb.org/download/", pdb_id, ".cif")
        temp_file <- tempfile(fileext = ".cif")

        tryCatch(
            {
                utils::download.file(cif_url, temp_file, quiet = TRUE, mode = "wb")
            },
            error = function(e) {
                cli::cli_abort("Failed to download CIF file for {.val {pdb_id}}: {e$message}")
            }
        )

        # Read CIF file
        pdb <- bio3d::read.cif(temp_file)
        unlink(temp_file)

        # Extract chain information from atom records
        if (!is.null(pdb$atom)) {
            chains <- unique(pdb$atom$chain)

            # Get detailed info per chain
            chain_summary <- lapply(chains, function(ch) {
                chain_atoms <- pdb$atom[pdb$atom$chain == ch, ]

                # Determine if it's protein or nucleic acid
                res_types <- unique(chain_atoms$resid)

                # Common RNA residues
                rna_residues <- c("A", "C", "G", "U", "DA", "DC", "DG", "DT")
                is_rna <- any(res_types %in% rna_residues)

                # Common protein residues
                protein_residues <- c(
                    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU",
                    "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
                    "PRO", "SER", "THR", "TRP", "TYR", "VAL"
                )
                is_protein <- any(res_types %in% protein_residues)

                mol_type <- if (is_rna) "RNA" else if (is_protein) "Protein" else "Other"

                data.frame(
                    chain = ch,
                    n_residues = length(unique(chain_atoms$resno)),
                    n_atoms = nrow(chain_atoms),
                    type = mol_type,
                    first_residue = min(chain_atoms$resno),
                    last_residue = max(chain_atoms$resno),
                    stringsAsFactors = FALSE
                )
            })

            result <- do.call(rbind, chain_summary)
            result <- result[order(result$type, result$chain), ]

            cli::cli_alert_success("Found {nrow(result)} chains in {.val {pdb_id}}")

            # Print summary
            rna_chains <- result[result$type == "RNA", ]
            protein_chains <- result[result$type == "Protein", ]

            if (nrow(rna_chains) > 0) {
                cli::cli_alert_info("RNA chains: {.val {rna_chains$chain}}")
            }
            if (nrow(protein_chains) > 0) {
                cli::cli_alert_info("Protein chains: {nrow(protein_chains)} total")
            }

            return(result)
        } else {
            cli::cli_abort("No atom records found in structure")
        }
    } else {
        # Try standard PDB format (will fail for large structures)
        if (!requireNamespace("bio3d", quietly = TRUE)) {
            cli::cli_abort("Package {.pkg bio3d} is required. Install with: install.packages('bio3d')")
        }

        pdb <- tryCatch(
            {
                bio3d::read.pdb(pdb_id)
            },
            error = function(e) {
                cli::cli_abort("Failed to read PDB file. Large structures require mmCIF format (set use_mmcif=TRUE): {e$message}")
            }
        )

        chain_info <- unique(pdb$atom[, c("chain", "resid", "type")])

        # Summary by chain
        summary <- do.call(rbind, lapply(split(chain_info, chain_info$chain), function(x) {
            data.frame(
                chain = unique(x$chain),
                n_residues = nrow(x),
                type = paste(unique(x$type), collapse = ","),
                stringsAsFactors = FALSE
            )
        }))

        return(summary)
    }
}


#' Quick function to identify rRNA chains in a ribosome structure
#'
#' @param pdb_id PDB ID to inspect
#' @return Named vector of rRNA chain assignments
#' @export
identify_rRNA_chains <- function(pdb_id) {
    chains <- view_pdb_chains(pdb_id, use_mmcif = TRUE)

    # Filter for RNA chains
    rna_chains <- chains[chains$type == "RNA", ]

    if (nrow(rna_chains) == 0) {
        cli::cli_warn("No RNA chains detected in {.val {pdb_id}}")
        return(NULL)
    }

    # For human ribosome, try to identify based on length
    # Approximate lengths:
    # 18S: ~1800 nt
    # 5.8S: ~160 nt
    # 5S: ~120 nt
    # 28S: ~5000 nt

    rna_chains <- rna_chains[order(-rna_chains$n_residues), ]

    assignments <- character(0)

    for (i in 1:nrow(rna_chains)) {
        chain_id <- rna_chains$chain[i]
        n_res <- rna_chains$n_residues[i]

        # Assign based on length
        if (n_res > 4000) {
            name <- "28S"
        } else if (n_res > 1500 && n_res < 2500) {
            name <- "18S"
        } else if (n_res > 140 && n_res < 200) {
            name <- "5.8S"
        } else if (n_res > 100 && n_res < 140) {
            name <- "5S"
        } else {
            name <- paste0("RNA_", chain_id)
        }

        assignments[name] <- chain_id
        cli::cli_alert_info("{name} -> chain {.val {chain_id}} ({n_res} residues)")
    }

    return(assignments)
}
