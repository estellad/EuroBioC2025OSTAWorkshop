getBoxRange <- function(obj){
  data <- spatialCoords(obj) %>% as.data.frame() %>% 
    summarise(xmin = min(pxl_col_in_fullres), xmax = max(pxl_col_in_fullres),
              ymin = min(pxl_row_in_fullres), ymax = max(pxl_row_in_fullres))
  return(data)
}

plotBox <- function(data, col, lty, lwd){
  geom_rect(data = data,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
            color = col, fill = NA, linetype = lty, linewidth = lwd)
}

#' Flip the Y-axis of cell or nucleus segmentations to align with H&E image. 
#'
#' @param sf a sf object read from a `.geojson` file.
#' @param type "POINT" for cell centroid, or "POLYGON" for cell segmentation 
#' mask. Default is "POINT".
#' @param img_height total length along the Y axis of the image. Obtained by 
#' reading in `hires` or `lowre`s `.png` under `/spatial` folder with 
#' `magick::image_read()`.
#' @param scalef scaling factor from `/spatial/scalefactors_json.json` file
#'
#' @returns a sf object with Y-axis of the points or polygons flipped.
#'
#' @examples
#' geo_data_flipped <- flip_sf_Y(sf = geo_data, type = "POLYGON", 
#'                               img_height = 3886, scalef = 0.079)
flip_sf_Y <- function(sf, type = "POINT", img_height, scalef){
  st_geometry(sf) <- st_sfc( 
    lapply(st_geometry(sf), function(geom) {
      coords <- st_coordinates(geom)
      coords[, 2] <- img_height/scalef - coords[, 2] 
      if(type == "POINT"){
        st_point(coords)
      }else{ # type == "POLYGON"
        st_polygon(list(matrix(coords[, 1:2], ncol = 2)))
      }
    }),
    crs = st_crs(sf)
  )
  
  return(sf)
}


#' Reader for Visium HD cell segmented output.
#'
#' @param td path to unzipped Visium HD data download with `/segmented_outputs`
#' @param res "hires" or "lowres" for what .png to read in. 
#'
#' @returns a `SpatialFeatureExperiment` object.
#'
#' @examples
#' vhdsfe <- readVisiumHDCellSeg(td = "~/Desktop", res = "hires")
readVisiumHDCellSeg <- function(td, res){
  if(res == "hires"){
    png_name = "tissue_hires_image.png"
    scalef_name = "tissue_hires_scalef"
  }else{ # res = "lowres"
    png_name = "tissue_lowres_image.png"
    scalef_name = "tissue_lowres_scalef"
  }
  
  ## Coord ---------------------------------------------------------------
  # Read the GeoJSON file of cell segmentation
  geo_data <- st_read(file.path(td, "segmented_outputs", 
                                "cell_segmentations.geojson"))
  scalef <- rjson::fromJSON(file = file.path(td, "segmented_outputs", "spatial", 
                                             "scalefactors_json.json"))
  image <- image_read(file.path(td, "segmented_outputs", 
                                "spatial", png_name))
  info <- image_info(image)
  
  # Get centroid
  st_crs(geo_data) <- NA
  rownames(geo_data) <- geo_data$cell_id
  centroids <- st_centroid(geo_data)
  centroids$cell_id <- as.character(centroids$cell_id)
  
  
  ## Countmat ------------------------------------------------------------
  countmat_file <- file.path(td, "segmented_outputs", 
                             "filtered_feature_cell_matrix.h5")
  vhdcellsce <- DropletUtils::read10xCounts(countmat_file, col.names = TRUE)
  
  colnames(vhdcellsce) <- sub("cellid_0*([0-9]+)-.*", "\\1", vhdcellsce$Barcode)
  
  rownames(vhdcellsce) <- rowData(vhdcellsce)$Symbol
  rownames(vhdcellsce) <- make.unique(rownames(vhdcellsce))
  
  
  ## Matching coordinates with count matrix ------------------------------
  # Subset to common cells
  common_cells <- intersect(centroids$cell_id, colnames(vhdcellsce))
  
  geo_data <- geo_data[geo_data$cell_id %in% common_cells, ]
  centroids <- centroids[centroids$cell_id %in% common_cells, ]
  vhdcellsce <- vhdcellsce[, colnames(vhdcellsce) %in% common_cells]
  
  # # Sanity check on ordering of cell ids
  # all(as.character(centroids$cell_id) == colnames(vhdcellsce))
  
  # Extract coordinates
  coords <- st_coordinates(centroids)
  colnames(coords) <- c("pxl_col_in_fullres", "pxl_row_in_fullres")
  
  
  ## Construct SPE -------------------------------------------------------
  vhd <- SpatialExperiment(
    assays = list(counts = as(counts(vhdcellsce), "dgCMatrix")),
    rowData = rowData(vhdcellsce),
    colData = cbind(colData(vhdcellsce), coords),
    spatialCoordsNames = colnames(coords),
    scaleFactors = scalef$tissue_hires_scalef, 
    imageSources = file.path(td, "segmented_outputs", "spatial", png_name), 
    loadImage = TRUE
  )
  imgData(vhd)$image_id <- res
  
  # vhd
  
  ## Coerce to SFE -------------------------------------------------------
  vhdsfe <- toSpatialFeatureExperiment(vhd)
  # vhdsfe
  
  ## Add polygons to colGeometries
  colGeometries(vhdsfe)$cellSeg <- geo_data
  
  
  # Flip Y in colGeometries -------------------------------------------------
  # Flip Centroid Y
  centroids <- colGeometries(vhdsfe)$centroids
  colGeometries(vhdsfe)$updatecentroids <- 
    flip_sf_Y(centroids, type = "POINT", img_height = info$height, 
              scalef = imgData(vhdsfe)$scaleFactor) 
  
  # Flip Cellseg Y
  cellSeg <- colGeometries(vhdsfe)$cellSeg
  colGeometries(vhdsfe)$updatecellSeg <- 
    flip_sf_Y(cellSeg, type = "POLYGON", img_height = info$height, 
              scalef = imgData(vhdsfe)$scaleFactor) 
  
  # Flip Centroid Y in `spatialCoords()`
  spatialCoords(vhdsfe)[, "pxl_row_in_fullres"] <- 
    info$height/imgData(vhdsfe)$scaleFactor - 
    spatialCoords(vhdsfe)[, "pxl_row_in_fullres"]
  
  # Add high res coords to `colData()`
  vhdsfe[[paste0("pxl_col_in_", res)]] <- 
    spatialCoords(vhdsfe)[, "pxl_col_in_fullres"] * imgData(vhdsfe)$scaleFactor
  vhdsfe[[paste0("pxl_row_in_", res)]] <- 
    spatialCoords(vhdsfe)[, "pxl_row_in_fullres"] * imgData(vhdsfe)$scaleFactor
  
  return(vhdsfe)
}



#' Subset a `SpatialFeatureExperiment` to a region with x and y bounding box, 
#' based on the cell centroid after scaling with scalefactor.
#'
#' @param vhdsfe a `SpatialFeatureExperiment` object.
#' @param xmin lower x coordinate of the subset region. 
#' @param xmax higher x coordinate of the subset region. 
#' @param ymin lower y coordinate of the subset region. 
#' @param ymax higher y coordinate of the subset region. 
#'
#' @returns a subsetted `SpatialFeatureExperiment` object.
#'
#' @examples
#' vhdsfe_subset <- subsetVisiumHD(vhdsfe, xmin = 4214.113, xmax = 4222.916, 
#'                                         ymin = 3062.505, ymax = 3135.374)
subsetVisiumHD <- function(vhdsfe, xmin, xmax, ymin, ymax){
  # Subset to roi --------------------------------------------------------
  # From OSTA Visium HD bin-level workflow: xmin: 4214.113 ymin: 3062.505 
                                          # xmax: 4222.916 ymax: 3135.374
  
  vhdsfe_subset <- vhdsfe[, vhdsfe$pxl_col_in_hires >= xmin & 
                            vhdsfe$pxl_col_in_hires <= xmax &
                            vhdsfe$pxl_row_in_hires >= ymin & 
                            vhdsfe$pxl_row_in_hires <= ymax]
  
  return(vhdsfe_subset)
}


#' Get rank change plot before and after decontamination
#'
#' @param target_gene target gene name that exists in data.frame of stateChanges 
#' and stateChangesCon. Default "CEACAM6".
#' @param n number of top DE genes in the rank of original 
#' @param stateChanges Original rank data.frame obtained with `calcStateChanges()`
#' @param stateChangesCon After decontamination with deconvolution proportion, 
#' the updated rank data.frame obtained with `calcStateChanges()`
#'
#' @returns a ggplot object
#'
#' @examples
#' plotRankChange(target_gene = "CEACAM6", n = 20, 
#'                stateChanges,            stateChangesCon)
plotRankChange <- function(target_gene = "CEACAM6", n = 20, 
                           stateChanges, stateChangesCon) {
  rnkbdf <- stateChanges |> 
    mutate(rnk = row_number(), status = "Before") |> 
    slice_head(n = n) |> 
    select(marker, rnk, status)
  
  rnkadf <- stateChangesCon |> 
    mutate(rnk = row_number(), status = "After") |> 
    filter(marker %in% rnkbdf$marker) |>
    select(marker, rnk, status)
  
  rnk <- bind_rows(rnkbdf, rnkadf) |>
    mutate(status     = factor(status, levels = c("Before", "After")),
           gene_group = if_else(marker == target_gene, target_gene,"Other genes"))
  
  # Highlight CEACAM6
  target_df <- rnk |> filter(marker == target_gene)
  
  seg_df <- target_df |>
    tidyr::pivot_wider(names_from = status, values_from = rnk) |>
    transmute(x = "Before", y = Before, xend = "After", yend = After)
  
  # Text annotation
  before_rank <- target_df$rnk[target_df$status == "Before"][1]  # 14
  after_rank  <- target_df$rnk[target_df$status == "After"][1]   # 173
  
  before_fdr <- stateChanges$fdr[stateChanges$marker == target_gene][1]      # 0.0001078
  after_fdr  <- stateChangesCon$fdr[stateChangesCon$marker == target_gene][1] # 0.4283528
  
  label_df <- data.frame(
    status = c("Before", "After"),
    rnk    = c(before_rank, after_rank),
    label  = c(
      paste0("Rank: ", before_rank, "\n", "FDR", 
             ifelse(before_fdr < 0.05, " < 0.05", 
                    paste0(" = ", round(before_fdr, 2)))),
      paste0("Rank: ", after_rank, "\n", "FDR = ", 
             round(after_fdr, 2))
    ),
    nudge_x = c(-0.2, 0.2)  # left of Before, right of After
  )
  
  p <- ggplot(rnk, aes(status, rnk, group = marker)) +
    geom_line(alpha = 0.25, color = "grey80", show.legend = FALSE) +
    geom_point(aes(color = gene_group, size = gene_group, alpha = gene_group)) +
    geom_point(data = target_df, aes(status, rnk), color = "darkorange",size=3) +
    geom_segment(data = seg_df, aes(x, y, xend = xend, yend = yend), 
                 arrow = arrow(length = grid::unit(0.3, "cm")), 
                 size = 1.2, color = "black", inherit.aes = FALSE) +
    geom_text_repel(data = label_df, aes(x = status, y = rnk, label = label),
                    nudge_x = label_df$nudge_x, size = 3, inherit.aes = FALSE,
                    color = "black", segment.size = 0.2, min.segment.length = 0)+
    scale_color_manual(values = c("Other genes"="grey60", "CEACAM6"="darkorange"),
                       breaks = c("Other genes", "CEACAM6"), name = NULL) +
    scale_size_manual(values = c("Other genes" = 1, "CEACAM6" = 3),
                      breaks = c("Other genes", "CEACAM6"), name = NULL) +
    scale_alpha_manual(values = c("Other genes"=0.8, "CEACAM6"=1), guide="none")+
    scale_y_reverse() +
    labs(x="", y="Rank (lower = better)", title="Rank Change for Top 20 Genes") +
    theme_classic() + 
    theme(panel.grid.major = element_line(colour = "grey90", linetype = 2),
          axis.text.x = element_text(size = 13),
          plot.title = element_text(hjust = 0.5))
  
  return(p)
}

rowScale <- function(mu){
  ## 1. Get a plain numeric matrix out of your object
  mu_mat <- as.matrix(assay(mu))   # 50 x 9 numeric
  
  ## 2. Row-wise scaling: z-score per pathway
  mu_row_scaled <- t(
    apply(mu_mat, 1, function(z) {
      m <- mean(z, na.rm = TRUE)
      s <- sd(z, na.rm = TRUE)
      if (is.na(s) || s == 0) {
        rep(0, length(z))       # avoid division by zero
      } else {
        (z - m) / s
      }
    })
  )
  
  ## keep dimnames
  rownames(mu_row_scaled) <- rownames(mu_mat)
  colnames(mu_row_scaled) <- colnames(mu_mat)
  
  return(mu_row_scaled)
}

getTopSigDF <- function(mu_row_scaled, top_n = 2){
  ## mu_row_scaled is a matrix: pathways x cell types
  mu_df <- as.data.frame(mu_row_scaled) |>
    rownames_to_column("pathway") |>
    pivot_longer(
      cols = -pathway,
      names_to = "cellType",
      values_to = "mu_row_scaled"
    )
  
  # Per cell type, keep top 2 signatures
  top_mu <- mu_df |>
    group_by(cellType) |>
    slice_max(order_by = abs(mu_row_scaled), n = top_n) |>
    ungroup() |>
    mutate(pathway = fct_reorder(pathway, abs(mu_row_scaled)))
  
  return(top_mu)
}


plotDivergentBar <- function(top_mu){
  # Waterfall / diverging bar plot
  p <- ggplot(top_mu, aes(x = pathway, y = mu_row_scaled, fill = cellType)) +
    geom_hline(yintercept = 0, color = "grey70") +
    geom_col() + coord_flip() +
    scale_fill_manual(values = unname(pals::trubetskoy())) + 
    labs(x = NULL, y = "Average enrichment score (scaled)",
         title = "Top 2 pathways per cell type") +
    theme_minimal() +
    theme(legend.key.size = grid::unit(1, "lines"),
          axis.text.y = element_text(size = 11),
          panel.grid.major.y = element_blank())
  
  return(p)
  
}



