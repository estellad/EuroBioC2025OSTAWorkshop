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





