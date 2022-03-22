#' Properties of oncoEnrichR annotation datasets
#'
#' @format A data.frame with four columns, indicating size and checksums of annotation datasets
#' \itemize{
#'   \item \emph{name} - Name of oncoEnrichR annotation dataset
#'   \item \emph{size} - Size of oncoEnrichR annotation dataset (bytes)
#'   \item \emph{checksum} - Checksum of oncoEnrichR annotation dataset (R.cache::getChecksum)
#'   \item \emph{version} - Version of oncoEnrichR
#' }
#'
"db_props"

#' Versions/attributes of oncoEnrichR annotation datasets
#'
#' @format A list of lists, eech item being a list of five variables
#' \itemize{
#'   \item \emph{url} - URL for annotation/software resource
#'   \item \emph{description} - Description of resource
#'   \item \emph{version} - version
#'   \item \emph{name} - Name of resource
#'   \item \emph{resource_type} - Type of resource (db/software)
#' }
#'
"release_notes"

