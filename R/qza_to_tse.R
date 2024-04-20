#' Generates a TreeSummarizedExperiment (TSE) object from .qza artifacts
#'\code{\link[TreeSummarizedExperiment:SummarizedExperiment-class]{TreeSummarizedExperiment}}
#'
#' Construct a TreeSummarizedExperiment object from multiple 
#' qiime2 artifacts (.qza). 
#' '\code{\link[TreeSummarizedExperiment:SummarizedExperiment-class]{TreeSummarizedExperiment}}
#' Embedded metadata for provenance is not maintained in this function and 
#' instead \code{read_qza()} should be used.
#' 
#' @param features file path for artifact containing a feature (OTU/SV) table
#' @param tree file path for artifact containing a tree
#' @param taxonomy file path for artifact containg taxonomy
#' @param metadata file path for a qiime2-compliant TSV metadata file
#' @param tmp a temporary directory that the object will be decompressed to.
#' @return a TreeSummarizedExperiment object
#' \code{\link[TreeSummarizedExperiment:SummarizedExperiment-class]{TreeSummarizedExperiment}}
#' 
#' @import TreeSummarizedExperiment
#' @importFrom S4Vectors SimpleList DataFrame make_zero_col_DFrame
#' @importFrom SummarizedExperiment colData colData<-
#' 
#' @examples 
#' (Data is from tutorial
#' \url{https://docs.qiime2.org/2020.2/tutorials/moving-pictures/)}
#' 
#' \donttest{tse <-qza_to_tse(
#'     features="path_to_table.qza",
#'     tree="path_to_rooted-tree.qza",
#'     taxonomy="path_to_taxonomy.qza",
#'     metadata = "path_to_sample-metadata.tsv"
#' )}
#' 
#' @export
#' 
#' @author
#' Leo Lahti
#' Noah de Gunst

qza_to_tse<- function(features,tree,taxonomy,metadata, tmp) {
    # Input check
    if(missing(features)){
        stop("No features have been specified. Please specify one to create a TSE from qza data.")
    }
    
    # If no temporary extract location is specified, use default
    if(missing(tmp)){tmp <- tempdir()}
    
    if(!missing(features)){
        # Read the qza features
        feat <- read_qza(features, tmp=tmp)
        # Transform the data to a matrix
        features <- as.matrix(feat$data)
        # Create a list of assays
        assays <- S4Vectors::SimpleList(counts=features)
        # Get the metadata
        meta <- feat
        meta$data <- NULL # remove data as it is put in features
        
    }
    
    if(!missing(taxonomy)){
        # Read the taxonomy data
        rowData <- read_qza(taxonomy, tmp=tmp)$data
        rowData <- parse_taxonomy(rowData)
        # Convert to S4 DataFrame
        rowData <- S4Vectors::DataFrame(data.frame(rowData))
    }
    else{
        # Creates empty table with zero columns if no taxonomy table is specified 
        rowData <- S4Vectors::make_zero_col_DFrame(nrow(assays$counts))
    }
    # Store the rownames
    rownames(rowData) <- rownames(assays$counts)
    
    if(!missing(metadata)){
        # Reads metadata
        if(is_q2metadata(metadata)){
            colData <- read_q2metadata(metadata)
        } else {
            colData <- read.table(metadata, row.names=1, sep='\t', quote="", header=TRUE)
        }
        # Converts to S4 DataFrame
        colData <- S4Vectors::DataFrame(data.frame(colData)) 
        
    }
    else{
        # If there's no metadata available then set the Coldata to zeros
        colData <- S4Vectors::make_zero_col_DFrame(ncol(assays$counts))
    }
    # Stores the colnames of the counts assay as colData
    rownames(colData) <- colnames(assays$counts) 
    
    if(!missing(tree)){
        # Reads the taxonomy tree if one is specified
        rowTree<-read_qza(tree, tmp=tmp)$data
    }
    else {
        rowTree <- NULL
    }

    # Creates a TSE from the given data
    TreeSummarizedExperiment::TreeSummarizedExperiment(assays = assays,
                             rowData = rowData,
                             colData = colData,
                             rowTree = rowTree,
                             referenceSeq = NULL,
                             metadata = meta)
}
