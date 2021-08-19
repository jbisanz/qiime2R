#' generates a phyloseq object from .qza artifacts
#'
#' Construct a phyloseq object from multiple qiime2 artifacts (.qza). Embedded metadata for provenance is not maintained in this function and instead read_qza() should be used.
#'
#' @param features file path for artifact containing a feature (OTU/SV) table
#' @param tree file path for  artifact containing a tree
#' @param taxonomy file path for artifact containg taxonomy
#' @param metadata file path for a qiime2-compliant TSV metadata file
#' @param tmp a temporary directory that the object will be decompressed to.
#' @return a phyloseq object
#'
#' @examples \dontrun{physeq<-qza_to_phyloseq(features="data/table.qza", tree="data/rooted-tree.qza", taxonomy="data/taxonomy.qza", metdata="data/sample-metadata.qza")}
#' @export
#'
#'
#'

qza_to_phyloseq<-function(features,tree,taxonomy,metadata, tmp){

   if(missing(features) & missing(tree) & missing(taxonomy) & missing(metadata)){
    stop("At least one required artifact is needed (features/tree/taxonomy/) or the metadata.")
   }
  
if(missing(tmp)){tmp <- tempdir()}

  argstring<-""

  if(!missing(features)){
    features<-as.data.frame(read_qza(features, tmp=tmp)$data)
    argstring<-paste(argstring, "otu_table(features, taxa_are_rows=T),")
  }

  if(!missing(taxonomy)){
    taxonomy<-read_qza(taxonomy, tmp=tmp)$data
    taxonomy<-parse_taxonomy(taxonomy)
    taxonomy<-as.matrix(taxonomy)
    argstring<-paste(argstring, "tax_table(taxonomy),")
  }

  if(!missing(tree)){
    tree<-read_qza(tree, tmp=tmp)$data
    argstring<-paste(argstring, "phy_tree(tree),")
  }

  if(!missing(metadata)){
    if(is_q2metadata(metadata)){
      metadata<-read_q2metadata(metadata)
      rownames(metadata)<-metadata$SampleID
      metadata$SampleID<-NULL
    } else{
      metadata<-read.table(metadata, row.names=1, sep='\t', quote="", header=TRUE)
    }
    argstring<-paste(argstring, "sample_data(metadata),")
  }

  argstring<-gsub(",$","", argstring) #remove trailing ","

  physeq<-eval(parse(text=paste0("phyloseq(",argstring,")")))

  return(physeq)
}

