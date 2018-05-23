#' generates a phyloseq (.qza)
#'
#' Construct a phyloseq object from multiple qiime2 artifacts (.qza). Embedded metadata for provenance is not maintained in this function and instead read_qza() should be used.
#'
#' @param features file path for artifact containing a feature (OTU/SV) table
#' @param tree file path for  artifact containing a tree
#' @param taxonomy file path for artifact containg taxonomy
#' @param metadata file path for a qiime2-compliant TSV metadata file
#' @return a phyloseq object
#' @usage physeq<-qza_to_phyloseq(features="table.qza", tree="rooted-tree.qza", taxonomy="taxonomy.qza", metadata="sample-metadata.tsv")
#' @export
#'
#'
#'

qza_to_phyloseq<-function(features,tree,taxonomy,metadata){

   if(missing(features) & missing(tree) & missing(taxonomy) & missing(metadata)){
    stop("At least one required artifact is needed (features/tree/taxonomy/) or the metadata.")
  }


  argstring<-""

  if(!missing(features)){
    features<-read_qza(features)$data
    argstring<-paste(argstring, "otu_table(features, taxa_are_rows=T),")
  }

  if(!missing(taxonomy)){
    taxonomy<-read_qza(taxonomy)$data
    taxt<-suppressWarnings(do.call(rbind, strsplit(as.character(taxonomy$Taxon),"\\; ")))
    rownames(taxt)<-taxonomy$Feature.ID
    argstring<-paste(argstring, "tax_table(taxt),")
  }

  if(!missing(tree)){
    tree<-read_qza(tree)$data
    argstring<-paste(argstring, "phy_tree(tree),")
  }

  if(!missing(metadata)){
    metadata<-read.table(metadata, row.names=1, sep='\t', quote="")
    argstring<-paste(argstring, "sample_data(metadata),")
    sample_data(metadata)
  }

  argstring<-gsub(",$","", argstring) #remove trailing ","

  physeq<-eval(parse(text=paste0("phyloseq(",argstring,")")))

  return(physeq)
}
