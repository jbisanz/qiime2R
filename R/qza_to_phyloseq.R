#' generates a phyloseq (.qza)
#'
#' Construct a phyloseq object from multiple qiime2 artifacts (.qza). Embedded metadata for provenance is not maintained in this function and instead read_qza() should be used.
#'
#' @param features file path for artifact containing a feature (OTU/SV) table
#' @param tree file path for  artifact containing a tree
#' @param taxonomy file path for artifact containg taxonomy
#' @param metadata file path fora qiime2-compliant TSV metadata file
#' @return a phyloseq object
#' @usage physeq<-qza_to_phyloseq(features="table.qza", tree="rooted-tree.qza", taxonomy="taxonomy.qza", metadata="sample-metadata.tsv")
#' @export
#'
#'
#'

#to do: be more flexible in which objects are put into the phyloseq object

qza_to_phyloseq<-function(features,tree,taxonomy,metadata){

   if(missing(features) | missing(tree) | missing(taxonomy) | missing (metadata)){
    stop("At least one required artifact was missing (features/tree/taxonomy/) or the metadata.")
  }

  tax<-read_qza(taxonomy)$data
  taxt<-do.call(rbind, strsplit(as.character(tax$Taxon),"\\; "))
  rownames(taxt)<-tax$Feature.ID

  physeq<-
    phyloseq(
      otu_table(read_qza(features)$data, taxa_are_rows=T),
      phy_tree(read_qza(tree)$data),
      tax_table(taxt),
      sample_data(read.table(metadata, row.names=1, sep='\t', quote=""))
    )

  return(physeq)
}
