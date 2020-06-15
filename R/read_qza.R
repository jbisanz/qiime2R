#' read qiime2 artifacts (.qza)
#'
#' extracts embedded data and object metadata into an R session
#'
#' @param file path to the input file, ex: file="~/data/moving_pictures/table.qza"
#' @param tmp a temporary directory that the object will be decompressed to (default="tempdir()")
#' @param rm should the decompressed object be removed at completion of function (T/F default=TRUE)
#' @return a named list of the following objects:
#' \itemize{
#' \item artifact$data - the raw data ex OTU table as matrix or tree in phylo format
#' \item artifact$uuid - the unique identifer of the artifact
#' \item artifact$type - the semantic type of the object (ex FeatureData[Sequence])
#' \item artifact$format - the format of the qiime artifact
#' \item artifact$provenance - information tracking how the object was created
#' \item artifact$contents - a table of all the files contained within the artifact and their file size
#' \item artifact$version - the reported version for the artifact, a warning error may be thrown if a new version is seen
#' }
#'
#' @examples \dontrun{SVs<-read_qza("data/table.qza")}
#' @export
#'
#'


read_qza <- function(file, tmp, rm) {

if(missing(tmp)){tmp <- tempdir()}
if(missing(file)){stop("Path to artifact (.qza) not provided.")}
if(!file.exists(file)){stop("Input artifact (",file,") not found. Please check path and/or use list.files() to see files in current working directory.")}
if(missing(rm)){rm=TRUE} #remove the decompressed object from tmp
if(!grepl("qza$", file)){stop("Provided file is not qiime2 artifact (.qza).")}
  
  
  
unzip(file, exdir=tmp) 
unpacked<-unzip(file, exdir=tmp, list=TRUE)

artifact<-read_yaml(paste0(tmp,"/", paste0(gsub("/..+","", unpacked$Name[1]),"/metadata.yaml"))) #start by loading in the metadata not assuming it will be first file listed
artifact$contents<-data.frame(files=unpacked)
artifact$contents$size=sapply(paste0(tmp, "/", artifact$contents$files), file.size)
artifact$version=read.table(paste0(tmp,"/",artifact$uuid, "/VERSION"))


  #get data dependent on format
if(grepl("BIOMV", artifact$format)){
  artifact$data<-read_q2biom(paste0(tmp, "/", artifact$uui,"/data/feature-table.biom"))
} else if (artifact$format=="NewickDirectoryFormat"){
  artifact$data<-read.tree(paste0(tmp,"/",artifact$uuid,"/data/tree.nwk"))
} else if (artifact$format=="DistanceMatrixDirectoryFormat") {
  artifact$data<-as.dist(read.table(paste0(tmp,"/", artifact$uuid, "/data/distance-matrix.tsv"), header=TRUE, row.names=1, fill= TRUE))
} else if (grepl("StatsDirFmt", artifact$format)) {
  if(paste0(artifact$uuid, "/data/stats.csv") %in% artifact$contents$files.Name){artifact$data<-read.csv(paste0(tmp,"/", artifact$uuid, "/data/stats.csv"), header=TRUE, row.names=1)}
  if(paste0(artifact$uuid, "/data/stats.tsv") %in% artifact$contents$files.Name){artifact$data<-read.table(paste0(tmp,"/", artifact$uuid, "/data/stats.tsv"), header=TRUE, row.names=1, sep='\t')} #can be tsv or csv
} else if (artifact$format=="TSVTaxonomyDirectoryFormat"){
  artifact$data<-read.table(paste0(tmp,"/", artifact$uuid, "/data/taxonomy.tsv"), sep='\t', header=TRUE, quote="", comment="")
} else if (artifact$format=="OrdinationDirectoryFormat"){
  artifact$data<-suppressWarnings(readLines(paste0(tmp,"/", artifact$uuid, "/data/ordination.txt")))
  artifact<-parse_ordination(artifact, tmp)
} else if (artifact$format=="DNASequencesDirectoryFormat") {
  artifact$data<-readDNAStringSet(paste0(tmp,"/",artifact$uuid,"/data/dna-sequences.fasta"))
} else if (artifact$format=="AlignedDNASequencesDirectoryFormat") {
  artifact$data<-readDNAMultipleAlignment(paste0(tmp,"/",artifact$uuid,"/data/aligned-dna-sequences.fasta"))
} else if (grepl("EMPPairedEndDirFmt|EMPSingleEndDirFmt|FastqGzFormat|MultiplexedPairedEndBarcodeInSequenceDirFmt|MultiplexedSingleEndBarcodeInSequenceDirFmt|PairedDNASequencesDirectoryFormat|SingleLanePerSamplePairedEndFastqDirFmt|SingleLanePerSampleSingleEndFastqDirFmt", artifact$format)) {
  artifact$data<-data.frame(files=list.files(paste0(tmp,"/", artifact$uuid,"/data")))
  artifact$data$size<-format(sapply(artifact$data$files, function(x){file.size(paste0(tmp,"/",artifact$uuid,"/data/",x))}, simplify = TRUE))
} else if (artifact$format=="AlphaDiversityDirectoryFormat") {
  artifact$data<-read.table(paste0(tmp, "/", artifact$uuid, "/data/alpha-diversity.tsv"))
} else if (artifact$format=="DifferentialDirectoryFormat") {
  defline<-suppressWarnings(readLines(paste0(tmp, "/", artifact$uuid, "/data/differentials.tsv"))[2])
  defline<-strsplit(defline, split="\t")[[1]]
  defline[grep("numeric", defline)]<-"double"
  defline[grep("categorical|q2:types", defline)]<-"factor"
  coltitles<-strsplit(suppressWarnings(readLines(paste0(tmp, "/", artifact$uuid, "/data/differentials.tsv"))[1]), split='\t')[[1]]
  artifact$data<-read.table(paste0(tmp, "/", artifact$uuid, "/data/differentials.tsv"), header=F, col.names=coltitles, skip=2, sep='\t', colClasses = defline, check.names = FALSE)
  colnames(artifact$data)[1]<-"Feature.ID"
} else {
  message("Format not supported, only a list of internal files and provenance is being imported.")
  artifact$data<-list.files(paste0(tmp,"/",artifact$uuid, "/data"))
}

#Add Provenance
pfiles<-paste0(tmp,"/", grep("..+provenance/..+action.yaml", unpacked$Name, value=TRUE))
artifact$provenance<-lapply(pfiles, read_yaml)
names(artifact$provenance)<-grep("..+provenance/..+action.yaml", unpacked$Name, value=TRUE)

if(rm==TRUE){unlink(paste0(tmp,"/", artifact$uuid), recursive=TRUE)}


return(artifact)
}
