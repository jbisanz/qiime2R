#' read qiime2 artifacts (.qza)
#'
#' extracts embedded data and object metadata into an R session
#'
#' @param file path to the input file, ex: file="~/data/moving_pictures/table.qza"
#' @param tmp a temporary directory that the object will be decompressed to (default="/tmp")
#' @param rm should the decompressed object be removed at completion of function (T/F default=T)
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
#' @export
#'
#'


read_qza <- function(file, tmp, rm) {

if(missing(tmp)){tmp="/tmp/"}
if(missing(file)){stop("Path to artifact (.qza) not provided")}
if(missing(rm)){rm=T} #remove the decompressed object from tmp

untar(file, exdir=tmp)
unpacked<-untar(file, exdir=tmp, list=T)

artifact<-read_yaml(paste0(tmp,"/", grep(paste0(gsub("/..+","", unpacked[1]),"/metadata.yaml"), unpacked, value=T) )) #start by loading in the metadata not assuming it will be first file listed
artifact$contents<-data.frame(files=unpacked)
artifact$contents$size=sapply(paste0(tmp, "/", artifact$contents$files), file.size)
artifact$version=read.table(paste0(tmp,"/",artifact$uuid, "/VERSION"))

if(sum(artifact$version$V2==c("2","4","2018.4.0"))!=3){warning("It is unclear if current artifact version is supported, errors may occur in import...")}#check version and throw warning if new format

  #get data dependent on type
if(grepl("FeatureTable\\[[A-z]+\\]", artifact$type)){ #FeatureTable[Balance]|FeatureTable[Composition]|FeatureTable[Frequency]|FeatureTable[PresenceAbsence]|FeatureTable[RelativeFrequency]
  suppressWarnings(artifact$data<-as(biom_data(read_biom(paste0(tmp, "/", artifact$uui,"/data/feature-table.biom"))),"matrix")) #suppressing warning about \n characters
} else if (grepl("Phylogeny\\[[A-z]+\\]", artifact$type)){ #Phylogeny[Rooted]|Phylogeny[Unrooted]
  artifact$data<-read.tree(paste0(tmp,"/",artifact$uuid,"/data/tree.nwk"))
} else if (artifact$type=="DistanceMatrix") { #DistanceMatrix
  artifact$data<-as.dist(read.table(paste0(tmp,"/", artifact$uuid, "/data/distance-matrix.tsv"), header=T, row.names=1))
} else if (artifact$type %in% c("DeblurStats", "QualityFilterStats")) { #DeblurStats|QualityFilterStats
  artifact$data<-read.csv(paste0(tmp,"/", artifact$uuid, "/data/stats.csv"), header=T, row.names=1)
} else if (artifact$type=="FeatureData[Taxonomy]"){ #FeatureData[Taxonomy]
  artifact$data<-read.table(paste0(tmp,"/", artifact$uuid, "/data/taxonomy.tsv"), sep='\t', header=T)
} else if (artifact$type=="PCoAResults"){

  results<-scan(file=paste0(tmp,"/", artifact$uuid, "/data/ordination.txt"), what="character", sep='\t') #deal with merged table by reading as one long string

  start<-grep("Proportion explained",results)+2
  len<-as.numeric(results[grep("Proportion explained",results)+1])
  artifact$data$ProportionExplained<-as.numeric(results[start:(start+(len-1))])
  names(artifact$data$ProportionExplained)<-paste0("PC",1:len)

  results<-results[(grep("Site", results)[1]+3):(grep("Biplot", results)[1]-1)]
  artifact$data$Vectors<-as.data.frame(t(matrix(results, ncol=len)))
  colnames(artifact$data$Vectors)<-c("SampleID", paste0("PC", 1:len))
  artifact$data$Vectors[,2:len]<-apply(artifact$data$Vectors[,2:len], 2, as.numeric)
} else if (artifact$type=="FeatureData[Sequence]") {
  artifact$data<-readDNAStringSet(paste0(tmp,"/",artifact$uuid,"/data/dna-sequences.fasta"))
} else if (artifact$type=="FeatureData[AlignedSequence]") {
  artifact$data<-readDNAMultipleAlignment(paste0(tmp,"/",artifact$uuid,"/data/aligned-dna-sequences.fasta"))
} else if (grepl("Sequences", artifact$type)) {
  artifact$data<-data.frame(files=list.files(paste0(tmp,"/", artifact$uuid,"/data")))
  artifact$data$size<-format(sapply(artifact$data$files, function(x){file.size(paste0(tmp,"/",artifact$uuid,"/data/",x))}, simplify = T))
} else {
  message("Semantic type not supported, only a list of data files and provenance is being imported.")
  artifact$data<-list.files(paste0(tmp,"/",artifact$uuid, "/data"))
}

pfiles<-paste0(tmp,"/", grep("..+provenance/..+action.yaml", unpacked, value=T))
artifact$provenance<-lapply(pfiles, read_yaml)
names(artifact$provenance)<-grep("..+provenance/..+action.yaml", unpacked, value=T)
if(rm==T){unlink(paste0(tmp,"/", artifact$uuid), recursive=T)}
return(artifact)
}
