#' Parse an imported q2 ordination text file
#'
#' @param artifact takes an artifact containing the unparsed ordination format (result of readLines) 
#' @param tmp the tmp dir for artifact 
#' 
#' @return a data.frame with feature IDs as row names and the columns: Kingdom, Phylum, Class, Order, Family, Genus, Species
#'
#' @examples \dontrun{ord<-parse_ordination(ordination_artifact, tmpdir)}
#' @export
#'


parse_ordination <- function(artifact, tmp){
  if(missing(artifact)){stop("Ordination not provided.")}
  if(missing(tmp)){stop("Temp directory not passed.")}
  
  artifact$data<-artifact$data[sapply(artifact$data, function(x) x!="")]
  
  for (i in 1:length(artifact$data)){
    if(grepl("^Eigvals\\t|^Proportion explained\\t|^Species\\t|^Site\\t|^Biplot\\t|^Site constraints\\t", artifact$data[i])){
      curfile=strsplit(artifact$data[i],"\t")[[1]][1]
    } else {
      write(artifact$data[i], paste0(tmp,"/", artifact$uuid, "/data/",curfile,".tmp"), append=TRUE)
    }
  }
  
  backup<-artifact$data
  artifact$data<-list()
  for (outs in list.files(paste0(tmp,"/", artifact$uuid,"/data"), full.names = TRUE, pattern = "\\.tmp")){
    NewLab<-gsub(" ", "", toTitleCase(gsub("\\.tmp", "", basename(outs))))
    artifact$data[[NewLab]]<-read.table(outs,sep='\t', header=FALSE)
    if(NewLab %in% c("Eigvals","ProportionExplained")){colnames(artifact$data[[NewLab]])<-paste0("PC",1:ncol(artifact$data[[NewLab]]))}
    if(NewLab %in% c("Site","SiteConstraints")){colnames(artifact$data[[NewLab]])<-c("SampleID", paste0("PC",1:(ncol(artifact$data[[NewLab]])-1)))}
    if(NewLab %in% c("Species")){colnames(artifact$data[[NewLab]])<-c("FeatureID", paste0("PC",1:(ncol(artifact$data[[NewLab]])-1)))}
  }
  artifact$data$Vectors<-artifact$data$Site #Rename Site to Vectors so this matches up with the syntax used in the tutorials
  artifact$data$Site<-NULL
  return(artifact)
}
