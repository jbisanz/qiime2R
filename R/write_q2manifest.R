#' Generates a  read manifest for importing sequencing data into qiime2
#'
#' Scans a directory for files with matching sequencing data (default: fastq.gz) and then generates a q2 compliant manifest.
#'
#' @param outfile filename for output (default: manifest_[timestamp].txt)
#' @param directory directory containing reads
#' @param extension file extension (default: .fastq.gz)
#' @param paired are reads in paired format? TRUE/FALSE (default=FALSE)
#' @param Fwd string used to denote a forward read (default= _R1)
#' @param Rev string used to denote a reverse read (default= _R2)
#'
#' @examples \dontrun{write_q2manifest("q2manifest.txt","/yourdirhere/reads/", extension=".fastq.gz", paired=TRUE)}
#' @export
#'

write_q2manifest<-function(outfile, directory, extension, paired, Fwd, Rev){
  if(missing(outfile)){outfile<-paste0("q2manifest_", gsub(" |:","", timestamp(prefix="", suffix="")),".txt")}
  if(missing(directory)){stop("Directory containing reads not provided.")}
  if(missing(extension)){extension=".fastq.gz"}
  if(missing(paired)){paired=FALSE}
  if(missing(Fwd)){Fwd="_R1"}
  if(missing(Rev)){Rev="_R2"}
  
  files<-list.files(directory, pattern=gsub("\\.", "\\.", extension))
  
  if(!paired){
    output<-data.frame(`sample-id`=gsub(gsub("\\.", "\\.", extension), "", files), `absolute-filepath`=paste0(getwd(), "/", files), check.names = FALSE)
    write.table(output, file=outfile, row.names=F, quote=F, sep="\t")
  } else {
    output<-data.frame(`sample-id`=gsub(gsub("\\.", "\\.", extension), "", files), file=paste0(getwd(), "/", files), check.names = FALSE)
    output$Read<-case_when(
      grepl(Fwd, output$file)~"forward-absolute-filepath",
      grepl(Rev, output$file)~"reverse-absolute-filepath",
      TRUE~"Error"
    )
    if(!sum(grep("Error", output$Read))==0){stop("Could not assign all reads to forward or reverse in paired mode.")}
    output$`sample-id`<-gsub(Fwd, "", output$`sample-id`)
    output$`sample-id`<-gsub(Rev, "", output$`sample-id`)
    output<-spread(output, key=Read, value=file)
    write.table(output, file=outfile, row.names=F, quote=F, sep="\t")
  }
  
  message("Manifest written to", outfile)
}

