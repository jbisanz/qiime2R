#' Functions for decoding golay barcodes
#'
#' This is a port of the original QIIME2 python code from q2_demux/_ecc.py into R made by Cecilia Noecker 2020/11/17 and then added into qiime2R by JB
#' @param barcodes a vector of barcodes to be corrected
#' @return a data.table with the following columns: original_bc, corrected_bc, n_errors
#'
#' @export
#'

decode_golay<-function(barcodes){
  if(missing(barcodes) ){
    stop("At least one barcode is needed for correction")
  }
  nt_bits_mapping = data.table(Base = c("A", "C", "T", "G"), Bits = c("11", "00", "10", "01"))
  
  # Parity submatrix
  default_p = matrix(   rbind(c(0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
                              c(1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0),
                              c(1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1),
                              c(1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1),
                              c(1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0),
                              c(1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1),
                              c(1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1),
                              c(1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1),
                              c(1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0),
                              c(1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0),
                              c(1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0),
                              c(1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1 )), nrow = 12)
  
  default_g = cbind(default_p, diag(12))
  
  #Generator matrix
  # matrix H satisfies G dot H.T = zeros (mod 2 arithmetic)
  # also satisfies syn = H dot rec = H dot err (rec is recieved 24 bits,
  # err is 24 bit error string added to transmitted 24 bit vec)
  # (all mod 2 arithmetic)
  default_h = cbind(diag(12), t(default_p))
  
  make_3bit_errors = function(veclen = 24){
    nvecs = veclen + choose(veclen, 2) + choose(veclen, 3)
    
    # +1 for an all zero vector
    errorvecs = matrix(rep(0), nrow = nvecs + 1, ncol = veclen)  
    # one 1
    offset = 2
    for(i in 1:veclen){
      vec = errorvecs[offset,]
      vec[i] = 1
      errorvecs[offset,] = vec
      offset = offset + 1
    }
    
    # two 1s
    for(i in 1:(veclen-1)){
      for(j in (i + 1):veclen){
        vec = errorvecs[offset,]
        vec[i] = 1
        vec[j] = 1
        errorvecs[offset,] = vec
        offset = offset + 1
      }
    }
    
    # three 1s
    for(i in 1:(veclen-2)){
      for(j in (i+1):(veclen-1)){
        for(k in (j+1):veclen){
          vec = errorvecs[offset,]
          vec[i] = 1
          vec[j] = 1
          vec[k] = 1
          errorvecs[offset,] = vec
          offset = offset + 1
        }
      }
    }
    return(errorvecs)
  }
  
  make_syndrome_lookup_table = function(){
    syndrome_lookup_table = data.table(SynID = rep("", nrow(all_3bit_errors)), Errvec = rep("", nrow(all_3bit_errors)), ErrvecNumeric = vector("list", length = nrow(all_3bit_errors)))
    # build syndrome lookup table
    for(j in 1:nrow(all_3bit_errors)){
      errvec = all_3bit_errors[j,]
      syn = (default_h %*% errvec) %% 2
      syn_id = paste(as.character(syn), collapse = "")
      syndrome_lookup_table[j, SynID:=syn_id] 
      #  syndrome_lookup_table[j, Syn:=t(syn)]
      syndrome_lookup_table[j, Errvec:=paste(as.character(errvec), collapse = "")]
      syndrome_lookup_table[j, ErrvecNumeric:=list(errvec)]
    }
    return(syndrome_lookup_table)
  }
  
  rev_comp = function(seq){
    if(length(seq)==1){
      return(as.character(reverseComplement(DNAString(seq))))
    } else {
      complement_tab = data.table(Base= c("A", "C", "T", "G"), Base2 = c("T", "G", "A", "C"))
      seq_nts = data.table(Base = unlist(sapply(seq, strsplit, split = "")), 
                           Seq = c(sapply(1:length(seq), rep , times=12)))
      seq_nts = merge(seq_nts,complement_tab, by = "Base", all.x = T, sort = F)
      return(seq_nts[,paste(rev(Base2), collapse = ""), by = Seq][,V1])
    }
  }
  
  seq_to_bits = function(seq){
    if(length(seq) == 1){
      seq_nts = data.table(Base = strsplit(seq, split = "")[[1]])
      seq_nts = merge(seq_nts, nt_bits_mapping, by = "Base", all.x = T, sort = F)
      return(paste(seq_nts[,Bits], collapse = ""))
    } else {
      seq_nts = data.table(Base = unlist(sapply(seq, strsplit, split = "")), Seq = c(sapply(1:length(seq), rep , times=12)))
      seq_nts = merge(seq_nts, nt_bits_mapping, by = "Base", all.x = T, sort = F)
      return(seq_nts[,paste(Bits, collapse = ""), by = Seq][,V1])
    }
  }
  
  bits_to_seq = function(bits){
    if(!is.list(bits)){
      nts = data.table(Bits = paste0(bits[c(TRUE, FALSE)], bits[c(FALSE, TRUE)]))
      nts = merge(nts, nt_bits_mapping, by =  "Bits", all.x = T, sort = F)
      seq = paste0(nts[,Base], collapse = "")
      
    } else {
      nts = data.table(Bits = sapply(bits, function(x) {paste0(x[c(TRUE, FALSE)], x[c(FALSE, TRUE)]) }), Seq = 1:length(bits))
      nts = merge(nts, nt_bits_mapping, by =  "Bits", all.x = T, sort = F)
      seq = nts[,paste0(nts[,Base], collapse = ""), by = Seq][,V1]
    }
    return(seq)
  }
  
  all_3bit_errors = make_3bit_errors()
  syndrome_lookup_table = make_syndrome_lookup_table()

  seq_table<-data.table(StartingBarcode=barcodes)
  #taking the core of CN's vectorized function
    seq_table[,GoodBarcode:=grepl("^[ACTG]*$", StartingBarcode)] #ignore barcodes with Ns
    if(nrow(seq_table[GoodBarcode == T]) == 0) { 
	stop("No valid barcodes found. Do your fastq files have index sequences in the sequence headers? They are required for Golay decoding of undetermined reads. Consider checking your bcl2fastq settings and re-running FASTQ generation.")
    }
    seq_table[GoodBarcode==T ,Bits:=seq_to_bits(StartingBarcode)]
    seq_table[GoodBarcode==T, BitsNumeric:=lapply(Bits, function(x){
      as.numeric(strsplit(x, split = "")[[1]])
    })]
    seq_table[GoodBarcode==T, SynID:=sapply(BitsNumeric, function(x) {
      paste(((default_h %*% x) %% 2), collapse = "") })]
    seq_table = merge(seq_table, syndrome_lookup_table, by = "SynID", all.x = T, sort = F)
    seq_table[!is.na(Errvec), Corrected:=lapply((Map(`+`, BitsNumeric, ErrvecNumeric)), function(x){x %%2 })]
    seq_table[!is.na(Errvec), CorrectedSeq:=sapply(Corrected, bits_to_seq)]
    seq_table[!is.na(Errvec), NumErrors:=sapply(ErrvecNumeric, function(x){ sum(x==1)})]

    seq_table<-seq_table %>% select(original_bc=StartingBarcode, corrected_bc=CorrectedSeq, n_errors=NumErrors)
  
    return(seq_table)
}














