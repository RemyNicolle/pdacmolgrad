
.getUGM=function (m, g, w)
{
  if (!all.equal(nrow(m), length(g), length(w))) {
    stop("nrow of m should be equal to lenght of g and w")
  }
  i = order(w, decreasing = T)
  oki = i[which(!duplicated(g[i]) & !g[i] %in% c("---", " ",
                                                 "", NA))]
  okm = m[oki, ]
  rownames(okm) = g[oki]
  okm
}




#' projectMolGrad
#'
#' @param newexp gene expression matrix or dataframe with gene in row and sample in columns.
#' @param geneSymbols vector of gene symbols for the newexp dataset (simply set to rownames if the nexzexp is already in single values per gene symbols)
#'
#' @return data frame of three projections
#' @export
#'
#' @examples
projectMolGrad=function(newexp,geneSymbols){

  if(nrow(newexp)!= length(geneSymbols)){
    stop("geneSymbols should be a vector of gene symbols exactly corresponding to each row of the newexp dataset")
  }
  expg=pdacmolgrad:::.getUGM(newexp,geneSymbols,rowSds(as.matrix(newexp)))



  comg=intersect(rownames(pdacmolgrad:::molgradeSystems$pdx),rownames(expg))
  if(length(comg)<5000){
    pdxproj=NA
  }else{
    invs=MASS::ginv(as.matrix(molgradeSystems$pdx[comg,]))
    scexp=scale(expg[comg,])
    pdxprojInvSc=t(t((t(	scexp )%*% t(invs))))[,1]
    pdxproj=(pdxprojInvSc-mean(pdxprojInvSc))/(3*sd(pdxprojInvSc))
  }

  comg=intersect(rownames(pdacmolgrad:::molgradeSystems$icgc),rownames(expg))
  scexp=scale(expg[comg,])
  if(length(comg)<5000){
    icgcproj=NA
  }else{
    invs=MASS::ginv(as.matrix(pdacmolgrad:::molgradeSystems$icgc[comg,]))
    icgcprojInvSc=t(t((t(	scexp )%*% t(invs))))[,1]
    icgcproj=(icgcprojInvSc-mean(icgcprojInvSc))/(3*sd(icgcprojInvSc))

  }
  comg=intersect(rownames(pdacmolgrad:::molgradeSystems$puleo),rownames(expg))
  scexp=scale(expg[comg,])
  if(length(comg)<5000){
    puleoproj=NA
  }else{
    invs=MASS::ginv(as.matrix(pdacmolgrad:::molgradeSystems$puleo[comg,]))
    puleoprojInvSc=t(t((t(	scexp )%*% t(invs))))[,2]
    puleoproj=(puleoprojInvSc-mean(puleoprojInvSc))/(3*sd(puleoprojInvSc))
  }
  data.frame(pdxproj=pdxproj,
             icgcproj=icgcproj,
             puleoproj=puleoproj,
             row.names=colnames(expg),stringsAsFactors = F)
}
