inedx_HelixSheet <- function(PDB){
  hs <- PDB$helix$start; he <- PDB$helix$end
  ss <- PDB$sheet$start; se <- PDB$sheet$end
  indhs<-c()
  if(!is.null(hs)|length(hs)!=0){
    for (i in 1:length(hs)) {
      h <- hs[i]:he[i]
      indhs <- c(indhs,h)
    }
  }
  if(!is.null(ss)|length(ss)!=0){
    for (i in 1:length(ss)) {
      s <- ss[i]:se[i]
      indhs <- c(indhs,s)
    }
  }
  return(indhs)
}
