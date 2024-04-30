amacid2numSC <- function(am2numi,numi){
  idxSC<-c(9, 14:16, 21:24, 29:32, 37:40, 45:47, 52:57,
           62:68, 73:80, 85:94, 99:100, 105:106, 111:113,
           118:121, 126:130, 135:138, 143:147, 152:156, 161:167)
  aaSC <- c(21, 22, 23, 24, 25,
                           26, 27, 28,
                           29, 30, 31, 32, 33,
                           34, 35, 36, 37, NA,38,
                           39)
  if(!any(numi==idxSC)){
    v <- am2numi
  }else if(any(numi==idxSC)){
    v <- aaSC[am2numi]
  }
  return(v)
}

