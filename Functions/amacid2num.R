amacid2num <- function(amname){
  v <- switch (amname,
       "PHE" = 1,
       "LEU" = 2,
       "ILE" = 3,
       "VAL" = 4,
       "TRP" = 5,
       "TYR" = 6,
       "MET" = 7,
       "CYH" = 8,
       "CYS" = 8,
       "CYSS" =8,
       "HIS" = 9,
       "THR" = 10,
       "ARG" = 11,
       "ALA" = 12,
       "ASN" = 13,
       "GLN" = 14,
       "PRO" = 15,
       "SER" = 16,
       "ASP" = 17,
       "GLY" = 18,
       "LYS" = 19,
       "LYZ" = 19,
       "GLU" = 20,
       21
       )
  v
  }

