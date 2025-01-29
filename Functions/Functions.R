# -----------------------------------------------------------------------------
rm(list = ls())
# -----------------------------------
# -------------------------------------------------------------
# ------------------------------
# ---- Load Knowledge_based potential from structure and sequence
# --------------------------------------------
# --------------------------------------------------------
# -------------------------------------------------------------------------------------
energy_dell_dunbrack <- read.csv("Data/csv/energy.csv",header = FALSE,sep ="," )
energy_dell_dunbrack <- energy_dell_dunbrack[,5:33]
# ------Pairwise Amino Acid names-------------
AA <- c("PHE","LEU","ILE","VAL","TRP","TYR","MET","CYS","HIS","THR",
        "ARG","ALA","ASN","GLN","PRO","SER","ASP","GLY","LYS","GLU")
letters_list <- bio3d::aa321(AA)
pair_inter <- character(0)
for (i in 1:20) {
  for (j in i:20) {
    pair_inter <- c(pair_inter, paste0(letters_list[i], letters_list[j]))
  }
} 
# -----------------------------------
# -----  Sequence Energy Estimator (predictor)
aaenergy <- read.csv("Data/csv/Pij.csv", header = F,sep = ";")
aaenergy <- (aaenergy + t(aaenergy))/2
colnames(aaenergy) <- rownames(aaenergy) <- letters_list
# -------------------------------------------------------------
rm(i,j,AA)
# -----------
# -------------------------
# ------ This function code Amino Acid name into number.
# ------------------------------------------------------------
amacid2num <- readRDS('Functions/amacid2num.rds')
# -------------------------------------------------
# -----------------------
# ----------- This function code pair of Atom types into numbers.
# -----------------------------------------------------
atomtype2num167 <- readRDS('Functions/atomtype2num167.rds')
# ------------------------------------------------------
# ---------------------------
# ------------- This function extract Chain ID and Start and End positions for each domain.
# -------------------------------------------------
split_position <- function(position){
  # position means the chain and residue numbers
  input_string <- position
  split_strings <- unlist(strsplit(input_string, ",")) # split the string by comma
  chains <- starts <- ends <-character(0)
  # Process each part of the string
  for (part in split_strings) {
    parts <- strsplit(part, ":")[[1]]
    chain <- parts[1]
    if (length(parts) > 1) {
      range <- unlist(strsplit(parts[2], "-"))
      start <- ifelse(length(range) == 1, NA, as.numeric(range[1]))
      end <- ifelse(length(range) == 1, NA, as.numeric(range[2]))
    } else {
      start <- NA
      end <- NA
    }
    chains <- c(chains, chain)
    starts <- c(starts, start)
    ends <- c(ends, end)
  }
  df <- data.frame(chain = chains, start = starts, end = ends)
  df[is.na(df)] <- "NA"
  return(df)
}
# ----------------------------------------------------
# -------------------------------
# ------------ This function convert PDB structure into edge list of contact graph. 
# -----------------------------------------------
NetworkFrame <- readRDS('Functions/NetworkFrame.rds')
# ------------------------------------------------
# --------------------------------------------------------
# --------------------------------
# ------------- This function convert the contact graph into the 210D profile energy.
# ----------------------------------------------------------
Energy_SPE210 <- readRDS('Functions/Energy_SPE210.rds')
# --------------------------------------------------
# --------------------------------
# ------------- This function calculate the 210D profile energy from Amino Acid sequence.
# ---------------------------------------------------------
Energy_CPE210 <- readRDS('Functions/Energy_CPE210.rds')
# ---------------------------------------------------------------
# ------------------------------------
# ------------ This function calculates the distance between two profiles of energy.
# -----------------------------------------------------------
manhat <- function(X, Y){
  x<-X; y<-Y
  xy <- abs(x - y)
  xy <- xy[xy != 0]
  if(length(xy)==0){
    d <- 0
  }else{
    d <- sum(xy)/length(xy)
  }
  return(d)
}
# ---------------------------------------------------
# -----------------------------------
# ---------- Customized theme for figures.
# ---------------------------------------------------------
mytheme <- function(){
  # ------------  theme ----------------
  mythem <- theme_bw(base_line_size = .2)+
    theme(axis.text.x = element_text(size=13, face = "bold"),
          axis.text.y = element_text(size=13, face = "bold"),
          axis.title.x = element_text(size = 20, face = "bold"),
          axis.title.y = element_text(size = 20, face = "bold"),
          legend.text = element_text(size = 15, colour = "black"),
          legend.title = element_blank(),
          title = element_text(size = 15)
    )
  return(mythem)
}
# --------------------------------
# Create or append to the file
MYwriteData <- function(file, sheet, data, row.Names = FALSE) {
  if (file.exists(file)) {
    wb <- loadWorkbook(file)
  } else {
    wb <- createWorkbook()
  }
  addWorksheet(wb, sheet)
  writeData(wb, sheet, data,rowNames = row.Names)
  saveWorkbook(wb, file, overwrite = TRUE)
}
