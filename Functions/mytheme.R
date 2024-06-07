mytheme <- function(){
#  # ------------- Execution function -------------
#  run_chunk <- function(chunk_label) {
#    file_name <- 'TargetedOlink_Mass_Proteomics.Rmd'
#    x <- readLines(paste0("./", file_name))
#    idxStart <- grepl(paste0("\\{r ", chunk_label), x)
#    if(length(which(idxStart))==0){
#      stop(paste0('there is no chunk called ‘',chunk_label,'’'))
#    }
#    chunk_start <- (which(idxStart) + 1)
#    chunk_signs <- grepl("```", x)
#    idxEnd <- which(which(chunk_signs) > chunk_start)[1]
#    chunk_end <- (which(chunk_signs)[idxEnd] - 1)
#    chunk_lines <- c(chunk_start:chunk_end)
#    eval(parse(text=strwrap(paste0(x[chunk_lines]))), envir = .GlobalEnv)
#  }
  # ------------------------------------
  # ------------  theme ----------------
  mythem <- theme_bw(base_line_size = .2)+
    theme(axis.text.x = element_text(size=13, face = "bold"),
          axis.text.y = element_text(size=13, face = "bold"),
          axis.title = element_text(size = 20, face = "bold"),
          legend.text = element_text(size = 15, colour = "black"),
          legend.title = element_blank(),
          title = element_text(size = 15)
    )
  return(mythem)
}