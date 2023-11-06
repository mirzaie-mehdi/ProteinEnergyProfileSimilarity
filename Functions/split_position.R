split_position <- function(position){# position means the chain and residue numbers
  input_string <- position
  split_strings <- unlist(strsplit(input_string, ","))# split the string by comma
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
