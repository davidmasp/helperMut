

download_signature_set <- function(type = c("cosmic","pcawg"),
                                   synapse_token = Sys.getenv("SYNAPSE_PAT")){
  if (length(type)>1){
    type = type[1] # defaults to "cosmic"
  }

  # URL / CODES
  cosmic_url = "https://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt"
  switch (type,
    cosmic = {
      dat = readr::read_tsv(file = cosmic_url)
      sigs = as.matrix(dat[,4:33])
      k1 = dat$Trinucleotide
      alt = stringr::str_extract(pattern = "(?<=\\>)[ACTG]",
                                 string = dat$`Substitution Type`)
      rownames(sigs) = glue::glue("{k1}>{alt}")
      sigs
    }
  ) -> res_sigs
  res_sigs
}

compare_signature_sets <- function(x,y,simplify_set = c("C","A"),sep = ">") {
  ## x and y need to be 2 matrix with signatures as columns
  ## rownames should contain the mut names
  names_x = rownames(x)
  names_y = rownames(y)

  rownames(x) = simplify_muts(rownames(x),simplify_set = simplify_set,sep = sep)
  rownames(y) = simplify_muts(rownames(y),simplify_set = simplify_set,sep = sep)

  y = y[rownames(x),]

  res_dist = matrix(rep(NA,ncol(x)*ncol(y)),
                    nrow = ncol(x),
                    ncol = ncol(y),
                    dimnames = list(colnames(x),colnames(y)))

  for (i in 1:ncol(x)){
    for (j in 1:ncol(y)){
      cosine = cos_sim_vect(x = x[,i],y = y[,j])
      res_dist[i,j] = cosine
    }
  }

  res_dist
}



cos_sim_vect = function (x, y) {
  stopifnot(length(x) == length(y))

  # reorder the vectors based on names
  base_names=names(x)
  y = y[base_names]

  # compute the cosine sim
  if (! is.vector(x) || ! is.vector(y) ) {
    return(NA);
  } else {
    return( ( crossprod(x, y)/sqrt(crossprod(x) * crossprod(y)) ) [1,1] );
  }
}

