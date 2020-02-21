# $$$$$$$\  $$$$$$$\   $$$$$$\  $$$$$$$$\ $$$$$$\ $$\       $$$$$$$$\  $$$$$$\
# $$  __$$\ $$  __$$\ $$  __$$\ $$  _____|\_$$  _|$$ |      $$  _____|$$  __$$\
# $$ |  $$ |$$ |  $$ |$$ /  $$ |$$ |        $$ |  $$ |      $$ |      $$ /  \__|
# $$$$$$$  |$$$$$$$  |$$ |  $$ |$$$$$\      $$ |  $$ |      $$$$$\    \$$$$$$\
# $$  ____/ $$  __$$< $$ |  $$ |$$  __|     $$ |  $$ |      $$  __|    \____$$\
# $$ |      $$ |  $$ |$$ |  $$ |$$ |        $$ |  $$ |      $$ |      $$\   $$ |
# $$ |      $$ |  $$ | $$$$$$  |$$ |      $$$$$$\ $$$$$$$$\ $$$$$$$$\ \$$$$$$  |
# \__|      \__|  \__| \______/ \__|      \______|\________|\________| \______/




#' Download cosmic signature set
#'
#' It downloads the matrix with mutation signatures profiles
#' from cosmic. Version 2 and 3 are available.
#'
#' For v3, signatures are stored in \url{synapse}{https://www.synapse.org/}.
#' Therefore, to download such signatures you should include an apiKey and
#' your email. To obtain the apiKey, log in, go to profile, settings and
#' API key
#'
#' @param version Either v2 or v3 indicating the version of the signatures
#' @param synapse_email Better use SYNAPSE_EMAIL environtment variable
#' @param synapse_apiKey Better use SYNAPSE_PAT environtment variable
#'
#' @return a matrix with selected set of mutational signatures from cosmic
#' @export
#'
#' @references
#'
#' Go to the following link for more information: \url{https://cancer.sanger.ac.uk/cosmic/signatures}
#'
download_signature_set <- function(version = c("v2","v3"),
                                   synapse_email = Sys.getenv("SYNAPSE_EMAIL"),
                                   synapse_apiKey = Sys.getenv("SYNAPSE_PAT")){

  if (length(version)>1){
    version = version[1] # defaults to "cosmic"
  }

  if (version == "v3" & is.null(synapse_apiKey)){
    stop("The v3 is stored in Synapse and a apiKey is required")
    if (!requireNamespace("synapser")) {
      stop("You need to install synapser to use this function.")
    }
  }

  # URL / CODES
  v2_url = "https://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt"

  switch (version,
    v2 = {
      dat = suppressMessages(suppressWarnings({
        readr::read_tsv(file = v2_url)
        }))
      sigs = as.matrix(dat[,4:33])
      k1 = dat$Trinucleotide
      alt = stringr::str_extract(pattern = "(?<=\\>)[ACTG]",
                                 string = dat$`Substitution Type`)
      rownames(sigs) = glue::glue("{k1}>{alt}")
      sigs
    },
    v3 = {
      synapser::synLogin(email = synapse_email,
                         apiKey = synapse_apiKey)
      # Obtain a pointer and download the data
      syn11738319 <- synapser::synGet(entity='syn11738319')

      dat = suppressMessages(suppressWarnings({
        readr::read_csv(file = syn11738319$path)
      }))

      sigs = as.matrix(dat[,3:ncol(dat)])
      k1 = dat$SubType
      alt = stringr::str_extract(pattern = "(?<=\\>)[ACTG]",
                                 string = dat$Type)
      rownames(sigs) = glue::glue("{k1}>{alt}")
      sigs
    }
  ) -> res_sigs

  res_sigs
}



#' Compare 2 sets of signatures
#'
#' Each signature of the 2 matrices is compared against each of the
#' signatures of the other matrix using cosine similarity.
#'
#'
#' @param x A matrix with signatures as columns and rows with mutation names
#' @param y A matrix with signatures as columns and rows with mutation names
#' @param simplify_set Reference bases to simplify the mutation names
#' @param sep separation between trinucleotide and alternative allele
#'
#' @return Returns a matrix with pair-wise cosine similarity between signatures
#' @export
#'
compare_signature_sets <- function(x,
                                   y,
                                   simplify_set = c("C","A"),
                                   sep = ">") {
  ## x and y need to be 2 matrix with signatures as columns
  ## rownames should contain the mut names
  rownames(x) = simplify_muts(rownames(x),
                              simplify_set = simplify_set,
                              sep = sep)
  rownames(y) = simplify_muts(rownames(y),
                              simplify_set = simplify_set,
                              sep = sep)

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



#' Cosine similarity between 2 named vectors
#'
#' It reorder the given vectors according to the names in x and performs a
#' cosine similarity computation.
#'
#' @param x A vector (optionally named)
#' @param y A vector (optionally named)
#'
#' @return A single double indicating the cosine distant
#' @export
#'
#' @references
#'
#' See here \url{https://en.wikipedia.org/wiki/Cosine_similarity}
#'
#' @seealso
#'
#' Similar function based on matrices available at \link[lsa]{cosine}
#'
#' @examples
#'
#' cos_sim_vect(x = c(1,2,3),c(2,3,4))
#'
cos_sim_vect = function (x, y) {

  if (! is.vector(x) || ! is.vector(y) ) {
    stop("Inputs are not vectors")
  }

  stopifnot(length(x) == length(y))

  # I can put doble because if the first is null then no nomes are considered.
  if(is.null(names(x)) || is.null(names(y))) {
    warning("No names considered, using given order")
  } else{
    # reorder the vectors based on names
    base_names=names(x)
    y = y[base_names]
  }

  # compute the cosine sim
  return( ( crossprod(x, y)/sqrt(crossprod(x) * crossprod(y)) ) [1,1] )
}

