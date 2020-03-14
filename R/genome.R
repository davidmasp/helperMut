#  $$$$$$\  $$$$$$$$\ $$\   $$\  $$$$$$\  $$\      $$\ $$$$$$$$\
# $$  __$$\ $$  _____|$$$\  $$ |$$  __$$\ $$$\    $$$ |$$  _____|
# $$ /  \__|$$ |      $$$$\ $$ |$$ /  $$ |$$$$\  $$$$ |$$ |
# $$ |$$$$\ $$$$$\    $$ $$\$$ |$$ |  $$ |$$\$$\$$ $$ |$$$$$\
# $$ |\_$$ |$$  __|   $$ \$$$$ |$$ |  $$ |$$ \$$$  $$ |$$  __|
# $$ |  $$ |$$ |      $$ |\$$$ |$$ |  $$ |$$ |\$  /$$ |$$ |
# \$$$$$$  |$$$$$$$$\ $$ | \$$ | $$$$$$  |$$ | \_/ $$ |$$$$$$$$\
#  \______/ \________|\__|  \__| \______/ \__|     \__|\________|


# You can use BSgenome::available.genomes(TRUE)

#' Genome Selector
#'
#' @description
#' Wrapper function for a BSgenome usage. It can handle interactive downloads,
#' listing of available and installed reference sequences.
#'
#' Main function is to retrieve a genome object without loading the
#' whole package, can be done through aliases or interactive selection.
#'
#' @param alias string containing the alias of the genome (ex: Hsapiens.NCBI.GRCh38)
#' @param install boolean if genome has to be installed in the system if required.
#' @param available List available genomes in bioConductor
#' @param script_mode Don't toggle the interactive mode.
#' @param installed List installed genomes in the local machine
#'
#' @return A BSgenome object.
#'
#' @export
#'
#' @examples
#'
#' \dontrun{
#' # list installed genomes in the system
#' # choose one to save it as an object
#' genome = genome_selector(installed = T)
#'
#' # how to install a genome package
#' # then press the desired genome package and press 1 to install.
#' # after that it will be saved as an object
#' genome = genome_selector(available = T)
#'
#'
#' # select the desired genome package
#' genome = genome_selector(alias = "Hsapiens.1000genomes.hs37d5")
#' }
#'
genome_selector <- function(alias="Hsapiens.UCSC.hg19",
                            install = FALSE,
                            available=FALSE,
                            script_mode=FALSE,
                            installed = FALSE){

  stopifnot(requireNamespace("glue",quietly = TRUE))
  stopifnot(requireNamespace("BSgenome",quietly = TRUE))

  interactive_mode = interactive() & !script_mode

  all_genomes = c(BSgenome::available.genomes(),BSgenome::installed.genomes())

  if (install){
    installed = FALSE
    available = FALSE
  }

  if(installed & interactive_mode){
    ch_g = BSgenome::installed.genomes()

    if (!is.na(Sys.getenv("RSTUDIO", unset = NA))) {
      # RStudio specific code
      choice_g = utils::menu(choices = ch_g,
                             graphics = FALSE,
                             title = "Choose one of the available genomes.")
      if(choice_g == 0)return(NULL)
    } else {
      choice_g = utils::menu(choices = ch_g,
                             graphics = TRUE,
                             title = "Choose one of the available genomes.")
      if(choice_g == 0)return(NULL)
    }

    query = ch_g[choice_g]

  } else if (available & interactive_mode){

    ch_g = BSgenome::available.genomes()

    if (!is.na(Sys.getenv("RSTUDIO", unset = NA))) {
      # RStudio specific code
      choice_g = utils::menu(choices = ch_g,
                             graphics = FALSE,
                             title = "Choose one of the available genomes.")
      if(choice_g == 0)return(NULL)
    } else {
      choice_g = utils::menu(choices = ch_g,
                             graphics = TRUE,
                             title = "Choose one of the available genomes.")
      if(choice_g == 0)return(NULL)
    }
    query = ch_g[choice_g]
  } else if (available | installed){
    stop("No interactive session. Use install=TRUE to automatically install the needed package.")
  } else {
    query = glue::glue("BSgenome.{alias}")
    stopifnot(query %in% all_genomes)
  }

  # sanity check
  stopifnot(query %in% all_genomes)

  if (query %in% BSgenome::installed.genomes()){
    genome = BSgenome::getBSgenome(query)
  } else {
    message(glue::glue("{query} not installed."))

    # enter the choice if to install
    if (interactive_mode & !install){
      choice = utils::menu(c("Yes", "No"),
                           title="Do you want to install the package now?")
    } else if (install){
      choice = 1
    } else {
      choice = 0
    }

    if (choice == 1){
      BiocManager::install(query)
      genome = BSgenome::getBSgenome(query)
    } else {
      stop("Use install=TRUE to automatically install the needed package.")
    }
  }

  return(genome)
}






#' Get chunks from regions
#'
#' From a genome file or a GR object generate a chunked verision
#'
#' @param gr a GRanges object or a BSgenome object
#' @param wl the window length of the resulting chunks
#' @param unlist if result should be returned as a single GR, default true
#'
#' @return A GRanges object chunked
#' @export
#'
#' @examples
get_region_chunks <- function(gr,wl = 10000,unlist = TRUE) {

  # deciding which type of object we have
  if (is(gr,"GRanges")){
    start_vector = BiocGenerics::start(gr)
    end_vector = BiocGenerics::end(gr)
  } else if (is(gr,"BSgenome")) {
    start_vector = rep(1,length(gr))
    end_vector = GenomeInfoDb::seqlengths(gr)
  } else {
    stop("gr class not supported")
  }

  start = purrr::map2(.x = start_vector,
                      .y = end_vector,
                      .f = seq,
                      by = wl)

  end = purrr::map2(.x = start,
                    .y = end_vector,
                    function(x,max_value){
                      end = dplyr::lead(x) - 1
                      end[length(end)] = max_value
                      return(end)
                    })

  lol = list(
    as.character(GenomeInfoDb::seqnames(gr)),
    end,
    start
  )

  #browser()

  purrr::pmap(.l = lol,.f = function(name,end,start){
    GenomicRanges::GRanges(
      seqnames = as.character(name),
      ranges = IRanges::IRanges(start = start,end = end)
    )
  }) -> res

  if (unlist){
    names(res) = NULL
    # the warning s due to multiple things don't share seqnmas
    # could be solved by explicitely seting seqnames
    # I think it wouldn't be worth the effort.
    res = suppressWarnings(do.call(what = "c",args = res))
  }

  return(res)
}



