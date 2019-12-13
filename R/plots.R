## plots


mutMB_preplot_VR <- function(vr,pivot_vector,total_mbases){
  
  #browser()
  stopifnot(.lenunique(sampleNames(vr))== length(pivot_vector))
  
  vr$project = pivot_vector[as.character(sampleNames(vr))]
  vr$sample = sampleNames(vr)
  
  tm = mcols(vr) %>% data.frame() %>% dplyr::group_by(sample,project) %>% 
    dplyr::summarise(n = n())
  
  tm = tm %>% dplyr::ungroup() %>% dplyr::group_by(project) %>% 
    dplyr::mutate(n_samples = n()) %>% 
    dplyr::mutate(project_label = glue::glue("{project} (n={n_samples})"))
  
  tm$mutMB = tm$n / total_mbases
  
  
  return(tm)
}


plot_mutMB <- function(vr,
                       total_mbases = 2000,
                       main = "",
                       pivot_vector = NULL,
                       cutoffs = c(1,10,100,1000)) {
  library(cowplot)
  theme_set(theme_cowplot(font_size = 10))
  library(forcats)
  
  #browser()
  if (is.null(pivot_vector)){
    pivot_vector = rep("all",.lenunique(sampleNames(vr)))
    names(pivot_vector) = unique(sampleNames(vr))
  }

  seqlevelsStyle(vr) <- "UCSC"
  
  tm = mutMB_preplot_VR(vr = vr,
                        pivot_vector = pivot_vector,
                        total_mbases = total_mbases)
  
  break_points = cutoffs
  type_plot = total_mbases
  
  
  tm$bp = findInterval(x = tm$mutMB,vec = break_points)
  
  tm$bp = break_points[tm$bp +1 ]
  bp = tm %>% dplyr::ungroup() %>% dplyr::group_by(bp,project_label) %>% 
    dplyr::summarise(n = n()) %>% dplyr::ungroup() %>% 
    dplyr::group_by(project_label) %>% 
    dplyr::mutate(per=scales::percent(n/sum(n)))
  #browser()
  
  
  fp = ggplot(tm,aes(x =fct_reorder(sample,-n,max),y = mutMB)) + 
    geom_point(aes(color=factor(bp))) +
    geom_text(inherit.aes = FALSE,
              data = bp,
              aes(x = Inf,
                  y=bp,
                  label = glue::glue("{per}")),
              hjust = 1.5,
              vjust = 1.1) +
    scale_y_log10(breaks=break_points) +
    facet_grid(.~project_label, scales="free",space = "free") +
    theme(axis.text.x = element_blank(),
          panel.grid.major.y = element_line(color = "gray",
                                            linetype="dashed"),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          legend.position = "none",
          axis.line.x = element_blank(),
          strip.background = element_rect(fill = NA)) +
    scale_color_brewer(palette = "Set1") +
    labs(y = glue::glue("Mutations / Mbp"),
         caption = glue::glue("{type_plot} Mbp used as estimation of length"),
         title = main)
  
  return(fp)

}

preprocess_rainPlot <- function(vr){
  stopifnot(requireNamespace("VariantAnnotation",
                             "glue",quietly = TRUE))
  stopifnot(.lenunique(VariantAnnotation::sampleNames(vr))==1)
  
  dat = vr
  
  dist = GenomicRanges::distanceToNearest(dat)
  
  dat$dist = -1
  GenomicRanges::mcols(dat[S4Vectors::queryHits(dist)])$dist = GenomicRanges::mcols(dist)$distance
  no_pair = dat$dist < 0
  dat = dat[!no_pair]
  warning(glue::glue("Removed {sum(no_pair)} muts because no pair was found."))
  
  return(dat)
}

plot_rainPlot <- function(vr,main = NULL) {
  
  stopifnot(requireNamespace("VariantAnnotation",quietly = TRUE))
  stopifnot(requireNamespace("glue",quietly = TRUE))
  stopifnot(requireNamespace("ggplot2",quietly = TRUE))
  library(ggplot2)
  
  #GenomeInfoDb::seqlevelsStyle(vr) = "UCSC"
  dat = preprocess_rainPlot(vr)
  if(is.null(main) ){main = unique(VariantAnnotation::sampleNames(vr))}
  
  mut_type = simplify_muts(paste(VariantAnnotation::ref(dat),
                                 VariantAnnotation::alt(dat),
                                 sep = ">"))
  
  # this should mantain the original order of the seqnames provided 
  # by the VR object
  sqNames = factor(as.character(GenomicRanges::seqnames(dat)),
                   levels = seqlevels(dat))
  
  df = tibble::tibble(
    start = GenomicRanges::start(dat),
    dist = dat$dist,
    seqnames = sqNames,
    mutation_type = mut_type
  )
  
  #browser()
  mMb = (nrow(df) / (3000* 2/3)) %>% round(digits = 1)
  # approximation counting alignability
  caption_text = glue::glue("{mMb} muts/mB (WGS)
                            Some random noise added to the raw values to help pair-mutation visualization")
  
  #df$seqnames = factor(df$seqnames,levels = chromosomes_UCSC_in)
  
  rainPlot = df %>%
    ggplot(aes(x = start,
               y=dist)) +
    facet_grid(.~seqnames,
               scales = "free",
               space="free_x",
               switch = "x") +
    geom_point(position=position_jitter(h=0.1,w=10),
               aes(color = mutation_type)) +
    scale_y_log10(labels = scales::comma) +
    scale_color_manual(values = tr_colors) +
    geom_hline(yintercept = 1000,
               linetype = "dashed") +
    geom_smooth(inherit.aes = F,
                method = "lm",
                mapping = aes(x = start,y=dist),
                se = F,
                color = "black") +
    theme_bw() +
    theme(legend.position = "none",
          axis.line.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.spacing = unit(0, "lines"),
          strip.text = element_text(angle = 90),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.border = element_rect(color = "gray50", fill = NA)) +
    labs(
      caption = caption_text,
      x = "Position",
      y = "Distance to nearest (bp)",
      title = main
    )
  
  return(rainPlot)
}
  

