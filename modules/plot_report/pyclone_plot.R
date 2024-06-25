library(dplyr)
library(ggplot2)
library(patchwork)

my_ggplot_theme <- function (cex = 1) 
{
  cex_opt = 1
  theme_light(base_size = 10 * cex_opt) + theme(legend.position = "bottom", 
                                                legend.key.size = unit(0.3 * cex_opt, "cm"), panel.background = element_rect(fill = "white"))
}
pyclone_smart_colors <- function (x, pl, colors)
{
  wh_col = unique(x$cluster_id)
  wh_col_missing = !(wh_col %in% names(colors))
  wh_col = wh_col[wh_col_missing]
  mycolors = colors
  new_col = NULL
  if (length(wh_col) < 9) {
    new_col = suppressWarnings(RColorBrewer::brewer.pal(length(wh_col), 
                                                        "Set1"))
  }else {new_col = rainbow(length(wh_col))}
  names(new_col) = sort(wh_col)
  return(c(mycolors, new_col[!is.na(names(new_col))]))
}

### MIXING PROP ###
pyclone_mixing_proportion <- function (x, colors = c(Tail = "gainsboro")) 
{
  mixing_prop <- x %>% dplyr::count(cluster_id) %>% 
    dplyr::mutate(cluster_id=as.character(cluster_id)) %>% 
    dplyr::mutate(prop=round(n/sum(n),2))
  p = ggplot(mixing_prop, aes(x = cluster_id, y = prop, fill = cluster_id)) + 
    geom_bar(stat = "identity") + my_ggplot_theme(cex) + 
    ylim(0, 1) + labs(title = bquote(bold("Mixing proportions"))) + 
    scale_color_manual(values = pyclone_smart_colors(x, p, colors))
  return(p)
  # mobster:::add_fill_color_pl(x, pl, colors)
}

### ELBO ###
pyclone_ELBO <- function (h5_file, cex = 1) 
{
  file_path <- h5_file
  # h5ls(file_path)
  # data <- h5read(file_path, "data")
  stats <- rhdf5::h5read(file_path, "stats")
  # stats$elbo
  # stopifnot(inherits(x, "vb_bmm"))
  ELBO = tidyr::tibble(step = 1:length(stats$elbo), ELBO = stats$elbo)
  plt <- ggplot(ELBO, aes(step, ELBO)) + geom_line(color = "steelblue") + 
    geom_point(size = 1.5 * cex) + #+ mobster:::my_ggplot_theme(cex) +
    labs(title = bquote(bold("ELBO"))) + my_ggplot_theme(cex)
}


##################### MARGINAL UNIVARIATE ###################### 
# x is the found in CNAqc2tsv/TEST/MSeq_Set06/*/joint_table.tsv
# y is the best_fit.tsv from pylone
###############################################################
pyclone_plot_1D <- function (x, y, colors = NA) 
{
  patient_id <- data.frame(do.call('rbind', strsplit(as.character(y$mutation_id),':',fixed=TRUE))) %>% 
    dplyr::select(X1) %>% 
    unique()
  x <- x %>% 
    dplyr::mutate(mutation_id=paste0(patient_id$X1,":",chr,":",from,":",alt)) %>% 
    dplyr::rename("sample_id"=Indiv)
  x$mutation_id <- gsub(pattern = "chr",replacement = "",x = x$mutation_id)
  
  joined <- dplyr::inner_join(x = x,y = y, by =c("mutation_id","sample_id"))
  sample_ids = joined$sample_id %>% unique()
  F_data_n <- joined %>% 
    dplyr::mutate(cluster_id=paste0("C",cluster_id)) %>% 
    dplyr::select(sample_id,VAF,cluster_id)
  ns = sample_ids %>% length()
  by_row = ceiling(sqrt(ns * (ns - 1)/2))
  if (by_row > sample_ids %>% length()) 
    by_row = sample_ids %>% length()
  myp = ggplot(F_data_n) + geom_histogram(aes(x = VAF, fill = cluster_id), 
                                          binwidth = 0.01) + facet_wrap(~sample_id) + 
    xlim(0.01, 1.01) + guides(fill = guide_legend("Cluster")) + my_ggplot_theme(cex)
  if (!all(is.na(colors))) 
    myp = myp + scale_fill_manual(values = colors)
  return(myp)
}

######################## MULTIVARIATE ########################
# x is the found in CNAqc2tsv/TEST/MSeq_Set06/*/joint_table.tsv
# y is the best_fit.tsv from pylone
###############################################################
pyclone_plot_2D <- function (x, y, d1, d2, cex = 1, alpha = 0.3, cut_zeroes = TRUE)
{
  caption = ""
  patient_id <- data.frame(do.call('rbind', strsplit(as.character(y$mutation_id),':',fixed=TRUE))) %>% 
    dplyr::select(X1) %>% 
    unique()
  x <- x %>% 
    dplyr::mutate(mutation_id=paste0(patient_id$X1,":",chr,":",from,":",alt)) %>% 
    dplyr::rename("sample_id"=Indiv)
  x$mutation_id <- gsub(pattern = "chr",replacement = "",x = x$mutation_id)
  joined <- dplyr::inner_join(x = x,y = y, by =c("mutation_id","sample_id"))
  data <- joined %>% dplyr::select(mutation_id,sample_id, VAF, cluster_id) %>% 
    dplyr::arrange(mutation_id) %>% 
    dplyr::mutate(VAF=as.numeric(VAF)) %>% 
    tidyr::pivot_wider(values_from = VAF,names_from = sample_id)
  
  p = ggplot(data = data, aes(x = eval(parse(text = d1)), y = eval(parse(text = d2)), colour = factor(cluster_id))) + 
    geom_point(alpha = alpha, size = 1 * cex) + 
    labs(title = bquote(bold(.(d1)) ~ "vs" ~ bold(.(d2))), 
         caption = caption, 
         x = d1, 
         y = d2) + 
    guides(color = guide_legend(title = "Cluster", 
                                override.aes = list(alpha = 1))) + 
    theme(legend.position = "bottom", 
          legend.key.size = unit(0.3 * cex, "cm"), 
          legend.text = element_text(size = 8 * cex)) + 
    geom_vline(xintercept = 0, colour = "darkgray", size = 0.3) + 
    geom_hline(yintercept = 0, colour = "darkgray", size = 0.3) + 
    guides(fill = "none") + my_ggplot_theme(cex)
  return(p)
}



pyclone_cluster_peaks <- function (x, y,cex = 1, colors = NA) 
{
  patient_id <- data.frame(do.call('rbind', strsplit(as.character(y$mutation_id), ':', fixed=TRUE))) %>% 
    dplyr::select(X1) %>% 
    unique()
  x <- x %>% 
    dplyr::mutate(mutation_id=paste0(patient_id$X1,":",chr,":",from,":",alt)) %>% 
    dplyr::rename("sample_id"=Indiv)
  x$mutation_id <- gsub(pattern = "chr",replacement = "",x = x$mutation_id)
  joined <- dplyr::inner_join(x = x,y = y, by =c("mutation_id","sample_id"))
  
  peaks <- joined %>% 
    dplyr::mutate(cluster_id=as.character(cluster_id)) %>% 
    group_by(cluster_id, sample_id) %>% 
    summarise(assignement_prob=mean(VAF), .groups = 'drop')
  colnames(peaks) <- c("cluster_id","Dimension","assignement_prob")
  # > peaks
  # Var1 Var2      value
  # 1 Set7_57   C1 0.01812786
  # 2 Set7_57   C2 0.45048328
  # 
  # colnames(peaks) = c("cluster_id", "assignement_prob","Dimension")
  # prop <- dplyr::inner_join(mixing_prop,peaks,"cluster_id")
  # prop$sample_id <-
  # prop$proportion = paste0("Mix. prop. ", prop$prop, 
  #                          "%")
  p = ggplot(peaks, aes(x = Dimension, y = assignement_prob, ymax = assignement_prob, ymin = 0, color = cluster_id)) + 
    geom_linerange() + 
    geom_point() + 
    my_ggplot_theme(cex) +
    facet_wrap(~cluster_id, nrow = 1) + 
    ylim(0, 1) +
    labs(title = bquote(bold("Binomial peaks")))
  # geom_text(data = prop, aes(label = proportion, x = 1, 
  #                            y = 1), inherit.aes = FALSE, hjust = 0, size = 2.5 * 
  #             cex) 
  return(p)
}


plot_summary_pyclone <- function(x, y, h5_file, d1, d2 = NULL, cex = 1, alpha = 0.3, cut_zeroes = TRUE)
{
  n_samples <- 2
  marginals <- pyclone_plot_1D(x = x, y = y, colors = NA)
  if (is.null(d2)){
    multivariate <- ggplot()
  } else{
    multivariate <- pyclone_plot_2D(x = x, y = y, d1 = d1, d2 = d2,cex = cex,alpha = 0.3, cut_zeroes = TRUE)
  }
  top_p = patchwork::wrap_plots(marginals, multivariate, design=ifelse(n_samples>2, "A\nB\nB", "AAB"))
  
  elbo <- pyclone_ELBO(h5_file)
  mix_p <- pyclone_mixing_proportion(x = y)
  binom <- pyclone_cluster_peaks(x, y,cex)
  bottom_p = patchwork::wrap_plots(mix_p, binom, elbo, design="ABBBC")
  report_fig = patchwork::wrap_plots(top_p, bottom_p, design=ifelse(n_samples>2, "A\nA\nA\nB", "A\nA\nB"))
  return(report_fig)
}

