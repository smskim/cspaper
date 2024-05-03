DEFAULT_SEED <- 1948505650

PAPER_WIDTH <- 5

POINTS_TO_MM <- 0.3528
GG_FAMILY <- 'CM Roman'
GG_SLIDE_FAMILY <- 'Fira Sans Light'

base_modifications <- theme(
  panel.grid.major=element_blank(),
  panel.grid.minor=element_blank(),
  axis.title=element_text(size=10),
  axis.text=element_text(size=10),
  legend.title=element_text(size=10, face='plain'),
  legend.text=element_text(size=10),
  plot.title=element_text(size=10, hjust=0.5),
  panel.background=element_rect(fill='transparent', color=NA),
  plot.background=element_rect(fill='transparent', color=NA),
  legend.background=element_rect(fill='transparent', color=NA),
  legend.box.background=element_rect(fill='transparent', color=NA),
  legend.key=element_rect(fill='transparent', color=NA)
)

report_theme <- theme_bw(base_family=GG_FAMILY) + base_modifications

COLORS <- brewer.pal(8, 'Dark2')

margins_pt <- function(top=5.5, right=5.5, bottom=5.5, left=5.5) {
    theme(plot.margin=unit(c(top, right, bottom, left), 'pt'))
}

# https://github.com/tidyverse/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
extract_legend <- function(plot) {
    tmp <- ggplot_gtable(ggplot_build(plot))
    legend_index <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[legend_index]]
    return(legend)
}

text_ <- function(label, ..., parse=FALSE, slides=TRUE, poster=FALSE) {
    if (!parse) {
        label <- gsub('-', '\uad', label)
    }
    annotate('text', ..., label=label, parse=parse, family=GG_FAMILY,
             size=10 * POINTS_TO_MM)
}

arrow_ <- function(...) {
    annotate('segment', arrow=arrow(angle=15, length=unit(2, 'mm')), ...)
}

save_plot <- function(plot, base_name, dims, tag=NULL, filename=NULL) {
    if (!is.null(tag)) {
        base_name <- sprintf('%s_%s', base_name, tag)
    }
    if (is.null(filename)) {
        filename <- sprintf('build/%s_paper.pdf', base_name)
    }
    cat(sprintf('Saving to %s...\n', filename))
    dir.create('build', showWarnings=FALSE)
    dir.create(dirname(filename), showWarnings=FALSE)
    ggsave(filename, plot, width=dims[1], height=dims[2])
    embed_fonts(filename)
}

lighten <- function(color, factor){
    col <- col2rgb(color)
    col <- 255 - (255 - col) / factor
    col <- rgb(t(pmin(255, col)), maxColorValue=255)
    col
}
