theme_custom <- function(font.size = 12) {
  ggplot2::theme_bw(base_size = font.size) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(colour = "grey90"),
      panel.grid.minor = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey85", colour = "grey50"),
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
}