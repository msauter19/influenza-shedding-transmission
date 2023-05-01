export_plot <- function(plot, f, width, height) {
  ggsave(
    filename = paste0("../flu-shedding/figures/", f, ".png"),
    plot = plot,
    width = width,
    height = height,
    units = "in"
  )
}