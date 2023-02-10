# Taken from https://stackoverflow.com/a/54443955/15196126
library(ggplot2)
library(gtable)
library(lemon)

shift_legend <- function(p) {
  # ...
  # to grob
  gp <- ggplotGrob(p)
  facet.panels <- grep("^panel", gp[["layout"]][["name"]])
  empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]))
  empty.facet.panels <- facet.panels[empty.facet.panels]

  # establish name of empty panels
  empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
  names <- empty.facet.panels$name
  # example of names:
  # [1] "panel-3-2" "panel-3-3"

  # now we just need a simple call to reposition the legend
  reposition_legend(p, "center", panel = names)
}
