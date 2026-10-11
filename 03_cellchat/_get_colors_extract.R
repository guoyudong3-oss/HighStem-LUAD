# Shared color map across Fig3A/B/C so node/bar/scatter colors are consistent.
get_fig3_group_colors <- function(groups) {
  base_colors <- c(
    "B lymphocytes" = "#E41A1C",
    "DTStem" = "#377EB8",
    "Endothelial cells" = "#4DAF4A",
    "Fibroblasts" = "#984EA3",
    "HighStem" = "#F18F01",
    "LowStem" = "#F781BF",
    "MAST cells" = "#B39BC8",
    "Myeloid cells" = "#A65628",
    "NK cells" = "#56B4E9",
    "NormalEpi" = "#1F2A7C",
    "T lymphocytes" = "#1B9E77",
    "T/NK cells" = "#9BCF6E",
    "Undetermined" = "#D9B300"
  )
  miss <- setdiff(groups, names(base_colors))
  if (length(miss) > 0) {
    extra <- grDevices::hcl.colors(length(miss), palette = "Dynamic")
    names(extra) <- miss
    base_colors <- c(base_colors, extra)
  }
  base_colors[groups]
}
