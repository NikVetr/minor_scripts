library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

dot_path <- "~/scripts/minor_scripts/postdoc/eGTEx-ASE_methods-workflow_digraph.txt"
dot_lines <- paste(readLines(dot_path), collapse = "\n")

graph_svg <- export_svg(grViz(dot_lines))
writeLines(graph_svg, "~/eGTEx-ASE_methods-workflow.svg")
rsvg_pdf("~/eGTEx-ASE_methods-workflow.svg", "~/eGTEx-ASE_methods-workflow.pdf")
