# Removed due to sendplot being removed from CRAN
# # Uses sendplot to annotate a cigraph object:
# cigraphAnnotate <- function(
#   graph, # graph object from cigraph
#   ..., # Named vectors indicating elements of the tooltip
#   fromcigraph = c("labels","nodeNames","tooltips","groups"), # Vector indicating which info should be extracted from cigraph object and plotted.
#   filename = "cigraph",
#   image.size = "600x600", 
#   window.size = image.size,
#   legend = FALSE # Overwries legend plotting
#   )
# {
#   if(!requireNamespace("sendplot")) stop("'sendplot' package needs to be installed.")
#   
#   # List containing the labels:
#   TooltipContents <- list(...)
# 
#   # Extract info from cigraph:
#   if ("labels" %in% fromcigraph && !is.null(graph$graphAttributes$Nodes$labels) && !is.logical(graph$graphAttributes$Nodes$labels))
#   {
#     TooltipContents$Label <- graph$graphAttributes$Nodes$labels
#   }
# 
#   if ("nodeNames" %in% fromcigraph && !is.null(graph$graphAttributes$Nodes$names) && !is.logical(graph$graphAttributes$Nodes$names))
#   {
#     TooltipContents$Name <- graph$graphAttributes$Nodes$names
#   }
#   
#   if ("tooltips" %in% fromcigraph && !is.null(graph$graphAttributes$Nodes$tooltips) && !is.logical(graph$graphAttributes$Nodes$tooltips))
#   {
#     TooltipContents$Tooltip <- graph$graphAttributes$Nodes$tooltips
#   }
#   
#   if ("groups" %in% fromcigraph && !is.null(graph$graphAttributes$Graph$groups) && length(graph$graphAttributes$Graph$groups) > 1)
#   {
#     gr <- graph$graphAttributes$Graph$groups
#     if (is.null(names(gr))) names(gr) <- paste("Group",seq_along(gr))
#         
#     TooltipContents$Group <- sapply(seq_len(graph$graphAttributes$Graph$nNodes), function(n)  paste(names(gr)[sapply(gr,function(g)n%in%g)], collapse = "; "))
#   }
#   
#   TooltipContents <- as.data.frame(TooltipContents)
#   
#   # Fix for legend:
#   graph$plotOptions$legend <- legend
# 
#   # Create plot:
# #   xy.send(paste0("cigraph:::plot.cigraph(",dput(graph),")"),
#   save(graph, file = tempfile(fileext = ".RData") -> gObj)
#   if (grepl("(windows)|(ming)",R.Version()$os,ignore.case=TRUE)){
#     gObj <- gsub("\\\\","\\\\\\\\",gObj)
#   }
# 
#   if (NROW(TooltipContents) > 0)
#   {
#     sendplot::xy.send(paste0("load('",gObj,"');cigraph:::plot.cigraph(graph)"),
#             x.pos = graph$layout[,1],
#             y.pos = graph$layout[,2],
#             xy.labels = TooltipContents,
#             fname.root = filename,
#             dir = paste0(getwd(),"/"),
#             image.size = image.size,
#             window.size = window.size)
#   } else {
#     sendplot::xy.send(paste0("load('",gObj,"');cigraph:::plot.cigraph(graph)"),
#             x.pos = -100,
#             y.pos = -100,
#             xy.labels = data.frame(` ` = ''),
#             fname.root = filename,
#             dir = paste0(getwd(),"/"),
#             image.size = image.size,
#             window.size = window.size)
#   }
# #     xy.send("plot.cigraph(graph)",
# 
#   
#   return(paste0(filename,".html"))
# }