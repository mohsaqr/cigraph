
cigraph.gui <- function(input,corMat,...)
{
  
  stop("This function has been removed in favor of our Shiny app: https://jolandakos.shinyapps.io/NetworkApp/")
#   
# 	if (!require("rpanel")) stop("Package 'rpanel' is required to use GUI functionality")
# 
#   ## CHECK FOR CORRELATION MATRIX ###:
#   if (missing(corMat))
#   {
#     if (is.matrix(input) && isSymmetric(input) && all(diag(input)==1) && all(abs(input)<=1))
#     {
#       corMat <- TRUE
#     } else corMat <- FALSE
#   }
#   
#   if (!is.logical(corMat)) stop("'corMat' must be logical")
#   
# 	if (any(grepl("RStudio", .libPaths(), ignore.case=TRUE)))
# 	{
# 	  if (grepl("win",Sys.info()["sysname"],ignore.case=TRUE))
# 	  {
# 	    windows()
# 	  } else X11()
# 	} else
# 	{
# 	  dev.new()
# 	}
#   # Dummies to fool R CMD check:
#   graph <- minimum <- maximum <- esize <- vsize <- asize <- graph <- cbox <- filename <- dimensions <- OnTheFly <- LatSize <- GraphType <- FAopts <- File <- Control <- NULL
# 	
#   ### Correlation matrix GUI:
#   if (corMat)
#   {
#     covMat <- input
#     if (any(diag(input)!=1)) input <- round(cov2cor(input),12)
#     
#     cigraph.setup <- function(panel) {
#       if (isTRUE(panel$OnTheFly)) cigraph.draw(panel)
#       panel
#     }
#     
#     cigraph.draw <- function(panel) {
#       panel$details <- panel$cbox[1]
#   	panel$bg <- panel$transparency <-panel$cbox[2]
#   	panel$overlay <- panel$cbox[3]
#   	panel$borders <- panel$cbox[4]
#   	panel$legend <- panel$cbox[5]
#   # 	panel$width <- as.numeric(panel$dimensions[1])
#   # 	panel$height <- as.numeric(panel$dimensions[2])
#       if (panel$GraphType=="EFA")
#       {
#         panel2 <- panel
#         panel2$vsize <- c(panel$vsize,panel$LatSize)
#         panel2$dat <- covMat
#         panel2$factors <- as.numeric(panel$FAopts[1])
#         panel2$rotation <- panel$FAopts[2]
#         do.call(cigraph.efa,panel2)
#         panel2$graph <- NULL
#       } else if (panel$GraphType == "PCA")
#       {
#         panel2 <- panel
#         panel2$vsize <- c(panel$vsize,panel$LatSize)
#         panel2$cor <- covMat
#         panel2$factors <- as.numeric(panel$FAopts[1])
#         panel2$rotation <- panel$FAopts[2]
#         panel2$graph <- NULL
#         do.call(cigraph.pca,panel2)
#       } else
#       {
#         panel$graph <- panel$GraphType
#         do.call(cigraph,panel)
#       }
#       panel
#     }
#     cigraph.newplot <- function(panel)
#     {
#       if (any(grepl("RStudio", .libPaths(), ignore.case=TRUE)))
#       {
#         if (grepl("win",Sys.info()["sysname"],ignore.case=TRUE))
#         {
#           windows()
#         } else X11()
#       } else
#       {
#         dev.new()
#       }
#   	cigraph.draw(panel)
#   	panel
#     }
#     
#     cigraph.save <- function(panel) {
#       panel2 <- panel
#       panel2$width <-  par("din")[1]
#       panel2$height <- par("din")[2]
#       panel2$filename <- panel$File[1]
#       panel2$filetype <- panel$File[2]
#       cigraph.draw(panel2)
#       panel
#     }
#     
#     cigraph.panel <- rp.control("cigraph GUI", input = input, ...)
#     
#     rp.checkbox(cigraph.panel, OnTheFly ,cigraph.setup, title="Plot on the fly", pos = list(column=0,row=0), initval=FALSE)
#     
#     rp.slider(cigraph.panel, minimum, 0, 1 , cigraph.setup, "Minimum", 	initval = 0, showvalue = TRUE, pos = list(column=0,row=1))
#     rp.slider(cigraph.panel, cut, 0, 1 , cigraph.setup, "Cutoff", initval = 0.4, showvalue = TRUE, pos = list(column=0,row=2))
#     rp.slider(cigraph.panel, maximum, 0, 1 , cigraph.setup, "Maximum", initval = 1, showvalue = TRUE, pos = list(column=0,row=3))
#   
#      rp.slider(cigraph.panel, esize, 0, 20 , cigraph.setup, "Edge width", initval = 4, showvalue = TRUE,  pos = list(column=1,row=1))
#     rp.slider(cigraph.panel, vsize, 0, 20, cigraph.setup, "Node size", initval = 2, showvalue = TRUE, pos = list(column=1,row=2))
#     rp.slider(cigraph.panel, LatSize, 1, 20 , cigraph.setup, "Latent size (FA)", initval = 5, showvalue = TRUE, pos = list(column=1,row=3))
#     
#     rp.radiogroup(cigraph.panel, GraphType, c("association", "concentration", "factorial", "EFA", "PCA"), title = "Graph", action = cigraph.setup, pos = list(column=0,row=4))
#     
#     
#   	rp.textentry(cigraph.panel, FAopts, cigraph.setup, initval = c("1","promax"),  pos = list(column=1,row=4),title="FA Options", labels = c("# Factors","Rotation"))  
#     
#     
#   	rp.radiogroup(cigraph.panel, layout, c("circular", "spring", "tree"), title = "Layout", action = cigraph.setup, pos = list(column=0,row=5))
#   
#   	rp.checkbox(cigraph.panel, cbox,cigraph.setup, labels = c("Details","Background", "Overlay","Borders","Legend"), title="Options", pos = list(column=1,row=5), initval=c(FALSE,FALSE,FALSE,TRUE,TRUE))
#   
#   	rp.textentry(cigraph.panel, File, cigraph.setup, initval = c("cigraph","pdf"),  pos = list(column=1,row=6),title="Output", labels = c("Name","Type"))
#   	
#   # 	rp.textentry(cigraph.panel, dimensions, cigraph.setup, initval = c(7,7),  pos = list(column=1,row=4),labels = c("Width","Height"),title="Dimensions (enter to confirm)")
#   	
#   	rp.button(cigraph.panel, action = cigraph.draw, title = "Plot", ,pos = list(column=0,row=7))	
#   	
#   	rp.button(cigraph.panel, action = cigraph.newplot, title = "New" ,pos = list(column=1,row=7))	
#   
#   	rp.button(cigraph.panel, action = cigraph.save, title = "Save", pos = list(column=0,row=6))
#     
#   } else {
#     # Default GUI:
#     
#     cigraph.setup <- function(panel) {
#       if (isTRUE(panel$OnTheFly)) cigraph.draw(panel)
#       panel
#     }
#     
#     cigraph.draw <- function(panel) {
#       if (panel$Control[1]!="") {
#         panel$minimum <- as.numeric(panel$Control[1])
#       } else {
#         panel$minimum <- NULL
#       }
#       
#       if (panel$Control[2]!="") {
#         panel$cut <- as.numeric(panel$Control[1])
#       } else {
#         panel$cut <- NULL
#       }      
#       
#       if (panel$Control[3]!="") {
#         panel$maximum <- as.numeric(panel$Control[3])
#       } else {
#         panel$maximum <- NULL
#       }
#       
#       if (panel$Control[4]!="") {
#         panel$esize <- as.numeric(panel$Control[4])
#       } else {
#         panel$esize <- NULL
#       }
#       
#       if (panel$Control[5]!="") {
#         panel$vsize <- as.numeric(panel$Control[5])
#       } else {
#         panel$vsize <- NULL
#       }
#       
#       if (panel$Control[6]!="") {
#         panel$asize <- as.numeric(panel$Control[6])
#       } else {
#         panel$asize <- NULL
#       }
#       
#       panel$details <- panel$cbox[1]
#       panel$bg <- panel$transparency <-panel$cbox[2]
#       panel$overlay <- panel$cbox[3]
#       panel$borders <- panel$cbox[4]
#       panel$legend <- panel$cbox[5]
#       do.call(cigraph,panel)
#       panel
#     }
#     cigraph.newplot <- function(panel)
#     {
#       if (any(grepl("RStudio", .libPaths(), ignore.case=TRUE)))
#       {
#         if (grepl("win",Sys.info()["sysname"],ignore.case=TRUE))
#         {
#           windows()
#         } else X11()
#       } else
#       {
#         dev.new()
#       }
#       cigraph.draw(panel)
#       panel
#     }
#     
#     cigraph.save <- function(panel) {
#       panel2 <- panel
#       panel2$width <-  par("din")[1]
#       panel2$height <- par("din")[2]
#       panel2$filename <- panel$File[1]
#       panel2$filetype <- panel$File[2]
#       cigraph.draw(panel2)
#       panel
#     }
#     
#     cigraph.panel <- rp.control("cigraph GUI", input = input, ...)
#     
#     rp.checkbox(cigraph.panel, OnTheFly ,cigraph.setup, title="Plot on the fly", pos = list(column=0,row=0), initval=FALSE)
#     
#     rp.textentry(cigraph.panel, Control, cigraph.setup, initval = rep("",6),  pos = list(column=1,row=1),title="Control", labels = c("Minimum","Cut","Maximum","Edge Width", "Node Size", "Arrow Size"))
#     
#     
#     rp.checkbox(cigraph.panel, cbox,cigraph.setup, labels = c("Details","Background", "Overlay","Borders","Legend"), title="Options", pos = list(column=0,row=1), initval=c(FALSE,FALSE,FALSE,TRUE,TRUE))
#     
#     rp.radiogroup(cigraph.panel, layout, c("circular", "spring"), title = "Layout", action = cigraph.setup, pos = list(column=1,row=2))
#     
#     
#     rp.textentry(cigraph.panel, File, cigraph.setup, initval = c("cigraph","pdf"),  pos = list(column=1,row=3),title="Output", labels = c("Name","Type"))
#     
#     # 	rp.textentry(cigraph.panel, dimensions, cigraph.setup, initval = c(7,7),  pos = list(column=1,row=4),labels = c("Width","Height"),title="Dimensions (enter to confirm)")
#     
#     rp.button(cigraph.panel, action = cigraph.draw, title = "Plot", ,pos = list(column=0,row=4))	
#     
#     rp.button(cigraph.panel, action = cigraph.newplot, title = "New" ,pos = list(column=1,row=4))	
#     
#     rp.button(cigraph.panel, action = cigraph.save, title = "Save", pos = list(column=0,row=3))
#   }
}
