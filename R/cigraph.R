# Create cigraph model:

cigraph <- function( input, ... )
{
  
  # OTHER INPUT MODES: 
  # if (any(class(input)=="factanal") )
  # {
  #   return(cigraph.efa(input,...))
  # } else if (any(class(input)=="principal") )
  # {
  #   return(cigraph.pca(input,...))
  # } else if (any(class(input)=="lavaan"))
  # {
  #   return(cigraph.lavaan(input,edge.labels=TRUE,include=8,filetype="",...))
  # } else if (any(class(input)=="sem"))
  # {
  #   return(cigraph.sem(input,edge.labels=TRUE,include=6,filetype="",...))
  # } else 
    if (is(input,"loadings"))
  {
    return(cigraph.loadings(input,...))
  # }  
  # else if (any(class(input)=="semmod"))
  # {
  #   return(cigraph.semModel(input,...))
  } else if (is.list(input) && identical(names(input),c("Bhat", "omega", "lambda1", "lambda2")))
  {
    layout(t(1:2))
    
    Q1 <- cigraph((input$omega + t(input$omega) ) / 2,...)
    Q2 <- cigraph(input$Bhat,...)
    
    return(list(Bhat = Q1, omega = Q2))
  }
  
  
  ### EMPTY QGRAPH OBJECT ####
  cigraphObject <- list(
    Edgelist = list(),
    Arguments = list(),
    plotOptions = list(),
    graphAttributes = list(
      Nodes = list(),
      Edges = list(),
      Graph = list()
    ),
    layout = matrix(),
    layout.orig = matrix()
  )
  
  class(cigraphObject) <- "cigraph"
  

  ### Extract nested arguments ###
  # if ("cigraph"%in%class(input)) cigraphObject$Arguments <- list(...,input) else cigraphObject$Arguments <- list(...)
  cigraphObject$Arguments <- list(...,input=input) 
  
  if (isTRUE(cigraphObject$Arguments[['gui']]) | isTRUE(cigraphObject$Arguments[['GUI']])) 
  {
    cigraphObject$Arguments$gui <- cigraphObject$Arguments$GUI <- FALSE
    return(invisible(do.call(cigraph.gui,c(list(input=input),cigraphObject$Arguments))))
  }
  
  if(!is.null(cigraphObject$Arguments$adj))
  {
    stop("'adj' argument is no longer supported. Please use 'input'")
  }
  
  # Import cigraphObject$Arguments:
  if (length(cigraphObject$Arguments) > 0) cigraphObject$Arguments <- getArgs(cigraphObject$Arguments)
  
  # Import default arguments:
  def <-  getOption("cigraph")
  if (!is.null(def$cigraph)) class(def$cigraph) <- "cigraph"
  if (any(sapply(def,function(x)!is.null(x))))
  {
    cigraphObject$Arguments <- getArgs(c(cigraphObject$Arguments,def))
  }
  
  # If cigraph object is used as input, recreate edgelist input:
  if (is(input,"cigraph")) 
  {
    # if (is.null(cigraphObject$Arguments$directed)) cigraphObject$Arguments$directed <- input$Edgelist$directed
    # if (is.null(cigraphObject$Arguments$bidirectional)) cigraphObject$Arguments$bidirectional <- input$Edgelist$bidirectional
    # if (is.null(cigraphObject$Arguments$nNodes)) cigraphObject$Arguments$nNodes <- input$graphAttributes$Graph$nNodes
    # 
    if (!is.null(cigraphObject$Arguments$input)){
      input <- cigraphObject$Arguments$input
    } else {
      if(input[['graphAttributes']][['Graph']][['weighted']])
      {
        input <- cbind(input$Edgelist$from,input$Edgelist$to,input$Edgelist$weight)
      } else
      {
        input <- cbind(input$Edgelist$from,input$Edgelist$to)
      }      
    }

    # cigraphObject$Arguments$edgelist <- TRUE
  }
  
  ### PCALG AND GRAPHNEL ###
  if (is(input,"pcAlgo") | is(input,"graphNEL"))
  {
    if (is(input,"pcAlgo")) graphNEL <- input@graph else graphNEL <- input
    cigraphObject$Arguments$directed <- graphNEL@graphData$edgemode == "directed"
    cigraphObject$Arguments$bidirectional <- TRUE
    TempLabs  <- graphNEL@nodes
    if (is.null(cigraphObject$Arguments$labels)) cigraphObject$Arguments$labels  <- graphNEL@nodes
    weights <- sapply(graphNEL@edgeData@data,'[[','weight')
    
    EL <- laply(strsplit(names(weights),split="\\|"),'[',c(1,2))
    #       EL <- apply(EL,2,as.numeric)
    EL[,1] <- match(EL[,1],TempLabs)
    EL[,2] <- match(EL[,2],TempLabs)
    mode(EL) <- "numeric"
    # Create mixed graph if pcAlgo:
    if (is(input,"pcAlgo"))
    {
      srtInput <- aaply(EL,1,sort)
      cigraphObject$Arguments$directed <- !(duplicated(srtInput)|duplicated(srtInput,fromLast=TRUE))
      rm(srtInput)
    }
    input <- EL
    rm(EL)
    if (any(weights!=1)) input <- cbind(input,weights)
    cigraphObject$Arguments$edgelist <- TRUE
  }
  ### bnlearn ###
  if (is(input,"bn"))
  {
    # browser()
    bnobject <- input
    input <- as.matrix(bnobject$arcs)
    TempLabs  <- names(bnobject$nodes)
    if (is.null(cigraphObject$Arguments$labels)) cigraphObject$Arguments$labels  <- TempLabs
    
    input[] <- as.numeric(match(c(input), TempLabs))
    mode(input) <- "numeric"
    
    srtInput <- aaply(input,1,sort)
    input <- input[!duplicated(srtInput),]
    cigraphObject$Arguments$directed <- !(duplicated(srtInput)|duplicated(srtInput,fromLast=TRUE))
    cigraphObject$Arguments$directed <- cigraphObject$Arguments$directed[!duplicated(srtInput)]
    cigraphObject$Arguments$edgelist <- TRUE
  }
  if (is(input,"bn.strength"))
  {
    bnobject <- input
    input <- as.matrix(bnobject[c("from","to","strength")])
    TempLabs  <- unique(c(bnobject$from,bnobject$to))
    if (is.null(cigraphObject$Arguments$labels)) cigraphObject$Arguments$labels  <- TempLabs
    
    input[,1:2] <- as.numeric(match(c(input[,1:2]), TempLabs))
    
    input <- as.matrix(input)
    mode(input) <- "numeric"
    
    if (is.null(cigraphObject$Arguments$directed))
    {
      if (is.null(bnobject$direction) || all(bnobject$direction %in% c(0,0.5)))
      { 
        cigraphObject$Arguments$directed <- FALSE
      } else cigraphObject$Arguments$directed <- TRUE
    }
    
    if (!is.null(bnobject$direction))
    {
      input[,3] <- input[,3] * ( 1 - cigraphObject$Arguments$directed * (1- bnobject$direction ))
    }
    
    # remove undirect duplicates:
    srt <- cbind( pmin(input[,1],input[,2]), pmax(input[,1],input[,2]))
    input <- input[!(duplicated(srt)&!cigraphObject$Arguments$directed),  ]
    rm(srt)
    
    #     srtInput <- aaply(input,1,sort)
    #     input <- input[!duplicated(srtInput),]
    #     cigraphObject$Arguments$directed <- !(duplicated(srtInput)|duplicated(srtInput,fromLast=TRUE))
    #     cigraphObject$Arguments$directed <- cigraphObject$Arguments$directed[!duplicated(srtInput)]
    
    cigraphObject$Arguments$directed <- TRUE
    
    cigraphObject$Arguments$probabilityEdges <- TRUE
    
    if (is.null( cigraphObject$Arguments$parallelEdge))  cigraphObject$Arguments$parallelEdge <- TRUE
    
  }
  
  ### BDgraph ####
  if (is(input,"bdgraph"))
  {
    
    # browser()
    # stop("BDgraph support has temporarily been removed")
    
        if(is.null(cigraphObject$Arguments[['BDgraph']])){
          BDgraph=c("phat","Khat")
        } else {
          BDgraph=cigraphObject$Arguments[['BDgraph']]
        } 
        if (all(c("Khat","phat")%in%BDgraph)) layout(t(1:2))

        if(is.null(cigraphObject$Arguments[['BDtitles']])) BDtitles <- TRUE else BDtitles <- cigraphObject$Arguments[['BDtitles']]


        Res <- list()

        if (isTRUE(which(BDgraph == "phat") < which(BDgraph == "Khat")))
        {
          if(!requireNamespace("BDgraph")) stop("'BDgraph' package needs to be installed.")
          # phat:
          W <- as.matrix(BDgraph::plinks(input))
          W <- W + t(W)
          Res[["phat"]] <- do.call(cigraph,c(list(input=W,probabilityEdges = TRUE),cigraphObject$Arguments))
          L <- Res[["phat"]]$layout

          if (BDtitles) text(mean(par('usr')[1:2]),par("usr")[4] - (par("usr")[4] - par("usr")[3])/40,"Posterior probabilities", adj = c(0.5,1))

          # Khat:
          W <- as.matrix(input$K_hat)
          # diag(W) <- -1*diag(W)
          # W <-  - W / sqrt(diag(W)%o%diag(W))
          W <- wi2net(W)
          Res[["Khat"]] <- do.call(cigraph,c(list(input = W,layout = L), cigraphObject$Arguments))
          L <- Res[["Khat"]]$layout
          if (BDtitles) text(mean(par('usr')[1:2]),par("usr")[4] - (par("usr")[4] - par("usr")[3])/40,"Mean partial correlations", adj = c(0.5,1))

        } else
        {
          if ("Khat" %in% BDgraph)
          {
            W <- as.matrix(input$K_hat)
            # diag(W) <- -1*diag(W)
            # W <-  - W / sqrt(diag(W)%o%diag(W))
            W <- wi2net(input$K_hat)
            Res[["Khat"]] <- do.call(cigraph,c(list(input=W),cigraphObject$Arguments))
            L <- Res[["Khat"]]$layout
            if (BDtitles) text(mean(par('usr')[1:2]),par("usr")[4],"Mean partial correlations", adj = c(0.5,1))
          } else L <- cigraphObject$Arguments$layout

          if ("phat" %in% BDgraph)
          {
            W <- as.matrix(BDgraph::plinks(input))
            W <- W + t(W)
            Res[["phat"]] <- do.call(cigraph,c(list(input = W,layout = L,probabilityEdges= TRUE), cigraphObject$Arguments))
            if (BDtitles) text(mean(par('usr')[1:2]),par("usr")[4],"Posterior probabilities", adj = c(0.5,1))
          }
        }

        if (length(Res)==1) Res <- Res[[1]]
        return(Res)
    
  }
  
  
  ### GLASSO ###
  # glasso has no class but is a list with elements w, wi, loglik, errflag, approx, del and niter:
  if (is(input, "list") && all(c('w', 'wi', 'loglik','errflag', 'approx', 'del',  'niter' ) %in% names(input)))
  {
    input <- wi2net(input$wi)
  }
  
  
  ### Check arguments list:
  allArgs <- c("input", "layout", "groups", "minimum", "maximum", "cut", "details", 
               "threshold", "palette", "theme", "graph", "threshold", "sampleSize", 
               "tuning", "refit", "countDiagonal", "alpha", "bonf", "FDRcutoff", 
               "mar", "filetype", "filename", "width", "height", "normalize", "res",
               "DoNotPlot", "plot", "rescale", "standAlone", "color", "vsize", 
               "vsize2", "node.width", "node.height", "borders", "border.color", 
               "border.width", "shape", "polygonList", "vTrans", "subplots", 
               "subpars", "subplotbg", "images", "noPar", "pastel", "rainbowStart", 
               "usePCH", "node.resolution", "title", "preExpression", "postExpression", 
               "diag", "labels", "label.cex", "label.color", "label.prop", "label.norm", 
               "label.scale", "label.scale.equal", "label.font", "label.fill.vertical", 
               "label.fill.horizontal", "esize", "edge.width", "edge.color", 
               "posCol", "negCol", "unCol", "probCol", "negDashed", "probabilityEdges", 
               "colFactor", "trans", "fade", "loop", "lty", "edgeConnectPoints", 
               "curve", "curveAll", "curveDefault", "curveShape", "curveScale", 
               "curveScaleNodeCorrection", "curvePivot", "curvePivotShape", 
               "parallelEdge", "parallelAngle", "parallelAngleDefault", "edge.labels", 
               "edge.label.cex", "edge.label.bg", "edge.label.position", "edge.label.font", "edge.label.color",
               "repulsion", "layout.par", "layout.control", "aspect", "rotation", 
               "legend", "legend.cex", "legend.mode", "GLratio", "layoutScale", 
               "layoutOffset", "nodeNames", "bg", "bgcontrol", "bgres", "pty", 
               "gray",  "font", "directed", 
               "arrows", "arrowAngle", "asize", "open", "bidirectional", "mode", 
               "alpha", "sigScale", "bonf", "scores", "scores.range", "mode", 
               "edge.color", "knots", "knot.size", "knot.color", "knot.borders", 
               "knot.border.color", "knot.border.width", "means", "SDs", "meanRange", 
               "bars", "barSide", "barColor", "barLength", "barsAtSide", "pie", 
               "pieBorder", "pieColor", "pieColor2", "pieStart", "pieDarken", 
               "piePastel", "BDgraph", "BDtitles", "edgelist", "weighted", "nNodes", 
               "XKCD", "Edgelist", "Arguments", "plotOptions", "graphAttributes", 
               "layout", "layout.orig","resid","factorCors","residSize","filetype","model",
               "crossloadings","gamma","lambda.min.ratio","loopRotation","edgeConnectPoints","residuals","residScale","residEdge","CircleEdgeEnd","title.cex",  
               "node.label.offset", "node.label.position", "pieCImid", "pieCIlower", "pieCIupper", "pieCIpointcex", "pieCIpointcol",
               "edge.label.margin")
  
  if (any(!names(cigraphObject$Arguments) %in% allArgs)){
    wrongArgs <- names(cigraphObject$Arguments)[!names(cigraphObject$Arguments) %in% allArgs]
    warning(paste0("The following arguments are not documented and likely not arguments of cigraph and thus ignored: ",paste(wrongArgs,collapse = "; ")))
  }
  
  ## Extract arguments
  if(is.null(cigraphObject$Arguments[['verbose']]))
  {
    verbose <- FALSE
  } else verbose <- cigraphObject$Arguments[['verbose']] 
  
  if(is.null(cigraphObject$Arguments[['tuning']]))
  {
    tuning <- 0.5
  } else tuning <- cigraphObject$Arguments[['tuning']]  
  
  if(!is.null(cigraphObject$Arguments[['gamma']]))
  {
    tuning <- cigraphObject$Arguments[['gamma']]
  }
  
  if(is.null(cigraphObject$Arguments[['lambda.min.ratio']]))
  {
    lambda.min.ratio <- 0.01
  } else lambda.min.ratio <- cigraphObject$Arguments[['lambda.min.ratio']]  
 
  
  # Refit:
  if(is.null(cigraphObject$Arguments[['refit']]))
  {
    refit <- FALSE
  } else refit <- cigraphObject$Arguments[['refit']]  
  
  
  
  if(is.null(cigraphObject$Arguments[['FDRcutoff']]))
  {
    FDRcutoff <- 0.9
  } else FDRcutoff <- cigraphObject$Arguments[['FDRcutoff']]  
  
  
  ### HUGE (select via EBIC):
  if (is(input,"huge"))
  {
    if (input$method != "glasso") stop("Only 'glasso' method is supported")
    if(!requireNamespace("huge")) stop("'huge' package needs to be installed.")
    input <- huge::huge.select(input, "ebic", ebic.gamma = tuning)
  }
  
  ### HUGE select ###
  if (is(input,"select"))
  {
    if (input$method != "glasso") stop("Only 'glasso' method is supported")
    input <- wi2net(forceSymmetric(input$opt.icov))
  }
  
  # Coerce input to matrix:
  input <- as.matrix(input)
  
  # Set mode:
  sigSign <- FALSE
  if(is.null(cigraphObject$Arguments[['graph']])) graph <- "default" else graph=cigraphObject$Arguments[['graph']]
  if (graph == "fdr")
  {
    graph <- "fdr.cor"
  }
  if (graph == "EBICglasso"){
    graph <- "glasso"
  }
  if (graph == "ggmModSelect"){
    graph <- "ggmModSelect"
  }
  
  if (!graph %in% c("default","cor","pcor","glasso","ggmModSelect","factorial")){
    stop("'graph' argument must be one of 'default', 'cor', 'pcor', 'glasso', 'ggmModSelect', or 'factorial'")
  }
  
  # Reset graph for replotting:
  cigraphObject$Arguments[['graph']] <- NULL
  
  if (graph %in% c("sig2","significance2"))
  {
    graph <- "sig"
    sigSign <- TRUE
  }
  if (graph %in% c("sig","significance"))
  {
    #     if (!require("fdrtool")) stop("`fdrtool' package not found, is it installed?") 
    cigraphObject$Arguments[['mode']] <- "sig"
  }
  
  ### SIGNIFICANCE GRAPH ARGUMENTS ###
  if(is.null(cigraphObject$Arguments[['mode']])) mode <- "strength" else mode <- cigraphObject$Arguments[['mode']]
  if(is.null(cigraphObject$Arguments$sigScale)) sigScale <- function(x)0.7*(1-x)^(log(0.4/0.7,1-0.05)) else sigScale <- cigraphObject$Arguments$sigScale
  if (!mode%in%c("strength","sig","direct")) stop("Mode must be 'direct', 'sig' or 'strength'")	
  if(is.null(cigraphObject$Arguments$bonf)) bonf=FALSE else bonf=cigraphObject$Arguments$bonf
  if(is.null(cigraphObject$Arguments$OmitInsig)) OmitInsig=FALSE else OmitInsig <- cigraphObject$Arguments$OmitInsig
  if(is.null(cigraphObject$Arguments[['alpha']]))
  {
    if (mode != "sig")
    {
      alpha <- 0.05
    } else alpha <- c(0.0001,0.001,0.01,0.05) 
  } else alpha <- cigraphObject$Arguments[['alpha']]
  if (length(alpha) > 4) stop("`alpha' can not have length > 4")
  
  
  #####
  # Settings for the edgelist
  if(is.null(cigraphObject$Arguments$edgelist)) 
  {
    if (nrow(input)!=ncol(input)) {
      # Check if it is an edgelist or break:
      if (ncol(input) %in% c(2,3) && ((is.character(input[,1]) || is.factor(input[,1])) || all(input[,1] %% 1 == 0)) &&
          ((is.character(input[,2]) || is.factor(input[,2])) || all(input[,2] %% 1 == 0))){
        edgelist <- TRUE
      } else {
        stop("Input is not a weights matrix or an edgelist.")
      }
    } else edgelist <- FALSE
  } else edgelist=cigraphObject$Arguments$edgelist
  
  if(is.null(cigraphObject$Arguments[['edgeConnectPoints']])) edgeConnectPoints <- NULL else edgeConnectPoints <- cigraphObject$Arguments[['edgeConnectPoints']]
  
  
  if(is.null(cigraphObject$Arguments[['label.color.split']])) label.color.split <- 0.25 else label.color.split <- cigraphObject$Arguments[['label.color.split']]
  
  if(is.null(cigraphObject$Arguments$labels))
  {
    labels <- TRUE
    if (!edgelist && !is.null(colnames(input)))
    {
      #       if (nrow(input) <= 20 & all(colnames(input)==rownames(input)))
      #       {
      labels <- abbreviate(colnames(input),3)
      if (any(is.na(labels))){
        warning("Some labels where not abbreviatable.")
        labels <- ifelse(is.na(labels), colnames(input), labels)
      }
      #       }
    }
  } else labels <- cigraphObject$Arguments$labels
  
  
  if (edgelist)
  {
    if (is.character(input))
    {
      if(!is.logical(labels)) allNodes <- labels else allNodes <- unique(c(input[,1:2]))
      input[,1:2] <- match(input[,1:2],allNodes)
      input <- as.matrix(input)
      mode(input) <- "numeric"
      if (is.logical(labels) && labels) labels <- allNodes
    }
  }
  
  if(is.null(cigraphObject$Arguments$nNodes)) 
  {
    if (edgelist)
    {
      if (!is.logical(labels)) nNodes <- length(labels) else nNodes <- max(c(input[,1:2])) 
    } else nNodes=nrow(input)
  } else nNodes=cigraphObject$Arguments$nNodes
  
  
  #####
  
  
  #### Arguments for pies with Jonas
  # Arguments for pies:
  if(is.null(cigraphObject$Arguments[['pieRadius']])){
    pieRadius <- 1
  } else {
    pieRadius <- cigraphObject$Arguments[['pieRadius']]
  }
  
  if(is.null(cigraphObject$Arguments[['pieBorder']])){
    pieBorder <- .15
    if (any(pieBorder < 0 | pieBorder > 1)){
      stop("Values in the 'pieBorder' argument must be within [0,1]")
    }
  } else {
    pieBorder <- cigraphObject$Arguments[['pieBorder']]
  }
  
  if(is.null(cigraphObject$Arguments[['pieStart']])){
    pieStart <- 0
    if (any(pieStart < 0 | pieStart > 1)){
      stop("Values in the 'pieStart' argument must be within [0,1]")
    }
  } else {
    pieStart <- cigraphObject$Arguments[['pieStart']]
  }
  
  if(is.null(cigraphObject$Arguments[['pieDarken']])){
    pieDarken <- 0.25
    if (any(pieDarken < 0 | pieDarken > 1)){
      stop("Values in the 'pieDarken' argument must be within [0,1]")
    }
  } else {
    pieDarken <- cigraphObject$Arguments[['pieDarken']]
  }
  
  if(is.null(cigraphObject$Arguments[['pieColor']])){
    # pieColor <- 'grey'
    pieColor <- NA
  } else {
    pieColor <- cigraphObject$Arguments[['pieColor']]
  }
  
  
  if(is.null(cigraphObject$Arguments[['pieColor2']])){
    pieColor2 <- 'white'
  } else {
    pieColor2 <- cigraphObject$Arguments[['pieColor2']]
  }
  
  # Make arguments vectorized:
  if (length(pieColor) == 1){
    pieColor <- rep(pieColor,length=nNodes)
  }
  if (length(pieColor) != nNodes){
    stop("Length of 'pieColor' argument must be 1 or number of nodes")
  }
  
  if (length(pieColor2) == 1){
    pieColor2 <- rep(pieColor2,length=nNodes)
  }
  if (length(pieColor2) != nNodes){
    stop("Length of 'pieColor2' argument must be 1 or number of nodes")
  }
  
  if (length(pieBorder) == 1){
    pieBorder <- rep(pieBorder,length=nNodes)
  }
  if (length(pieBorder) != nNodes){
    stop("Length of 'pieBorder' argument must be 1 or number of nodes")
  }
  
  if (length(pieStart) == 1){
    pieStart <- rep(pieStart,length=nNodes)
  }
  if (length(pieStart) != nNodes){
    stop("Length of 'pieStart' argument must be 1 or number of nodes")
  }
  
  if (length(pieDarken) == 1){
    pieDarken <- rep(pieDarken,length=nNodes)
  }
  if (length(pieDarken) != nNodes){
    stop("Length of 'pieDarken' argument must be 1 or number of nodes")
  }
  
  
  
  
  
  if(is.null(cigraphObject$Arguments[['pie']])){
    drawPies <- FALSE
    pie <- NULL
  } else {
    # Obtain pie values:
    pie <- cigraphObject$Arguments[['pie']]
    
    # Check values:
    if (length(pie) != nNodes){
      stop("Length of 'pie' argument must be equal to number of nodes.")
    }
    #     if (any(pie < 0 | pie > 1)){
    #       stop("Values in the 'pie' argument must be within [0,1]")
    #     }
    
    
    # Dummy subplots (to be filed later)
    # subplots <- vector("list", nNodes)
    
    # Overwrite subplotbg to NA:
    # subplotbg <- NA
    
    # Overwrite borders to FALSE:
    # borders <- FALSE
    
    # Overwrite shape to circle:
    # shape <- "circle"
    
    # Logical:
    drawPies <- TRUE
  }
  
  # Pie CI:
  # Pie CI args:
  # "pieCIlower", "pieCIupper", "pieCIpointcex", "pieCIpointcol"
  
  pieCIs <- FALSE
  
  # Check if pieCIs are drawn:
  if(!is.null(cigraphObject$Arguments[['pieCIlower']])){
    pieCIs <- TRUE
    pieCIlower <- cigraphObject$Arguments[['pieCIlower']]
    if(is.null(cigraphObject$Arguments[['pieCIupper']])){
      pieCIupper <- 1
    }
  } 
  if(!is.null(cigraphObject$Arguments[['pieCIupper']])){
    pieCIs <- TRUE
    pieCIupper <- cigraphObject$Arguments[['pieCIupper']]
    if(is.null(cigraphObject$Arguments[['pieCIlower']])){
      pieCIlower <- 0
    }
  } 
  
  # Set up the pieCIs:
  if (isTRUE(pieCIs)){
    drawPies <- TRUE
    
    if (!is.null(cigraphObject$Arguments[['pieCImid']])){
      pieCImid <- cigraphObject$Arguments[['pieCImid']]
    } else stop("'pieCImid' may not be missing when pieCIs are used")
      
    
    # Vectorize:
    pieCIlower <- rep(pieCIlower,length=nNodes) 
    pieCIupper <- rep(pieCIupper,length=nNodes) 
    pieCImid <- rep(pieCImid, length=nNodes)
    
    # Check mid:
    if (any(pieCIlower > pieCImid | pieCIupper < pieCImid)){
      stop("'pieCImid' is not between 'pieCIlower' and 'pieCIupper'")
    }
    
    # Check bounds:
    if (any(pieCIlower < 0) | any(pieCImid > 1)){
      stop("pieCI range should be between 0 and 1")
    }
    
    # If pie argument is used, give an error:
    if(!is.null(cigraphObject$Arguments[['pie']])){
       stop("'pie' argument cannot be used in combination with pieCIs")
    }
    
    # get the size of the point:
    # "pieCIlower", "pieCIupper", "pieCIpointcex", "pieCIpointcol"
    
    if(!is.null(cigraphObject$Arguments[['pieCIpointcex']])){
      pieCIpointcex <- cigraphObject$Arguments[['pieCIpointcex']]
    } else {
      pieCIpointcex <- 0.01
    }
    pieCIpointcex <-  rep(pieCIpointcex, length=nNodes)
    
    if(!is.null(cigraphObject$Arguments[['pieCIpointcol']])){
      pieCIpointcol <- cigraphObject$Arguments[['pieCIpointcol']]
    } else {
      pieCIpointcol <- "black"
    }
    pieCIpointcol <- rep(pieCIpointcol, length = nNodes)
    
    # Now form the pie argument:
    pieTab <- cbind(
      0,
      pieCIlower,
      pmax(pieCIlower, pieCImid - (pieCIpointcex/2)),
      pmin(pieCIupper, pieCImid + (pieCIpointcex/2)),
      pieCIupper,
      1)
    pie <- list()
    pieColor2 <- list()
    
    for (i in seq_len(nrow(pieTab))){
      pie[[i]] <- diff(pieTab[i,])
      pieColor2[[i]] <- c("white",pieColor[i],pieCIpointcol[i],pieColor[i],"white")
    }
    pieColor <- pieColor2
  }
  
  
  
  #####
  
  
  if (is.expression(labels)) labels <- as.list(labels)
  
  if(is.null(cigraphObject$Arguments[['background']])) background <- NULL else background <- cigraphObject$Arguments[['background']]
  if(is.null(cigraphObject$Arguments[['label.prop']])){
    
    label.prop <- 0.9*(1-ifelse(pieBorder < 0.5,pieBorder,0))
    
  } else {
    label.prop <- cigraphObject$Arguments[['label.prop']]
  }
  
  if(is.null(cigraphObject$Arguments[['label.norm']])) label.norm <- "OOO" else label.norm <- cigraphObject$Arguments[['label.norm']]
  #   if(is.null(cigraphObject$Arguments[['label.cex']])) label.cex <- NULL else label.cex <- cigraphObject$Arguments[['label.cex']]
  #   
  
  if(is.null(cigraphObject$Arguments[['nodeNames']])) nodeNames <- NULL else nodeNames <- cigraphObject$Arguments[['nodeNames']]
  
  if(is.null(cigraphObject$Arguments[['subplots']])) {
    # if (!drawPies){
    subplots <- NULL       
    # }
  } else {
    #       if (drawPies){
    #         warning("'subplots' argument ignored if 'pie' argument is used.")     
    #       } else {
    subplots <- cigraphObject$Arguments[['subplots']]        
    # }
  }
  if(is.null(cigraphObject$Arguments[['subpars']])) subpars <- list(mar=c(0,0,0,0)) else subpars <- cigraphObject$Arguments[['subpars']]
  
  
  if(is.null(cigraphObject$Arguments[['subplotbg']])) {
    # if (!drawPies){
    subplotbg <- NULL       
    # }
  } else {
    #       if (drawPies){
    #         warning("'subplotbg' argument ignored if 'pie' argument is used.")
    #       } else {
    subplotbg <- cigraphObject$Arguments[['subplotbg']]        
    # }
  }
  
  if(is.null(cigraphObject$Arguments[['images']])) images <- NULL else images <- cigraphObject$Arguments[['images']]
  
  if(is.null(cigraphObject$Arguments[['noPar']])) noPar <- FALSE else noPar <- cigraphObject$Arguments[['noPar']]
  
  
  
  # Knots:
  if(is.null(cigraphObject$Arguments[['knots']])) knots <- list() else knots <- cigraphObject$Arguments[['knots']]
  if(is.null(cigraphObject$Arguments[['knot.size']])) knot.size <- 1 else knot.size <- cigraphObject$Arguments[['knot.size']]
  if(is.null(cigraphObject$Arguments[['knot.color']])) knot.color <- NA else knot.color <- cigraphObject$Arguments[['knot.color']]
  if(is.null(cigraphObject$Arguments[['knot.borders']])) knot.borders <- FALSE else knot.borders <- cigraphObject$Arguments[['knot.borders']]
  if(is.null(cigraphObject$Arguments[['knot.border.color']])) knot.border.color <- "black" else knot.border.color <- cigraphObject$Arguments[['knot.border.color']]
  if(is.null(cigraphObject$Arguments[['knot.border.width']])) knot.border.width <- 1 else knot.border.width <- cigraphObject$Arguments[['knot.border.width']]
  
  
  #####
  
  
  
  
  if(is.null(cigraphObject$Arguments$shape))  {
    # if (!drawPies){
    shape <- rep("circle",nNodes) 
    if (!is.null(subplots))
    {
      # Get which nodes become a subplot:
      whichsub <- which(sapply(subplots,function(x)is.expression(x)|is.function(x)))
      
      shape[whichsub][!shape[whichsub]%in%c("square","rectangle")] <- "square"
    }      
    # }
  } else {
    #     if (drawPies){
    #       warning("'shape' argument ignored if 'pie' argument is used.")
    #     } else {
    shape <- cigraphObject$Arguments[['shape']]        
    # }
  }
  
  
  if(is.null(cigraphObject$Arguments[['usePCH']])) 
  {
    if (nNodes > 50 && !drawPies) usePCH <- TRUE else usePCH <- NULL 
  } else usePCH <- cigraphObject$Arguments[['usePCH']]
  
  if(is.null(cigraphObject$Arguments[['node.resolution']])) node.resolution <- 100 else node.resolution <- cigraphObject$Arguments[['node.resolution']]
  
  
  # Default for fact cut and groups
  if (graph=="factorial") fact=TRUE else fact=FALSE
  if (fact & edgelist) stop('Factorial graph needs a correlation matrix')
  #     if (graph=="concentration") partial=TRUE else partial=FALSE
  
  #   if(is.null(cigraphObject$Arguments$cutQuantile)) cutQuantile <- 0.9 else cutQuantile <- cigraphObject$Arguments$cutQuantile
  defineCut <- FALSE
  if(is.null(cigraphObject$Arguments[['cut']])) 
  {
    cut=0
    #       if (nNodes<50) 
    if (nNodes>=20 | fact) 
    {
      cut=0.3
      defineCut <- TRUE
    }
    if (mode=="sig") cut <- ifelse(length(alpha)>1,sigScale(alpha[length(alpha)-1]),sigScale(alpha[length(alpha)]))
  } else if (mode != "sig") cut <- ifelse(is.na(cigraphObject$Arguments[['cut']]),0,cigraphObject$Arguments[['cut']]) else cut <- ifelse(length(alpha)>1,sigScale(alpha[length(alpha)-1]),sigScale(alpha[length(alpha)]))
  
  if(is.null(cigraphObject$Arguments$groups)) groups=NULL else groups=cigraphObject$Arguments$groups
  
  if (is.factor(groups) | is.character(groups)) groups <- tapply(1:length(groups),groups,function(x)x)
  
  
  
  
  
  # Factorial graph:
  if(is.null(cigraphObject$Arguments$nfact))
  {
    nfact=NULL
  } else nfact=cigraphObject$Arguments$nfact
  
  if (fact)
  {
    if (is.null(nfact)) 
    {
      if (is.null(groups)) nfact=sum(eigen(input)$values>1) else nfact=length(groups)
    }
    
    loadings=loadings(factanal(factors=nfact,covmat=input,rotation="promax"))
    
    loadings=loadings[1:nrow(loadings),1:ncol(loadings)]
    
    loadings[loadings<cut]=0
    loadings[loadings>=cut]=1
    
    input=(loadings%*%t(loadings)>0)*1
    
    diag(input)=0
  }
  
  # Glasso arguments:
  if(is.null(cigraphObject$Arguments[['sampleSize']]))
  {
    sampleSize <- NULL
  } else sampleSize <- cigraphObject$Arguments[['sampleSize']]
  
  if(is.null(cigraphObject$Arguments[['countDiagonal']]))
  {
    countDiagonal <- FALSE
  } else countDiagonal <- cigraphObject$Arguments[['countDiagonal']]
  
  
  
  # SET DEFAULT cigraphObject$Arguments:
  # General cigraphObject$Arguments:
  if(is.null(cigraphObject$Arguments$DoNotPlot)) DoNotPlot=FALSE else DoNotPlot=cigraphObject$Arguments$DoNotPlot
  if(is.null(cigraphObject$Arguments[['layout']])) layout=NULL else layout=cigraphObject$Arguments[['layout']]
  if(is.null(cigraphObject$Arguments$maximum)) maximum=0 else maximum=cigraphObject$Arguments$maximum
  if(is.null(cigraphObject$Arguments$minimum))
  {
    #     if (nNodes<50)  minimum=0
    #     if (nNodes>=50)  minimum=0.1
    minimum <- 0
    if (mode=="sig") minimum <- ifelse(length(alpha)>1,sigScale(alpha[length(alpha)]),0)
  } else 
  {
    if (mode!="sig") minimum=cigraphObject$Arguments$minimum else minimum <- ifelse(length(alpha)>1,sigScale(alpha[length(alpha)]),0)
    if (is.character(minimum))
    {
      if (grepl("sig",minimum,ignore.case = TRUE))
      {
        if (is.null(sampleSize))
        {
          stop("'sampleSize' argument must be assigned to use significance as minimum")
        }
        if (graph == "default")
        {
          warning("'graph' argument did not specify type of graph. Assuming correlation graph (graph = 'cor')")
          graph <- "cor"
        }
        if (graph %in% c("cor","pcor")) {
          # Find threshold for significance!
          # difference between cor and pcor is in df:
          if (graph == "cor")
          {
            df <- sampleSize - 2
          } else {
            df <- sampleSize - 2 - (nNodes - 2)
          }
          siglevel <- max(alpha)/2
          if (bonf)
          {
            siglevel <- siglevel / (nNodes*(nNodes-1)/2)
          }
          t <- abs(qt(siglevel, df, lower.tail=TRUE))
          minimum <- t/sqrt(t^2+df) 
        } else stop("minimum = 'sig' is not supported with this 'graph' argument")
        
      } else stop("Minimum is specified a string which is not 'sig'.")
    }
  }
  if (minimum < 0)
  {
    warning("'minimum' set to absolute value")
    minimum <- abs(minimum)
  }
  
  # Threshold argument removes edges from network:
  if(is.null(cigraphObject$Arguments[['threshold']]))
  {
    threshold <- 0
  } else {
    threshold <- cigraphObject$Arguments[['threshold']]
  }
  
  if(is.null(cigraphObject$Arguments$weighted)) weighted=NULL else weighted=cigraphObject$Arguments$weighted
  if(is.null(cigraphObject$Arguments$rescale)) rescale=TRUE else rescale=cigraphObject$Arguments$rescale
  if(is.null(cigraphObject$Arguments[['edge.labels']])) edge.labels=FALSE else edge.labels=cigraphObject$Arguments[['edge.labels']]
  if(is.null(cigraphObject$Arguments[['edge.label.bg']])) edge.label.bg=TRUE else edge.label.bg=cigraphObject$Arguments[['edge.label.bg']]
  
  if (identical(FALSE,edge.label.bg)) plotELBG <- FALSE else plotELBG <- TRUE

  if(is.null(cigraphObject$Arguments[['edge.label.margin']])) edge.label.margin=0 else edge.label.margin=cigraphObject$Arguments[['edge.label.margin']]
  
  
    
  ### Themes ###
  # Default theme:
  posCol <- c("#009900","darkgreen")
  negCol <- c("#BF0000","red")
  bcolor <- NULL
  bg <- FALSE
  negDashed <- FALSE
  parallelEdge <- FALSE
  fade <- NA 
  border.width <- 1 
  font <- 1 
  unCol <- "#808080" 
  # if (length(groups) < 8){
  #   palette <- "colorblind"
  # } else {
    palette <- "rainbow"
  # }
  
  if(!is.null(cigraphObject$Arguments[['theme']])){
    theme <- cigraphObject$Arguments[['theme']]
    if (length(theme) > 1) stop("'theme' must be of lenght 1")
    if (!theme %in% c("classic","Hollywood","Leuven","Reddit","TeamFortress","Fried",
                      "Borkulo","colorblind","gray","gimme","GIMME","neon","pride")){
      stop(paste0("Theme '",theme,"' is not supported."))
    }
 
    # Themes:
    if (theme == "classic"){
      posCol <- c("#009900","darkgreen")
      negCol <- c("#BF0000","red")
    } else if (theme == "Leuven"){
      dots <- list(...)
      dots$DoNotPlot <- TRUE
      dots$theme <- "classic"
      dots$input <- input
      return(getWmat(do.call(cigraph,dots )))
    } else if (theme == "Hollywood"){
      negCol <- "#FFA500"
      posCol <- "#005AFF"
    } else if (theme == "Reddit"){
      posCol <- "#CCCCFF"
      negCol <- "#FF4500"
    } else if (theme == "TeamFortress"){
      negCol <- "#B8383B"
      posCol <- "#5885A2"
    } else if (theme == "Fried"){
      posCol <- "black"
      negCol <- "black"
      bg <- "gray"
      palette <- "gray"
    } else if (theme == "Borkulo"){
      posCol <- "darkblue"
      negCol <- "red"
      bcolor <- "darkblue"
    } else if (theme == "colorblind"){
      posCol <- c("#0000D5","darkblue")
      negCol <- c("#BF0000","red")
      palette <- "colorblind"
    } else if (theme == "gray" | theme == "grey"){
      posCol <- negCol <- c("gray10","black")
      palette <- "gray"
      negDashed <- TRUE
    } else if (theme == "gimme" | theme == "GIMME"){
      posCol <- "red"
      negCol <- "blue"
      parallelEdge <- TRUE
      fade <- FALSE
    } else if(theme == "neon"){
      bg <- "black"
      label.color <- "#8ffcff"
      bcolor <- "#8ffcff"
      border.width <- 4
      font <- 2
      posCol <- "#f3ea5f"
      negCol <- "#c04df9"
      unCol <- "#8ffcff"
      palette <- "neon"
    }else if(theme == "pride"){
      posCol <- "#1AB3FF"
      negCol <- "#FF1C8D"
      unCol <- "#613915"
      palette <- "pride"
    }
  }
  
  # Overwrite:
    if(!is.null(cigraphObject$Arguments[['parallelEdge']]))  parallelEdge <- cigraphObject$Arguments[['parallelEdge']]
    
    if(!is.null(cigraphObject$Arguments[['fade']])) fade <- cigraphObject$Arguments[['fade']]
    
  if(!is.null(cigraphObject$Arguments[['negDashed']])) negDashed <- cigraphObject$Arguments[['negDashed']]
  if(!is.null(cigraphObject$Arguments[['posCol']])) posCol <- cigraphObject$Arguments[['posCol']]
  if(!is.null(cigraphObject$Arguments[['negCol']])) negCol <- cigraphObject$Arguments[['negCol']]
  if(!is.null(cigraphObject$Arguments[['border.width']])) border.width <- cigraphObject$Arguments[['border.width']]
    
  if(!is.null(cigraphObject$Arguments[['font']])) font <- cigraphObject$Arguments[['font']]
  if(is.null(cigraphObject$Arguments[['edge.label.font']])) edge.label.font=font else edge.label.font=cigraphObject$Arguments[['edge.label.font']]
  if(is.null(cigraphObject$Arguments[['label.font']])) label.font <- font else label.font <- cigraphObject$Arguments[['label.font']]
   if(!is.null(cigraphObject$Arguments[['unCol']])) unCol <- cigraphObject$Arguments[['unCol']] 
    
  
  if(is.null(cigraphObject$Arguments[['probCol']])) probCol <- "black" else probCol <- cigraphObject$Arguments[['probCol']]
  if(!is.null(cigraphObject$Arguments[['probabilityEdges']])) 
  {
    if (isTRUE(cigraphObject$Arguments[['probabilityEdges']]))
    {
      posCol <- probCol
    }
  }
    
  if (length(posCol)==1) posCol <- rep(posCol,2)
  if (length(posCol)!=2) stop("'posCol' must be of length 1 or 2.")
  if (length(negCol)==1) negCol <- rep(negCol,2)
  if (length(negCol)!=2) stop("'negCol' must be of length 1 or 2.")
  
  # border color:
  if(!is.null(cigraphObject$Arguments[['border.color']])) {
    bcolor <- cigraphObject$Arguments[['border.color']]
  }
  # Alias?
  if(!is.null(cigraphObject$Arguments[['border.colors']])) {
    bcolor <- cigraphObject$Arguments[['border.colors']]
  }
  
  # BG:
  if(!is.null(cigraphObject$Arguments$bg)) bg <- cigraphObject$Arguments$bg
  
  # Palette:
  
  # PALETTE either one of the defaults or a function
  if(!is.null(cigraphObject$Arguments[['palette']])){
    palette <- cigraphObject$Arguments[['palette']]
  }
  
  # Check palette:
  if (!is.function(palette)){
    if (length(palette) != 1 && !is.character(palette)){
      stop("'palette' must be a single string.")
    }
    if (!palette %in% c("rainbow","colorblind","R","ggplot2","gray","grey","pastel","neon","pride")){
      stop(paste0("Palette '",palette,"' is not supported."))
    }
  }
  
  
  ###
  

  if(is.null(cigraphObject$Arguments[['colFactor']])) colFactor <- 1 else colFactor <- cigraphObject$Arguments[['colFactor']]
  
  if(is.null(cigraphObject$Arguments[['edge.color']])) edge.color <- NULL else edge.color=cigraphObject$Arguments[['edge.color']]
  if(is.null(cigraphObject$Arguments[['edge.label.cex']])) edge.label.cex=1 else edge.label.cex=cigraphObject$Arguments[['edge.label.cex']]
  if(is.null(cigraphObject$Arguments[['edge.label.position']])) edge.label.position <- 0.5 else edge.label.position=cigraphObject$Arguments[['edge.label.position']]
  
  
  if(is.null(cigraphObject$Arguments$directed))
  {
    if (edgelist) directed=TRUE else directed=NULL 
  } else directed=cigraphObject$Arguments$directed
  if(is.null(cigraphObject$Arguments[['legend']]))
  {
    if ((!is.null(groups) & !is.null(names(groups))) | !is.null(nodeNames)) legend <- TRUE else legend <- FALSE
  } else legend <- cigraphObject$Arguments[['legend']]
  
  stopifnot(is.logical(legend))
  
  #     if (is.null(groups)) legend <- FALSE
  if(is.null(cigraphObject$Arguments$plot)) plot=TRUE else plot=cigraphObject$Arguments$plot
  if(is.null(cigraphObject$Arguments$rotation)) rotation=NULL else rotation=cigraphObject$Arguments$rotation
  if(is.null(cigraphObject$Arguments[['layout.control']])) layout.control=0.5 else layout.control=cigraphObject$Arguments[['layout.control']]
  
  # repulsion controls the repulse.rad argument
  if(is.null(cigraphObject$Arguments[['repulsion']])) repulsion=1 else repulsion=cigraphObject$Arguments[['repulsion']]
  if(is.null(cigraphObject$Arguments[['layout.par']])) {
    if (is.null(layout) || identical(layout,"spring")) layout.par <- list(repulse.rad = nNodes^(repulsion * 3))  else layout.par <- list()
  } else layout.par=cigraphObject$Arguments[['layout.par']]
  
  if(is.null(cigraphObject$Arguments[['layoutRound']])){
    layoutRound <- TRUE
  } else { 
    layoutRound <- cigraphObject$Arguments[['layoutRound']]
  }
  layout.par$round <- layoutRound
  
  if(is.null(cigraphObject$Arguments$details)) details=FALSE else details=cigraphObject$Arguments$details
  if(is.null(cigraphObject$Arguments$title)) title <- NULL else title <- cigraphObject$Arguments$title
  
  if(is.null(cigraphObject$Arguments[['title.cex']])) title.cex <- NULL else title.cex <- cigraphObject$Arguments[['title.cex']]
  
  if(is.null(cigraphObject$Arguments$preExpression)) preExpression <- NULL else preExpression <- cigraphObject$Arguments$preExpression
  if(is.null(cigraphObject$Arguments$postExpression)) postExpression <- NULL else postExpression <- cigraphObject$Arguments$postExpression
  
  
  # Output cigraphObject$Arguments:
  
  if(is.null(cigraphObject$Arguments[['edge.label.color']])) ELcolor <- NULL else ELcolor <- cigraphObject$Arguments[['edge.label.color']]
  
  
  
  # if(is.null(cigraphObject$Arguments[['border.width']])) border.width <- 1 else border.width <- cigraphObject$Arguments[['border.width']]
  #if (!DoNotPlot & !is.null(dev.list()[dev.cur()]))
  #{
  #	par(mar=c(0,0,0,0), bg=background)
  #	if (plot)
  #	{
  #		plot(1, ann = FALSE, axes = FALSE, xlim = c(-1.2, 1.2), ylim = c(-1.2 ,1.2),type = "n", xaxs = "i", yaxs = "i")
  #		plot <- FALSE
  #	}
  #}
  
  PlotOpen <- !is.null(dev.list()[dev.cur()])
  
  if(is.null(cigraphObject$Arguments$filetype)) filetype="default" else filetype=cigraphObject$Arguments$filetype
  if(is.null(cigraphObject$Arguments$filename)) filename="cigraph" else filename=cigraphObject$Arguments$filename
  if(is.null(cigraphObject$Arguments$width)) width <- 7 else width <- cigraphObject$Arguments[['width']]
  if(is.null(cigraphObject$Arguments$height)) height <- 7 else height <- cigraphObject$Arguments[['height']]
  if(is.null(cigraphObject$Arguments$pty)) pty='m' else pty=cigraphObject$Arguments$pty
  if(is.null(cigraphObject$Arguments$res)) res=320 else res=cigraphObject$Arguments$res
  if(is.null(cigraphObject$Arguments[['normalize']])) normalize <- TRUE else normalize <- cigraphObject$Arguments[['normalize']]
  
  # Graphical cigraphObject$Arguments
  #     defNodeSize <- max((-1/72)*(nNodes)+5.35,1) ### Default node size, used as standard unit.
  if(is.null(cigraphObject$Arguments[['mar']])) mar <- c(3,3,3,3)/10 else mar <- cigraphObject$Arguments[["mar"]]/10
  if(is.null(cigraphObject$Arguments[['vsize']])) 
  {
    vsize <- 8*exp(-nNodes/80)+1
    #     vsize <- max((-1/72)*(nNodes)+5.35,1)
    if(is.null(cigraphObject$Arguments[['vsize2']])) vsize2 <- vsize else vsize2 <- vsize * cigraphObject$Arguments[['vsize2']]
  } else {
    vsize <- cigraphObject$Arguments[['vsize']]
    if(is.null(cigraphObject$Arguments[['vsize2']])) vsize2 <- vsize else vsize2 <- cigraphObject$Arguments[['vsize2']]
  }
  
  if(!is.null(cigraphObject$Arguments[['node.width']])) 
  {
    vsize <- vsize * cigraphObject$Arguments[['node.width']]
  }
  
  if(!is.null(cigraphObject$Arguments[['node.height']])) 
  {
    vsize2 <- vsize2 * cigraphObject$Arguments[['node.height']]
  }
  
  if(is.null(cigraphObject$Arguments$color)) color=NULL else color=cigraphObject$Arguments$color
  
  if(is.null(cigraphObject$Arguments[['gray']])) gray <- FALSE else gray <- cigraphObject$Arguments[['gray']]
  
  if (gray) {
    posCol <- negCol <- c("gray10","black")
    warning("The 'gray' argument is deprecated, please use theme = 'gray' instead.")
  }
  
  if(is.null(cigraphObject$Arguments[['pastel']])){
    pastel <- FALSE 
  } else {
    warning("The 'pastel' argument is deprecated, please use palette = 'pastel' instead.")
    palette <- "pastel"
    pastel <- cigraphObject$Arguments[['pastel']]
  }
  
  
  
  if(is.null(cigraphObject$Arguments[['piePastel']])) piePastel <- FALSE else piePastel <- cigraphObject$Arguments[['piePastel']]
  if(is.null(cigraphObject$Arguments[['rainbowStart']])) rainbowStart <- 0 else rainbowStart <- cigraphObject$Arguments[['rainbowStart']]
  
  if(is.null(cigraphObject$Arguments$bgcontrol)) bgcontrol=6 else bgcontrol=cigraphObject$Arguments$bgcontrol
  if(is.null(cigraphObject$Arguments$bgres)) bgres=100 else bgres=cigraphObject$Arguments$bgres
  if(is.null(cigraphObject$Arguments[['trans',exact=FALSE]])) transparency <- NULL else transparency <- cigraphObject$Arguments[['trans',exact=FALSE]]
  if (is.null(transparency))
  {
    if (isTRUE(bg)) transparency <- TRUE else transparency <- FALSE
  }

  
  # Automatic fading?
  # autoFade <- isTRUE(fade)
  # if (isTRUE(fade)){
  #   fade <- NA
  # }
  # 
  if (is.logical(fade)){
    fade <- ifelse(fade,NA,1)
  }
  
  if (identical(fade,FALSE)){
    fade <- 1
  }
  
  if(is.null(cigraphObject$Arguments[['loop']])) loop=1 else loop=cigraphObject$Arguments[['loop']]
  if(is.null(cigraphObject$Arguments[['loopRotation']]))
  {
    loopRotation <- NA
  } else {
    loopRotation=cigraphObject$Arguments[['loopRotation']]
  }
  
  if(is.null(cigraphObject$Arguments[['residuals']])) residuals=FALSE else residuals=cigraphObject$Arguments[['residuals']]
  if(is.null(cigraphObject$Arguments[['residScale']])) residScale=1 else residScale=cigraphObject$Arguments[['residScale']]
  if(is.null(cigraphObject$Arguments[['residEdge']])) residEdge=FALSE else residEdge=cigraphObject$Arguments[['residEdge']]
  
  if(is.null(cigraphObject$Arguments[['bars']])) bars <- list() else bars <- cigraphObject$Arguments[['bars']]
  if(is.null(cigraphObject$Arguments[['barSide']])) barSide <- 1 else barSide <- cigraphObject$Arguments[['barSide']]
  if(is.null(cigraphObject$Arguments[['barLength']])) barLength <- 0.5 else barLength <- cigraphObject$Arguments[['barLength']]
  if(is.null(cigraphObject$Arguments[['barColor']])) barColor <- 'border' else barColor <- cigraphObject$Arguments[['barColor']]
  if(is.null(cigraphObject$Arguments[['barsAtSide']])) barsAtSide <- FALSE else barsAtSide <- cigraphObject$Arguments[['barsAtSide']]
  
  # Means and SDs:
  if(is.null(cigraphObject$Arguments[['means']])) means <- NA else means <- cigraphObject$Arguments[['means']]
  if(is.null(cigraphObject$Arguments[['SDs']])) SDs <- NA else SDs <- cigraphObject$Arguments[['SDs']]
  if(is.null(cigraphObject$Arguments[['meanRange']])) {
    if (all(is.na(means))) meanRange <- c(NA,NA) else meanRange <- range(means,na.rm=TRUE) 
  }else meanRange <- cigraphObject$Arguments[['meanRange']]
  
  
  if (!is.list(bars)) bars <- as.list(bars)
  
  if(is.null(cigraphObject$Arguments[['CircleEdgeEnd']])) CircleEdgeEnd=FALSE else CircleEdgeEnd=cigraphObject$Arguments[['CircleEdgeEnd']]
  if(is.null(cigraphObject$Arguments[['loopAngle']])) loopangle=pi/2 else loopAngle=cigraphObject$Arguments[['loopAngle']]
  if(is.null(cigraphObject$Arguments[['legend.cex']])) legend.cex=0.6 else legend.cex=cigraphObject$Arguments[['legend.cex']]
  if(is.null(cigraphObject$Arguments[['legend.mode']]))
  {
    if (!is.null(nodeNames) && !is.null(groups)){
      legend.mode <- "style1" # or style2
    } else if (!is.null(nodeNames)) legend.mode <- "names" else legend.mode <- "groups"
  }  else legend.mode=cigraphObject$Arguments[['legend.mode']]
  
  if(is.null(cigraphObject$Arguments$borders)){
    # if (!drawPies){
    borders <- TRUE       
    # }
  } else {
    #     if (drawPies){
    #       warning("'borders' argument ignored if 'pie' argument is used.")     
    #     } else {
    borders <- cigraphObject$Arguments[['borders']]        
    # }
  }
  
  
  ### Polygon lookup list:
  polygonList = list(
    ellipse = ELLIPSEPOLY,
    heart  = HEARTPOLY,
    star = STARPOLY,
    crown = CROWNPOLY
  )
  
  if(!is.null(cigraphObject$Arguments[['polygonList']])) polygonList  <- c( polygonList, cigraphObject$Arguments[['polygonList']])
  
  # Rescale to -1 - 1 and compute radians per point:
  for (i in seq_along(polygonList))
  {
    polygonList[[i]]$x <- (polygonList[[i]]$x - min(polygonList[[i]]$x)) / (max(polygonList[[i]]$x) - min(polygonList[[i]]$x)) * 2 - 1
    polygonList[[i]]$y <- (polygonList[[i]]$y - min(polygonList[[i]]$y)) / (max(polygonList[[i]]$y) - min(polygonList[[i]]$y)) * 2 - 1
  }
  
  if(is.null(cigraphObject$Arguments[['label.scale']])) label.scale=TRUE else label.scale=cigraphObject$Arguments[['label.scale']]
  
  if(is.null(cigraphObject$Arguments[['label.cex']])){ 
    if (label.scale){
      label.cex <- 1  
    } else {
      label.cex <- 1
    }
  } else label.cex <- cigraphObject$Arguments[['label.cex']]
  
  if(is.null(cigraphObject$Arguments$label.scale.equal)) label.scale.equal=FALSE else label.scale.equal=cigraphObject$Arguments$label.scale.equal
  
  if(is.null(cigraphObject$Arguments$label.fill.horizontal)) label.fill.horizontal<-1 else label.fill.horizontal <- cigraphObject$Arguments$label.fill.horizontal
  if(is.null(cigraphObject$Arguments$label.fill.vertical)) label.fill.vertical<-1 else label.fill.vertical <- cigraphObject$Arguments$label.fill.vertical
  if(is.null(cigraphObject$Arguments$node.label.offset)) node.label.offset<-c(0.5, 0.5) else node.label.offset <- cigraphObject$Arguments$node.label.offset
  if(is.null(cigraphObject$Arguments$node.label.position)) node.label.position<-NULL else node.label.position <- cigraphObject$Arguments$node.label.position
  
  
  
  if(is.null(cigraphObject$Arguments$scores)) scores=NULL else scores=cigraphObject$Arguments$scores
  if(is.null(cigraphObject$Arguments$scores.range)) scores.range=NULL else scores.range=cigraphObject$Arguments$scores.range
  if(is.null(cigraphObject$Arguments$lty)) lty=1 else lty=cigraphObject$Arguments$lty
  if(is.null(cigraphObject$Arguments$vTrans)) vTrans=255 else vTrans=cigraphObject$Arguments$vTrans
  # if(is.null(cigraphObject$Arguments[['overlay']])) overlay <- FALSE else overlay <- cigraphObject$Arguments[['overlay']]
  # if(is.null(cigraphObject$Arguments[['overlaySize']])) overlaySize <- 0.5 else overlaySize <- cigraphObject$Arguments[['overlaySize']]
  if(is.null(cigraphObject$Arguments[['GLratio']])) GLratio <- 2.5 else GLratio <- cigraphObject$Arguments[['GLratio']]
  if(is.null(cigraphObject$Arguments$layoutScale)) layoutScale <- 1 else layoutScale <- cigraphObject$Arguments$layoutScale
  if(is.null(cigraphObject$Arguments[['layoutOffset']])) layoutOffset <- 0 else layoutOffset <- cigraphObject$Arguments[['layoutOffset']]
  
  # Aspect ratio:
  if(is.null(cigraphObject$Arguments[['aspect']])) aspect=FALSE else aspect=cigraphObject$Arguments[['aspect']]
  
  # cigraphObject$Arguments for directed graphs:
  if(is.null(cigraphObject$Arguments[['curvePivot']])) curvePivot <- FALSE else curvePivot <- cigraphObject$Arguments[['curvePivot']]
  if (isTRUE(curvePivot)) curvePivot <- 0.1
  if(is.null(cigraphObject$Arguments[['curveShape']])) curveShape <- -1 else curveShape <- cigraphObject$Arguments[['curveShape']]
  if(is.null(cigraphObject$Arguments[['curvePivotShape']])) curvePivotShape <- 0.25 else curvePivotShape <- cigraphObject$Arguments[['curvePivotShape']]
  if(is.null(cigraphObject$Arguments[['curveScale']])) curveScale <- TRUE else curveScale <- cigraphObject$Arguments[['curveScale']]
  
  if(is.null(cigraphObject$Arguments[['curveScaleNodeCorrection']])) curveScaleNodeCorrection <- TRUE else curveScaleNodeCorrection <- cigraphObject$Arguments[['curveScaleNodeCorrection']]
  
  
  
  
  

  
  if(is.null(cigraphObject$Arguments[['parallelAngle']])) parallelAngle <- NA else parallelAngle <- cigraphObject$Arguments[['parallelAngle']]
  
  if(is.null(cigraphObject$Arguments[['parallelAngleDefault']])) parallelAngleDefault <- pi/6 else parallelAngleDefault <- cigraphObject$Arguments[['parallelAngleDefault']]
  
  if(is.null(cigraphObject$Arguments[['curveDefault']])) curveDefault <- 1 else curveDefault <- cigraphObject$Arguments[['curveDefault']]
  
  if(is.null(cigraphObject$Arguments[['curve']]))
  {
    if (any(parallelEdge))
    { 
      curve <- ifelse(parallelEdge,0,NA)
    } else curve <- NA 
  } else {      
    curve <- cigraphObject$Arguments[['curve']]
    if (length(curve)==1) 
    {
      curveDefault <- curve
      curve <- NA
    }
  }
  if(is.null(cigraphObject$Arguments[['curveAll']])) curveAll <- FALSE else curveAll <- cigraphObject$Arguments[['curveAll']]
  if (curveAll)
  {
    curve[is.na(curve)] <- curveDefault
  }
  if(is.null(cigraphObject$Arguments$arrows)) arrows=TRUE else arrows=cigraphObject$Arguments$arrows
  #     asize=asize*2.4/height
  if(is.null(cigraphObject$Arguments$open)) open=FALSE else open=cigraphObject$Arguments$open
  if(is.null(cigraphObject$Arguments$bidirectional)) bidirectional=FALSE else bidirectional=cigraphObject$Arguments$bidirectional
  
  # cigraphObject$Arguments for SVG pictures:
  # if(is.null(cigraphObject$Arguments$tooltips)) tooltips=NULL else tooltips=cigraphObject$Arguments$tooltips
  # if(is.null(cigraphObject$Arguments$SVGtooltips)) SVGtooltips=NULL else SVGtooltips=cigraphObject$Arguments$SVGtooltips
  if(is.null(cigraphObject$Arguments$hyperlinks)) hyperlinks=NULL else hyperlinks=cigraphObject$Arguments$hyperlinks
  
  # cigraphObject$Arguments for TEX:
  if(is.null(cigraphObject$Arguments$standAlone)) standAlone=TRUE else standAlone=cigraphObject$Arguments$standAlone
  
  ### EASTER EGGS ###
  if(is.null(cigraphObject$Arguments[['XKCD']])) XKCD <- FALSE else XKCD <- TRUE
  
  #     # Legend setting 1
  #     if (is.null(legend))
  #     {
  #       if (is.null(groups)) legend=FALSE else legend=TRUE
  #     }
  #if ((legend & filetype!='pdf' & filetype!='eps') | filetype=="svg")
  if ((legend&is.null(scores))|(identical(filetype,"svg")))
  {
    width=width*(1+(1/GLratio))
  }
  
  #     if (!DoNotPlot)
  #     {
  #       
  #       # Start output:
  #       if (filetype=='default') if (is.null(dev.list()[dev.cur()])) dev.new(rescale="fixed",width=width,height=height)
  #       if (filetype=='R') dev.new(rescale="fixed",width=width,height=height)
  #       if (filetype=='X11' | filetype=='x11') x11(width=width,height=height)
  #       if (filetype=='eps') postscript(paste(filename,".eps",sep=""),height=height,width=width, horizontal=FALSE)
  #       if (filetype=='pdf') pdf(paste(filename,".pdf",sep=""),height=height,width=width)
  #       if (filetype=='tiff') tiff(paste(filename,".tiff",sep=""),units='in',res=res,height=height,width=width)
  #       if (filetype=='png') png(paste(filename,".png",sep=""),units='in',res=res,height=height,width=width)
  #       if (filetype=='jpg' | filetype=='jpeg') jpeg(paste(filename,".jpg",sep=""),units='in',res=res,height=height,width=width)
  #       if (filetype=="svg")
  #       {
  #         if (R.Version()$arch=="x64") stop("RSVGTipsDevice is not available for 64bit versions of R.")
  #         require("RSVGTipsDevice")
  #         devSVGTips(paste(filename,".svg",sep=""),width=width,height=height,title=filename)
  #       }
  #       if (filetype=="tex")
  #       {
  #         # 	# Special thanks to Charlie Sharpsteen for supplying these tikz codes on stackoverflow.com !!!
  #         # 	
  #         # 	if (!suppressPackageStartupMessages(require(tikzDevice,quietly=TRUE))) stop("tikzDevice must be installed to use filetype='tex'")
  #         # 	opt= c( 
  #         # 	getOption('tikzLatexPackages'),  
  #         #     "\\def\\tooltiptarget{\\phantom{\\rule{1mm}{1mm}}}",
  #         #     "\\newbox\\tempboxa\\setbox\\tempboxa=\\hbox{}\\immediate\\pdfxform\\tempboxa \\edef\\emptyicon{\\the\\pdflastxform}",
  #         #     "\\newcommand\\tooltip[1]{\\pdfstartlink user{/Subtype /Text/Contents  (#1)/AP <</N \\emptyicon\\space 0 R >>}\\tooltiptarget\\pdfendlink}"
  #         # 	)
  #         # 	
  #         # 	place_PDF_tooltip <- function(x, y, text)
  #         # 	{
  #         # 
  #         # 		# Calculate coordinates
  #         # 		tikzX <- round(grconvertX(x, to = "device"), 2)
  #         # 		tikzY <- round(grconvertY(y, to = "device"), 2)
  #         # 		# Insert node
  #         # 		tikzAnnotate(paste(
  #         # 		"\\node at (", tikzX, ",", tikzY, ") ",
  #         # 		"{\\tooltip{", text, "}};",
  #         # 		sep = ''
  #         # 		))
  #         # 	  invisible()
  #         # 	}
  #         # 	
  #         # 	print("NOTE: Using 'tex' as filetype will take longer to run than other filetypes")
  #         # 	
  #         # 	tikzDevice:::tikz(paste(filename,".tex",sep=""), standAlone = standAlone, width=width, height=height, packages=opt)
  #         
  #         stop("Tikz device no longer supported due to removal from CRAN. Please see www.sachaepskamp.com/cigraph for a fix")
  #       }
  #     }	
  #if (!filetype%in%c('pdf','png','jpg','jpeg','svg','R','eps','tiff')) warning(paste("File type",filetype,"is not supported")) 
  
  # Specify background:
  if (is.null(background) && !DoNotPlot){
    background <- par("bg")
    if (background == "transparent") background <- "white"    
  } else {
    background <- "white"
  }
  if (isColor(bg)) background <- bg
  # Remove alpha:
  background <- col2rgb(background, alpha = TRUE)
  background <- rgb(background[1],background[2],background[3],background[4],maxColorValue=255)
  
  if (is.null(subplotbg)) subplotbg <- background
  
  if (isTRUE(edge.label.bg)) edge.label.bg <- background
  if(is.null(cigraphObject$Arguments[['label.color']])) {
    # if(is.null(cigraphObject$Arguments$lcolor)) lcolor <- ifelse(mean(col2rgb(background)/255) > 0.5,"black","white") else lcolor <- cigraphObject$Arguments$lcolor
    if(is.null(cigraphObject$Arguments$lcolor)) lcolor <- NA else lcolor <- cigraphObject$Arguments$lcolor
  } else lcolor <- cigraphObject$Arguments[['label.color']]
  
  # Legend setting 2
  if (legend & !is.null(scores))
  {
    layout(t(1:2),widths=c(GLratio,1))
  }
  
  # Weighted settings:
  if (is.null(weighted))
  {
    if (edgelist)
    {
      if (ncol(input)==2) weighted=FALSE else weighted=TRUE
    }
    if (!edgelist)
    {
      if (all(unique(c(input)) %in% c(0,1)) & !grepl("sig",mode)) weighted <- FALSE else weighted <- TRUE
    }
  }		
  if (!weighted) cut=0
  
  # par settings:
  #parOrig <- par(no.readonly=TRUE)
  if (!DoNotPlot)
  {
    par(pty=pty)
  }
  
  if (!edgelist)
  {
    if (!is.logical(directed)) if (is.null(directed))
    {
      if (!isSymmetric(unname(input))) directed=TRUE else directed=FALSE
    }
  }
  
  
  # Set default edge width:
  if(is.null(cigraphObject$Arguments[["esize"]])) 
  {
    if (weighted)
    {
      #       esize <- max((-1/72)*(nNodes)+5.35,2) 
      esize <- 15*exp(-nNodes/90)+1
    } else {
      esize <- 2
    }
    if (any(directed)) esize <- max(esize/2,1)
  } else esize <- cigraphObject$Arguments$esize
  
  # asize default:
  if(is.null(cigraphObject$Arguments[["asize"]]))
  {
    #       asize <- max((-1/10)*(nNodes)+4,1)
    #     asize <- ifelse(nNodes>10,2,3)
    asize <- 2*exp(-nNodes/20)+2
  } else asize <- cigraphObject$Arguments[["asize"]]
  
  if(!is.null(cigraphObject$Arguments[["edge.width"]]))
  {
    esize <- esize * cigraphObject$Arguments[["edge.width"]]
    asize <- asize * sqrt(cigraphObject$Arguments[["edge.width"]])
  }
  
  ## arrowAngle default:
  if(is.null(cigraphObject$Arguments[["arrowAngle"]])) 
  {
    if (weighted) arrowAngle <- pi/6 else arrowAngle <- pi/8
  } else {
    arrowAngle <- cigraphObject$Arguments[["arrowAngle"]]
  }
  
  
  ########### GRAPHICAL MODEL SELECTION #######
  
  if (graph == "cor") {
    if(!all(eigen(input)$values > 0))  {
      warning("Correlation/covariance matrix is not positive definite. Finding nearest positive definite matrix")
      
      input <- as.matrix(Matrix::nearPD(input, keepDiag = TRUE, ensureSymmetry = TRUE)$mat)
    }
  }
  
  
  # Partial graph:
  if (graph != "default")
  {
    if (edgelist) stop("Graph requires correlation or covariance matrix")
    
    # Check for symmetric matrix:
    if (!isSymmetric(input))
    {
      stop("Input matrix is not symmetric, thus can not be a correlation or covariance matrix.")
    }
    
    # Check for positive definiteness (glasso does its own check):
    if (graph != "glasso")
    {
      if(!all(eigen(input)$values > 0))  {
        warning("Correlation/covariance matrix is not positive definite. Finding nearest positive definite matrix")
        
        input <- as.matrix(Matrix::nearPD(input, keepDiag = TRUE, ensureSymmetry = TRUE)$mat)
      }
    }
    
    # Association graph:
    if (graph == "cor")
    {
      if (!all(diag(input) == 1)){
        input <- cov2cor(input)
      }
    }
    
    # Concentration graph:
    if (graph=="pcor") 
    {
      coln <- colnames(input)
      rown <- rownames(input)
      input <- cor2pcor(input)
      rownames(input) <- rown
      colnames(input) <- coln
    } 
    
    #     # FDR:
    #     if (tolower(graph)=="fdr.cor") 
    #     {
    #       if (!all(diag(input) == 1)){
    #         input <- cov2cor(input)
    #       }
    #       input <- FDRnetwork(input, FDRcutoff)
    #     } 
    #     
    #     if (tolower(graph)=="fdr.pcor") 
    #     {
    #       input <- cor2pcor(input)
    #       input <- FDRnetwork(input, FDRcutoff)
    #     } 
    #     
    #     if (tolower(graph) == "fdr")
    #     {
    #       input <- cor2pcor(input)
    #       testResult <- GeneNet::ggm.test.edges(input, fdr = TRUE, plot = FALSE)
    #       net <- GeneNet::extract.network(testResult)
    #       input <- matrix(0, nrow(input), ncol(input))
    #       for (i in seq_len(nrow(net)))
    #       {
    #         input[net$node1[i],net$node2[i]] <- input[net$node2[i],net$node1[i]] <- net$pcor[i]
    #       }
    #     }
    
    # Glasso graph:
    if (graph == "glasso")
    {
      if (edgelist) stop("Concentration graph requires correlation matrix")
      if (is.null(sampleSize)) stop("'sampleSize' argument is needed for glasso estimation")
      input <- EBICglasso(input, sampleSize, gamma = tuning,
                          refit=refit, lambda.min.ratio = lambda.min.ratio,
                          threshold = isTRUE(threshold))
    }
    
    if (graph == "ggmModSelect")
    {
      if (edgelist) stop("Concentration graph requires correlation matrix")
      if (is.null(sampleSize)) stop("'sampleSize' argument is needed for ggmModSelect estimation")
      input <- ggmModSelect(input, sampleSize, gamma = tuning, lambda.min.ratio = lambda.min.ratio)$graph
    }
    
    diag(input) <- 1
    input <- as.matrix(forceSymmetric(input))
  }
  
  
  ## Thresholding ####
  # If threshold is TRUE and graph = "glasso" or "ggmModSelect", set to FALSE:
  if (isTRUE(threshold) && (graph == "glasso" || graph == "ggmModSelect")){
    threshold <- cigraphObject$Arguments[['threshold']] <- 0
  }
    
    
  if (is.character(threshold))
  {    
    if (graph == "default")
    {
      if (verbose) message("'threshold' is assigned a string but 'graph' is not assigned. Detecting if input could be a correlation matrix.")
      # Detect if graph could be correlations or covariance matrix:
      
      # Check if input was a matrix:
      if (!is.matrix(input) | edgelist) stop(paste0("'",threshold,"' threshold requires a (partial) correlation/covariance matrix as input"))
      
      # Check if input is square matrix:
      if (!isSymmetric(input)) stop(paste0("'",threshold,"' threshold requires a (partial) correlation/covariance matrix as input: input was not a square matrix."))
      
      # Check if input is positive semi definite:
      if (any(eigen(input)$values < 0)) stop(paste0("'",threshold,"' threshold requires a (partial) correlation/covariance matrix as input: input was not a positive semi-definite matrix"))
      
      # If these checks are passed assume matrix is correlation or covariance: 
    } else {
      if (!graph %in% c("cor","pcor"))
      {
        stop("Thresholding by significance level only supported for graph = 'cor' or graph = 'pcor'")
      }
    }
    
    # Stop for incorrect threshold:
    if (!threshold %in% c('sig','holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none', 'locfdr'))
    {
      stop("'threshold' argument must be number or 'sig','holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none' or 'locfdr'")
    }
    
    # Significance:
    if (threshold != "locfdr")
    {
      if (grepl("sig",threshold,ignore.case=TRUE))
      {
        threshold <- "none"
      }
      
      if (is.null(sampleSize))
      {
        stop("'sampleSize' argument is needed for all thresholding with significance except 'locfdr'")
      }
      nadj <- sampleSize
      if (graph == "pcor")
      {
        nadj <- nadj - (nNodes - 2)
      }

      # Fix for col/row names bugs:
      if (is.null(colnames(input))){
        colnames(input) <- paste0("V",seq_len(ncol(input)))
      }
      if (is.null(rownames(input))){
        rownames(input) <- paste0("V",seq_len(ncol(input)))
      }
      
      # Compute p-values:
      if (all(diag(input)==1)) 
      {
        pvals <- psych::corr.p(input,n = nadj, adjust = threshold, alpha = max(alpha))$p
      } else {
        pvals <- psych::corr.p(cov2cor(input), n = nadj, adjust = threshold, alpha = max(alpha))$p
      }
      
      # Symmetrize:
      pvals[lower.tri(pvals)] <- t(pvals)[lower.tri(pvals)]
      
      # Remove insignificant edges:
      input <- input * (pvals < max(alpha))
    } else {
      input <- FDRnetwork(input, FDRcutoff)
    }
    
    threshold <- 0
  }
  
  
  #######################3
  
  ## diag default:
  if(is.null(cigraphObject$Arguments[['diag']])) 
  {
    if (edgelist) diag <- FALSE  else diag <- length(unique(diag(input))) > 1
  } else { 
    diag <- cigraphObject$Arguments$diag
  }
  
  # Diag:
  diagCols=FALSE
  diagWeights=0
  if (is.character(diag)) 
  {
    if (diag=="col" & !edgelist)
    {
      diagWeights=diag(input)
      diagCols=TRUE
      diag=FALSE
    }
  }
  if (is.numeric(diag))
  {
    if (length(diag)==1) diag=rep(diag,nNodes)
    if (length(diag)!=nNodes) stop("Numerical assignment of the 'diag' argument must be if length equal to the number of nodes")
    diagWeights=diag
    diagCols=TRUE
    diag=FALSE
  }
  if (is.logical(diag)) if (!diag & !edgelist) diag(input)=0
  
  # CREATE EDGELIST:
  
  E <- list()
  
  # Remove nonfinite weights:
  if (any(!is.finite(input)))
  {
    input[!is.finite(input)] <- 0
    warning("Non-finite weights are omitted")
  }
  
  if (edgelist)
  {
    E$from=input[,1]
    E$to=input[,2]
    if (ncol(input)>2) E$weight=input[,3] else E$weight=rep(1,length(E$from))
    if (length(directed)==1) directed=rep(directed,length(E$from))
    if (graph %in% c("sig","significance"))
    {
      if (sigSign)
      {
        E$weight <- sign0(E$weight) * fdrtool(E$weight,"correlation",plot=FALSE, color.figure=FALSE, verbose=FALSE)$pval
      } else E$weight <- fdrtool(E$weight,"correlation",plot=FALSE, color.figure=FALSE, verbose=FALSE)$pval
    }
    if (bonf)
    {
      if (mode=="sig") 
      {
        E$weight <- E$weight * length(E$weight)
        E$weight[E$weight > 1] <- 1
        E$weight[E$weight < -1] <- -1
      } # else warning("Bonferonni correction is only applied if mode='sig'")
    }
    if (mode=="sig" & any(E$weight < -1 | E$weight > 1))
    {
      warning("Weights under -1 set to -1 and weights over 1 set to 1")
      E$weight[E$weight< -1] <- -1
      E$weight[E$weight>1] <- 1
    }
    
    if (mode=="sig") 
    {
      Pvals <- E$weight
      E$weight <- sign0(E$weight) * sigScale(abs(E$weight))
    }
    if (OmitInsig)
    {
      #       if (!require("fdrtool")) stop("`fdrtool' package not found, is it installed?")
      if (mode != "sig") Pvals <- fdrtool(E$weight,"correlation",plot=FALSE, color.figure=FALSE, verbose=FALSE)$pval
      E$weight[abs(Pvals) > alpha[length(alpha)]] <- 0
    }
    
  } else
  {
    if (is.matrix(directed))
    {
      incl <- directed|upper.tri(input,diag=TRUE)
    } else
    {
      if (length(directed)>1) 
      {
        stop("'directed' must be TRUE or FALSE or a matrix containing TRUE or FALSE for each element of the input matrix") 
      } else
      { 
        if (directed)
        {
          incl <- matrix(TRUE,nNodes,nNodes)
        } else
        {
          if (isSymmetric(unname(input)))
          {
            
            incl <- upper.tri(input,diag=TRUE)
          } else 
          {
            incl <- matrix(TRUE,nNodes,nNodes)
          }
        }  
        directed <- matrix(directed,nNodes,nNodes)
      }
    }
    directed <- directed[incl]
    
    E$from=numeric(0)
    E$to=numeric(0)
    E$weight=numeric(0)
    
    E$from=rep(1:nrow(input),times=nrow(input))
    E$to=rep(1:nrow(input),each=nrow(input))
    E$weight=c(input)
    
    
    E$from <- E$from[c(incl)]
    E$to <- E$to[c(incl)]
    E$weight <- E$weight[c(incl)]
    if (graph %in% c("sig","significance"))
    {
      if (sigSign)
      {
        E$weight <- sign0(E$weight) * fdrtool(E$weight,"correlation",plot=FALSE, color.figure=FALSE, verbose=FALSE)$pval
      } else E$weight <- fdrtool(E$weight,"correlation",plot=FALSE, color.figure=FALSE, verbose=FALSE)$pval
    }
    if (bonf)
    {
      if (mode=="sig") 
      {
        E$weight <- E$weight * length(E$weight)
        E$weight[E$weight > 1] <- 1
        E$weight[E$weight < -1] <- -1
      } # else warning("Bonferonni correction is only applied if mode='sig'")
    }
    if (mode=="sig" & any(E$weight < -1 | E$weight > 1))
    {
      warning("Weights under -1 inputusted to -1 and weights over 1 input adjusted to 1")
      E$weight[E$weight < -1] <- -1
      E$weight[E$weight > 1] <- 1
    }
    
    if (mode=="sig") 
    {
      Pvals <- E$weight
      E$weight <- sign0(E$weight) * sigScale(abs(E$weight))
    }
    
    if (OmitInsig)
    {
      #       if (!require("fdrtool")) stop("`fdrtool' package not found, is it installed?")
      if (mode != "sig") Pvals <- fdrtool(E$weight,"correlation",plot=FALSE, color.figure=FALSE, verbose=FALSE)$pval
      E$weight[abs(Pvals) > alpha[length(alpha)]] <- 0
    }	
    if (is.list(knots))
    {
      knotList <- knots
      knots <- matrix(0,nNodes,nNodes)
      for (k in seq_along(knotList))
      {
        knots[knotList[[k]]] <- k
      }
      # If undirected, symmetrize:
      if (all(incl[upper.tri(incl,diag=TRUE)]) & !any(incl[lower.tri(incl)]))
      {
        knots <- pmax(knots,t(knots))
      }
    }
    if (is.matrix(knots))
    {
      knots <- knots[c(incl)]
      #       knots <- knots[E$weight!=0]
    }
    if (is.matrix(curve))
    {
      curve <- curve[c(incl)]
      #       curve <- curve[E$weight!=0]
    }
    if (is.matrix(parallelEdge))
    {
      parallelEdge <- parallelEdge[c(incl)]
      #       parallelEdge <- parallelEdge[E$weight!=0]
    }
    if (is.matrix(parallelAngle))
    {
      parallelAngle <- parallelAngle[c(incl)]
      #       parallelAngle <- parallelAngle[E$weight!=0]
    }
    if (is.matrix(bidirectional))
    {
      bidirectional <- bidirectional[c(incl)]
      #       bidirectional <- bidirectional[E$weight!=0]
    }
    if (is.matrix(residEdge))
    {
      residEdge <- residEdge[c(incl)]
      #       residEdge <- residEdge[E$weight!=0]
    }
    if (is.matrix(CircleEdgeEnd))
    {
      CircleEdgeEnd <- CircleEdgeEnd[c(incl)]
      #       CircleEdgeEnd <- CircleEdgeEnd[E$weight!=0]
    }      
    if (is.matrix(edge.labels))
    {
      edge.labels <- edge.labels[c(incl)]
      #       edge.labels <- edge.labels[E$weight!=0]
    }
    if (is.matrix(edge.color))
    {
      edge.color <- edge.color[c(incl)]
      #       edge.color <- edge.color[E$weight!=0]
    }
    if (is.matrix(edge.label.bg))
    {
      edge.label.bg <- edge.label.bg[c(incl)]
      #       edge.label.bg <- edge.label.bg[E$weight!=0]
    }
    if (is.matrix(edge.label.margin))
    {
      edge.label.margin <- edge.label.margin[c(incl)]
      #       edge.label.bg <- edge.label.bg[E$weight!=0]
    }
    if (is.matrix(edge.label.font))
    {
      edge.label.font <- edge.label.font[c(incl)]
      #       edge.label.font <- edge.label.font[E$weight!=0]
    }
    if (is.matrix(fade))
    {
      fade <- fade[c(incl)]
      #       edge.color <- edge.color[E$weight!=0]
    }
    if (!is.null(ELcolor))
    {
      if (is.matrix(ELcolor))
      {
        ELcolor <- ELcolor[c(incl)]
        #         ELcolor <- ELcolor[E$weight!=0]
      }      
    }
    
    # if (!is.null(edge.color)) if (length(edge.color) == length(E$weight)) edge.color <- edge.color[E$weight!=0]
    
    if (is.matrix(lty))
    {
      lty <- lty[c(incl)]
      #       lty <- lty[E$weight!=0]
    }
    
    if (!is.null(edgeConnectPoints))
    {
      if (is.array(edgeConnectPoints) && isTRUE(dim(edgeConnectPoints)[3]==2))
      {
        edgeConnectPoints <- matrix(edgeConnectPoints[c(incl,incl)],,2)
        #         edgeConnectPoints <- edgeConnectPoints[E$weight!=0,,drop=FALSE]
      }
    }
    
    if (is.matrix(edge.label.position))
    {
      edge.label.position <- edge.label.position[c(incl)]
      #       edge.label.position <- edge.label.position[E$weight!=0]
    }
  }	

  keep <- abs(E$weight)>threshold
  
  ######
  if (length(loopRotation)==1) loopRotation <- rep(loopRotation,nNodes)
  
  if (length(directed)==1) 
  {
    directed <- rep(directed,length(E$from))
  }
  directed <- directed[keep]
  
  if (!is.null(edge.color) && length(edge.color) != sum(keep)) 
  {
    edge.color <- rep(edge.color,length=length(E$from))
    if (length(edge.color) != length(keep)) stop("'edge.color' is wrong length")
    edge.color <- edge.color[keep]
  }
  
  if (!is.logical(edge.labels))
  {
    if (length(edge.labels) == 1) edge.labels <- rep(edge.labels,length(E$from))
    if (length(edge.labels) != length(keep) & length(edge.labels) != sum(keep)) stop("'edge.label.bg' is wrong length")
    if (length(edge.labels)==length(keep)) edge.labels <- edge.labels[keep]
    
    # edge.labels <- rep(edge.labels,length=length(E$from))
  }
  
  #     if (is.logical(edge.label.bg))
  #     {
  #       edge.label.bg <- "white"
  #     }
  if (length(edge.label.bg) == 1) edge.label.bg <- rep(edge.label.bg,length(E$from))
  if (length(edge.label.bg) != length(keep) & length(edge.label.bg) != sum(keep)) stop("'edge.label.bg' is wrong length")
  if (length(edge.label.bg)==length(keep)) edge.label.bg <- edge.label.bg[keep]
  
  #     }
  if (length(edge.label.margin) == 1) edge.label.margin <- rep(edge.label.margin,length(E$from))
  if (length(edge.label.margin) != length(keep) & length(edge.label.margin) != sum(keep)) stop("'edge.label.margin' is wrong length")
  if (length(edge.label.margin)==length(keep)) edge.label.margin <- edge.label.margin[keep]
  
  
  
  if (length(edge.label.font) == 1) edge.label.font <- rep(edge.label.font,length(E$from))
  if (length(edge.label.font) != length(keep) & length(edge.label.font) != sum(keep)) stop("'edge.label.font' is wrong length")
  if (length(edge.label.font)==length(keep)) edge.label.font <- edge.label.font[keep]
  
  
  if (length(lty) == 1) lty <- rep(lty,length(E$from))
  if (length(lty) != length(keep) & length(lty) != sum(keep)) stop("'lty' is wrong length")
  if (length(lty)==length(keep)) lty <- lty[keep]
  
  if (length(fade) == 1) fade <- rep(fade,length(E$from))
  if (length(fade) != length(keep) & length(fade) != sum(keep)) stop("'fade' is wrong length")
  if (length(fade)==length(keep)) fade <- fade[keep]
  
  
  if (!is.null(edgeConnectPoints))
  {
    if (length(edgeConnectPoints) == 1) edgeConnectPoints <- matrix(rep(edgeConnectPoints,2*length(E$from)),,2)
    if (nrow(edgeConnectPoints) != length(keep) & nrow(edgeConnectPoints) != sum(keep)) stop("Number of rows in 'edgeConnectPoints' do not match number of edges")
    if (nrow(edgeConnectPoints)==length(keep)) edgeConnectPoints <- edgeConnectPoints[keep,,drop=FALSE]
  }
  
  if (length(edge.label.position) == 1) edge.label.position <- rep(edge.label.position,length(E$from))
  if (length(edge.label.position) != length(keep) & length(edge.label.position) != sum(keep)) stop("'edge.label.position' is wrong length")
  if (length(edge.label.position)==length(keep)) edge.label.position <- edge.label.position[keep]  
  
  
  if (!is.null(ELcolor))
  {
    ELcolor <- rep(ELcolor,length = length(E$from))
    ELcolor <- ELcolor[keep]    
  }
  
  
  if (is.list(knots))
  {
    knotList <- knots
    knots <- rep(0,length(E$from))
    for (k in seq_along(knotList))
    {
      knots[knotList[[k]]] <- k
    }
  }
  if (length(knots)==length(keep)) knots <- knots[keep]
  
  if (length(bidirectional)==1) 
  {
    bidirectional <- rep(bidirectional,length(E$from))
  }
  if (length(bidirectional)==length(keep)) bidirectional <- bidirectional[keep]
  if (length(residEdge)==1) 
  {
    residEdge <- rep(residEdge,length(E$from))
  }
  if (length(residEdge)==length(keep)) residEdge <- residEdge[keep]    
  
  if (length(CircleEdgeEnd)==1) 
  {
    CircleEdgeEnd <- rep(CircleEdgeEnd,length(E$from))
  }
  if (length(CircleEdgeEnd)==length(keep)) CircleEdgeEnd <- CircleEdgeEnd[keep]    
  
  if (!is.logical(edge.labels))
  {
    if (length(edge.labels)==length(keep))
    {
      edge.labels <- edge.labels[keep]
    }
  }
  
  if (length(curve)==1) 
  {
    curve <- rep(curve,length(E$from))
  }
  if (length(curve)==length(keep)) curve <- curve[keep]   
  
  
  if (length(parallelEdge)==1) 
  {
    parallelEdge <- rep(parallelEdge,length(E$from))
  }
  if (length(parallelEdge)==length(keep)) parallelEdge <- parallelEdge[keep]    
  
  
  if (length(parallelAngle)==1) 
  {
    parallelAngle <- rep(parallelAngle,length(E$from))
  }
  if (length(parallelAngle)==length(keep)) parallelAngle <- parallelAngle[keep]    
  
  E$from=E$from[keep]
  E$to=E$to[keep]
  if (mode=="sig") Pvals <- Pvals[keep]
  E$weight=E$weight[keep]
  
  
  ## Define cut:
  if (defineCut)
  {
    
    if (length(E$weight) > 3*nNodes)
    {
      #       cut <- median(sort(E$weight,decreasing=TRUE)[seq_len(nNodes)])
      cut <- max(sort(abs(E$weight),decreasing=TRUE)[2*nNodes], quantile(abs(E$weight),0.75))
    } else if (length(E$weight) > 1) cut <- quantile(abs(E$weight),0.75) else cut <- 0
    #     cut <- quantile(abs(E$weight), cutQuantile)
  }
  
  
  if (length(E$from) > 0)
  {
    maximum=max(abs(c(maximum,max(abs(E$weight)),cut,abs(diagWeights))))
  } else maximum = 1
  if (cut==0)
  {
    avgW=(abs(E$weight)-minimum)/(maximum-minimum)
  } else if (maximum>cut) avgW=(abs(E$weight)-cut)/(maximum-cut) else avgW=rep(0,length(E$from))
  avgW[avgW<0]=0
  
  
  edgesort=sort(abs(E$weight),index.return=TRUE)$ix
  edge.width=rep(1,length(E$weight))
  
  
  # lty and curve settings:
  
  if (length(lty)==1) lty=rep(lty,length(E$from))
  
  if (length(edge.label.position)==1) edge.label.position=rep(edge.label.position,length(E$from))
  
  
  # Make bidirectional vector:
  if (length(bidirectional)==1) bidirectional=rep(bidirectional,length(E$from))
  if (length(bidirectional)!=length(E$from)) stop("Bidirectional vector must be of length 1 or equal to the number of edges")
  
  srt <- cbind(pmin(E$from,E$to), pmax(E$from,E$to) , knots, abs(E$weight) > minimum)
  
  if (!curveAll | any(parallelEdge))
  {
    dub <- duplicated(srt)|duplicated(srt,fromLast=TRUE)
    
    if (!curveAll)
    {
      if (length(curve)==1) curve <- rep(curve,length(E$from))
      curve <- ifelse(is.na(curve),ifelse(knots==0&dub&!bidirectional&is.na(curve),ifelse(E$from==srt[,1],1,-1) * ave(1:nrow(srt),srt[,1],srt[,2],bidirectional,FUN=function(x)seq(curveDefault,-curveDefault,length=length(x))),0),curve)
    }
    
    if (any(parallelEdge))
    {
      # Set parallelAngle value:   
      parallelAngle <- ifelse(is.na(parallelAngle),ifelse(knots==0&dub&!bidirectional&is.na(parallelAngle),ifelse(E$from==srt[,1],1,-1) * ave(1:nrow(srt),srt[,1],srt[,2],bidirectional,FUN=function(x)seq(parallelAngleDefault,-parallelAngleDefault,length=length(x))),0),parallelAngle) 
    }
    
    rm(dub)
  }
  
  parallelAngle[is.na(parallelAngle)] <- 0
  
  
  # Layout settings:
  if (nNodes == 1 & isTRUE(rescale))
  {
    layout <- matrix(0,1,2)
  } else {
    if (is.null(layout)) layout="default"
    
    if (!is.matrix(layout))
    {
      # If function, assume igraph function (todo: check this)
      if (is.function(layout))
      {
        Graph <- graph.edgelist(as.matrix(cbind(E$from,E$to)), any(directed))
        E(Graph)$weight <- E$weight
        
        # set roots:
        if (deparse(match.call()[['layout']]) == "layout.reingold.tilford" && is.null(layout.par[['root']]))
        {
          sp <- shortest.paths(Graph, mode = "out")
          diag(sp) <- Inf
          
          # Find root nodes:
          roots <- which(colSums(sp==Inf) == nrow(sp))
          # Find roots with longest outgoing paths:
          maxs <- sapply(roots,function(x)max(sp[x,sp[x,]!=Inf]))
          
          layout.par[['root']] <- roots[maxs==max(maxs)]
        }
        
        layout <- do.call(layout,c(list(graph = Graph),layout.par))
      } else {
        
        if (length(layout) > 1) stop("Incorrect specification of layout.")
        if (layout=="default" & (any(directed) | !weighted)) layout="spring"
        if (layout=="default" | layout=="circular" | layout=="circle" | layout=="groups") 
        {
          if (is.null(groups) | layout == "circle")
          {
            layout=matrix(0,nrow=nNodes,ncol=2)
            tl=nNodes+1
            layout[,1]=sin(seq(0,2*pi, length=tl))[-tl]
            layout[,2]=cos(seq(0,2*pi, length=tl))[-tl] 
          } else
          {
            if (is.null(rotation)) rotation=rep(0,length=length(groups))
            
            l1=matrix(0,nrow=length(groups),ncol=2)
            tl=nrow(l1)+1
            l1[,1]=sin(seq(0,2*pi, length=tl))[-tl]
            l1[,2]=cos(seq(0,2*pi, length=tl))[-tl]
            l1=l1*length(groups)*layout.control
            
            layout=matrix(0,nrow=nNodes,ncol=2)
            for (i in 1:length(groups)) 
            {
              tl=length(groups[[i]])+1
              layout[groups[[i]],1]=repulsion*sin(seq(rotation[i],rotation[i]+2*pi, length=tl))[-tl]+l1[i,1]
              layout[groups[[i]],2]=repulsion*cos(seq(rotation[i],rotation[i]+2*pi, length=tl))[-tl]+l1[i,2] 
            }
          }
        } else if (layout=="spring")
        {
          
          if (length(E$weight) > 0)
          {
            if (mode != "sig")
            {
              layout=cigraph.layout.fruchtermanreingold(cbind(E$from,E$to),abs(E$weight/max(abs(E$weight)))^2,nNodes,rotation=rotation,layout.control=layout.control,
                                                       niter=layout.par$niter,max.delta=layout.par$max.delta,area=layout.par$area,cool.exp=layout.par$cool.exp,repulse.rad=layout.par$repulse.rad,init=layout.par$init,
                                                       constraints=layout.par$constraints)
            } else
            {
              layout=cigraph.layout.fruchtermanreingold(cbind(E$from,E$to),abs(E$weight),nNodes,rotation=rotation,layout.control=layout.control,
                                                       niter=layout.par$niter,max.delta=layout.par$max.delta,area=layout.par$area,cool.exp=layout.par$cool.exp,repulse.rad=layout.par$repulse.rad,init=layout.par$init,
                                                       constraints=layout.par$constraints)
            }
          } else
          {
            if (mode != "sig")
            {
              layout=cigraph.layout.fruchtermanreingold(cbind(E$from,E$to),numeric(0),nNodes,rotation=rotation,layout.control=layout.control,
                                                       niter=layout.par$niter,max.delta=layout.par$max.delta,area=layout.par$area,cool.exp=layout.par$cool.exp,repulse.rad=layout.par$repulse.rad,init=layout.par$init,
                                                       constraints=layout.par$constraints)
            } else
            {
              layout=cigraph.layout.fruchtermanreingold(cbind(E$from,E$to),numeric(0),nNodes,rotation=rotation,layout.control=layout.control,
                                                       niter=layout.par$niter,max.delta=layout.par$max.delta,area=layout.par$area,cool.exp=layout.par$cool.exp,repulse.rad=layout.par$repulse.rad,init=layout.par$init,
                                                       constraints=layout.par$constraints)
            }
          }
        } 
      }
    }
    # Layout matrix:
    if (is.matrix(layout)) if (ncol(layout)>2)
    {
      layout[is.na(layout)] <- 0
      # If character and labels exist, replace:
      if (is.character(layout) && is.character(labels))
      {
        layout[] <- match(layout,labels)
        layout[is.na(layout)] <- 0
        mode(layout) <- 'numeric'
      }
      
      # Check:
      if (!all(seq_len(nNodes) %in% layout)) stop("Grid matrix does not contain a placement for every node.")
      if (any(sapply(seq_len(nNodes),function(x)sum(layout==x))>1)) stop("Grid matrix contains a double entry.")
      
      Lmat=layout
      LmatX=seq(-1,1,length=ncol(Lmat))
      LmatY=seq(1,-1,length=nrow(Lmat))
      layout=matrix(0,nrow=nNodes,ncol=2)
      
      loc <- t(sapply(1:nNodes,function(x)which(Lmat==x,arr.ind=T)))
      layout <- cbind(LmatX[loc[,2]],LmatY[loc[,1]])
      
    }
  }
  
  # Rescale layout:
  l=original.layout=layout
  if (rescale) {
    if (aspect)
    {
      # center:
      l[,1] <- l[,1] - mean(l[,1])
      l[,2] <- l[,2] - mean(l[,2])
      lTemp <- l
      
      if (length(unique(lTemp[,1]))>1)
      {
        l[,1]=(lTemp[,1]-min(lTemp))/(max(lTemp)-min(lTemp))*2-1
      } else l[,1] <- 0
      if (length(unique(lTemp[,2]))>1)
      {
        l[,2]=(lTemp[,2]-min(lTemp))/(max(lTemp)-min(lTemp))*2-1 
      } else l[,2] <- 0
      
      rm(lTemp)
      
      
      #     # Equalize white space:
      #     if (diff(range(l[,1])) < 2)
      #     {
      #       l[,1] <- diff(range(l[,1]))/2 + l[,1]
      #     }
      #     if (diff(range(l[,2])) < 2)
      #     {
      #       l[,2] <- (2-diff(range(l[,2])))/2 + l[,2]
      #     }
      
      layout=l    
    } else
    {
      if (length(unique(l[,1]))>1)
      {
        l[,1]=(l[,1]-min(l[,1]))/(max(l[,1])-min(l[,1]))*2-1
      } else l[,1] <- 0
      if (length(unique(l[,2]))>1)
      {
        l[,2]=(l[,2]-min(l[,2]))/(max(l[,2])-min(l[,2]))*2-1 
      } else l[,2] <- 0
      layout=l
    }
  }
  
  ## Offset and scale:
  if (length(layoutScale) == 1) layoutScale <- rep(layoutScale,2)
  if (length(layoutOffset) == 1) layoutOffset <- rep(layoutOffset,2)
  layout[,1] <- layout[,1] * layoutScale[1] + layoutOffset[1]
  layout[,2] <- layout[,2] * layoutScale[2] + layoutOffset[2]
  l <- layout
  
  
  # Set Edge widths:
  if (mode=="direct")
  {
    edge.width <- abs(E$weight)
  } else
  {
    if (weighted)
    {
      edge.width <- avgW*(esize-1)+1
      edge.width[edge.width<1]=1
    } else {
      edge.width <- rep(esize,length(E$weight))
    }
  }
  
  #     # Set edge colors:
  #     if (is.null(edge.color) || (any(is.na(edge.color)) || fade))
  #     {
  #       if (!is.null(edge.color))
  #       {
  #         repECs <- TRUE
  #         ectemp <- edge.color
  #       } else  repECs <- FALSE
  #       
  #       col <- rep(1,length(E$from))
  #       
  #       if (weighted) 
  #       {
  #         #Edge color:
  #         edge.color=rep("#00000000",length(E$from))
  #         
  #         
  #         if (mode=="strength"|mode=="direct")
  #         {
  #           if (cut==0) 
  #           {
  #             col=(abs(E$weight)-minimum)/(maximum-minimum)
  #           } else 
  #           {
  #             col=(abs(E$weight)-minimum)/(cut-minimum)
  #           }
  #           col[col>1]=1
  #           col[col<0]=0
  #           if (!gray)
  #           {
  #             if (transparency) 
  #             {
  #               col=col^(2)
  #               neg=col2rgb(rgb(0.75,0,0))/255
  #               pos=col2rgb(rgb(0,0.6,0))/255
  #               
  #               # Set colors for edges over cutoff:
  #               edge.color[E$weight< -1* minimum] <- rgb(neg[1],neg[2],neg[3],col[E$weight< -1*minimum])
  #               edge.color[E$weight> minimum] <- rgb(pos[1],pos[2],pos[3],col[E$weight> minimum])
  #             } else 
  #             {
  #               edge.color[E$weight>minimum]=rgb(1-col[E$weight > minimum],1-(col[E$weight > minimum]*0.25),1-col[E$weight > minimum])
  #               edge.color[E$weight< -1*minimum]=rgb(1-(col[E$weight < (-1)*minimum]*0.25),1-col[E$weight < (-1)*minimum],1-col[E$weight < (-1)*minimum])
  #             }	
  #           } else
  #           {
  #             if (transparency) 
  #             {
  #               col=col^(2)
  #               neg="gray10"
  #               pos="gray10"
  #               
  #               # Set colors for edges over cutoff:
  #               edge.color[E$weight< -1* minimum] <- rgb(neg[1],neg[2],neg[3],col[E$weight< -1*minimum])
  #               edge.color[E$weight> minimum] <- rgb(pos[1],pos[2],pos[3],col[E$weight> minimum])
  #             } else 
  #             {
  #               edge.color[E$weight>minimum]=rgb(1-col[E$weight > minimum],1-(col[E$weight > minimum]),1-col[E$weight > minimum])
  #               edge.color[E$weight< -1*minimum]=rgb(1-(col[E$weight < (-1)*minimum]),1-col[E$weight < (-1)*minimum],1-col[E$weight < (-1)*minimum])
  #             }
  #           }
  #         }
  #         if (mode == "sig")
  #         {	
  #           
  #           if (!gray)
  #           {
  #             
  #             # Set colors for edges over sig > 0.01 :
  #             if (length(alpha) > 3) edge.color[Pvals > 0 & Pvals < alpha[4]  & E$weight > minimum] <- "cadetblue1"	
  #             # Set colors for edges over sig > 0.01 :
  #             if (length(alpha) > 2) edge.color[Pvals > 0 & Pvals < alpha[3]  & E$weight > minimum] <- "#6495ED"
  #             # Set colors for edges over sig > 0.01 :
  #             if (length(alpha) > 1) edge.color[Pvals > 0 & Pvals < alpha[2]  & E$weight > minimum] <- "blue"				
  #             # Set colors for edges over sig < 0.01 :
  #             edge.color[Pvals > 0 & Pvals < alpha[1]  & E$weight > minimum] <- "darkblue"
  #             
  #             # Set colors for edges over sig > 0.01 :
  #             if (length(alpha) > 3) edge.color[Pvals < 0 & Pvals > (-1 * alpha[4])  & E$weight < -1 * minimum] <- rgb(1,0.8,0.4) 	
  #             # Set colors for edges over sig > 0.01 :
  #             if (length(alpha) > 2) edge.color[Pvals < 0 & Pvals > (-1 * alpha[3])  & E$weight < -1 * minimum] <- "orange"
  #             # Set colors for edges over sig > 0.01 :
  #             if (length(alpha) > 1) edge.color[Pvals < 0 & Pvals > (-1 * alpha[2])  & E$weight < -1 * minimum] <- "darkorange"				
  #             # Set colors for edges over sig < 0.01 :
  #             edge.color[Pvals < 0 & Pvals > (-1 * alpha[1])  & E$weight < -1 * minimum] <- "darkorange2"
  #             
  #             
  #             
  #             
  #           } else
  #           {
  #             Pvals <- abs(Pvals)
  #             # Set colors for edges over sig < 0.01 :
  #             if (length(alpha) > 3) edge.color[Pvals > 0 & Pvals < alpha[4]  & E$weight > minimum] <- rgb(0.7,0.7,0.7)
  #             if (length(alpha) > 2) edge.color[Pvals > 0 & Pvals < alpha[3]  & E$weight > minimum] <- rgb(0.5,0.5,0.5)
  #             if (length(alpha) > 1) edge.color[Pvals > 0 & Pvals < alpha[2]  & E$weight > minimum] <- rgb(0.3,0.3,0.3)
  #             edge.color[Pvals > 0 & Pvals < alpha[1]  & E$weight > minimum] <- "black"
  #             
  #           }
  #         }
  #         if (cut!=0)
  #         {
  #           if (!gray & (mode=="strength"|mode=="direct"))
  #           {
  #             # Set colors for edges over cutoff:
  #             edge.color[E$weight<= -1*cut] <- "red"
  #             edge.color[E$weight>= cut] <- "darkgreen"
  #           } else if (gray)
  #           {
  #             # Set colors for edges over cutoff:
  #             edge.color[E$weight<= -1*cut] <- "black"
  #             edge.color[E$weight>= cut] <- "black"
  #             
  #           }
  #         }
  #         
  #       } else
  #       {
  #         if (!is.logical(transparency)) Trans=transparency else Trans=1
  #         edge.color=rep(rgb(0.5,0.5,0.5,Trans),length(edgesort))
  #       }
  #       if (repECs)
  #       {
  #         ## Add trans:
  #         if (fade & any(!is.na(ectemp)))
  #         {
  #           if (!is.logical(transparency)) col <- rep(transparency,length(col))
  #           edge.color[!is.na(ectemp)] <- addTrans(ectemp[!is.na(ectemp)],round(255*col[!is.na(ectemp)]))
  #         } else {
  #           edge.color[!is.na(ectemp)] <- ectemp[!is.na(ectemp)]
  #         }
  #         rm(ectemp)
  #       }
  #     } else {
  #       if (length(edge.color) == 1) edge.color <- rep(edge.color,length(E$from))
  #       if (length(edge.color) != length(E$from)) stop("Number of edge colors not equal to number of edges")
  #     }
  
  
  # Set edge colors:
  if (is.null(edge.color) || (any(is.na(edge.color)) || any(is.na(fade)) || any(fade != 1)))
  {
    if (!is.null(edge.color))
    {
      repECs <- TRUE
      ectemp <- edge.color
    } else  repECs <- FALSE
    
    # col vector will contain relative strength:
    col <- rep(1,length(E$from))
    
    if (weighted) 
    {
      # Dummmy vector containing invisible edges:
      edge.color <- rep("#00000000",length(E$from))
      
      # Normal color scheme (0 is invisible, stronger is more visible)
      if (mode=="strength"|mode=="direct")
      {
        # Set relative strength:
        if (cut==0) 
        {
          col <- (abs(E$weight)-minimum)/(maximum-minimum)
        } else 
        {
          if (cut > minimum){
            col <- (abs(E$weight)-minimum)/(cut-minimum)  
          } else {
            col <- ifelse(abs(E$weight) > minimum, 1, 0)
          }
          
        }
        col[col>1] <- 1
        col[col<0] <- 0
        col <- col^colFactor      
        
        # Set edges between minimum and cut:

        # if (autoFade)
        # {
          if (isTRUE(transparency))
          {
            edge.color[E$weight > minimum] <- addTrans(posCol[1],round(ifelse(is.na(fade),col,fade)[E$weight > minimum]*255))
            edge.color[E$weight < -1*minimum] <- addTrans(negCol[1],round(ifelse(is.na(fade),col,fade)[E$weight < -1*minimum]*255))
          } else {
            edge.color[E$weight > minimum] <- Fade(posCol[1],ifelse(is.na(fade),col,fade)[E$weight > minimum], background)
            edge.color[E$weight < -1*minimum] <- Fade(negCol[1],ifelse(is.na(fade),col,fade)[E$weight < -1*minimum], background)
          }
        # }
      # else {
      #     if (isTRUE(transparency))
      #     {
      #       edge.color[E$weight > minimum] <- addTrans(posCol[1],round(fade[E$weight > minimum]*255))
      #       edge.color[E$weight < -1*minimum] <- addTrans(negCol[1],round(fade[E$weight < -1*minimum]*255))
      #     } else {
      #       edge.color[E$weight > minimum] <- Fade(posCol[1],fade[E$weight > minimum], background)
      #       edge.color[E$weight < -1*minimum] <- Fade(negCol[1],fade[E$weight < -1*minimum], background)
      #     }
          
          # edge.color[E$weight > minimum] <- posCol[1]
          # edge.color[E$weight < -1*minimum] <- negCol[1]
        # }
        
        # Set colors over cutoff if cut != 0:
        if (cut!=0)
        {
          # Old code:
          # if (posCol[1]!=posCol[2]) edge.color[E$weight >= cut] <- posCol[2]
          # if (negCol[1]!=negCol[2]) edge.color[E$weight <= -1*cut] <- negCol[2]
          
          # New code (1.9.7)
          edge.color[E$weight >= cut & abs(E$weight) >= minimum] <- posCol[2]
          edge.color[E$weight <= -1*cut & abs(E$weight) >= minimum] <- negCol[2]
        }
      } 
      
      if (mode == "sig")
      {	
        if (!gray)
        {
          
          # Set colors for edges over sig > 0.01 :
          if (length(alpha) > 3) edge.color[Pvals >= 0 & Pvals < alpha[4]  & E$weight > minimum] <- "cadetblue1"	
          # Set colors for edges over sig > 0.01 :
          if (length(alpha) > 2) edge.color[Pvals >= 0 & Pvals < alpha[3]  & E$weight > minimum] <- "#6495ED"
          # Set colors for edges over sig > 0.01 :
          if (length(alpha) > 1) edge.color[Pvals >= 0 & Pvals < alpha[2]  & E$weight > minimum] <- "blue"				
          # Set colors for edges over sig < 0.01 :
          edge.color[Pvals >= 0 & Pvals < alpha[1]  & E$weight > minimum] <- "darkblue"
          
          # Set colors for edges over sig > 0.01 :
          if (length(alpha) > 3) edge.color[Pvals < 0 & Pvals > (-1 * alpha[4])  & E$weight < -1 * minimum] <- rgb(1,0.8,0.4) 	
          # Set colors for edges over sig > 0.01 :
          if (length(alpha) > 2) edge.color[Pvals < 0 & Pvals > (-1 * alpha[3])  & E$weight < -1 * minimum] <- "orange"
          # Set colors for edges over sig > 0.01 :
          if (length(alpha) > 1) edge.color[Pvals < 0 & Pvals > (-1 * alpha[2])  & E$weight < -1 * minimum] <- "darkorange"				
          # Set colors for edges over sig < 0.01 :
          edge.color[Pvals < 0 & Pvals > (-1 * alpha[1])  & E$weight < -1 * minimum] <- "darkorange2"
          
          
          
          
        } else
        {
          Pvals <- abs(Pvals)
          # Set colors for edges over sig < 0.01 :
          if (length(alpha) > 3) edge.color[Pvals > 0 & Pvals < alpha[4]  & E$weight > minimum] <- rgb(0.7,0.7,0.7)
          if (length(alpha) > 2) edge.color[Pvals > 0 & Pvals < alpha[3]  & E$weight > minimum] <- rgb(0.5,0.5,0.5)
          if (length(alpha) > 1) edge.color[Pvals > 0 & Pvals < alpha[2]  & E$weight > minimum] <- rgb(0.3,0.3,0.3)
          edge.color[Pvals > 0 & Pvals < alpha[1]  & E$weight > minimum] <- "black"
          
        }
      }
    } else
    {
      if (!is.logical(transparency)) Trans <- transparency else Trans <- 1
      edge.color <- rep(addTrans(unCol,round(255*Trans)),length(edgesort))
    }
    if (repECs)
    {
      # Colors to fade:
      indx <- !is.na(ectemp) & is.na(fade)
      
      ## Add trans:
      if (any(is.na(fade)) & any(!is.na(ectemp)))
      {
        # Replace all edge colors:
        edge.color[!is.na(ectemp)] <- ectemp[!is.na(ectemp)]
        

        
        if (!is.logical(transparency)) col <- rep(transparency,length(col))
        if (isTRUE(transparency))
        {
          edge.color[indx] <- addTrans(ectemp[indx],round(255*col[indx]))
        } else {
          edge.color[indx] <- Fade(ectemp[indx],col[indx], background)
        }
      } else {
        edge.color[indx] <- ectemp[indx]
      }
      rm(ectemp)
      rm(indx)
    }
  } else {
    if (length(edge.color) == 1) edge.color <- rep(edge.color,length(E$from))
    if (length(edge.color) != length(E$from)) stop("Number of edge colors not equal to number of edges")
  }
  
  # Vertex color:
  # if (is.null(color) & !is.null(groups))
  # {
  #   if (!gray) 
  #   {
  #     if (pastel)
  #     {
  #       color <- rainbow_hcl(length(groups), start = rainbowStart * 360, end = (360 * rainbowStart + 360*(length(groups)-1)/length(groups)))
  #     } else {
  #       color <- rainbow(length(groups), start = rainbowStart, end = (rainbowStart + (max(1.1,length(groups)-1))/length(groups)) %% 1)   
  #     }
  #   }
  #   if (gray) color <- sapply(seq(0.2,0.8,length=length(groups)),function(x)rgb(x,x,x))
  # }
  if (is.null(color) & !is.null(groups))
  {
    if (is.function(palette)){
      color <- palette(length(groups))
    } else if (palette == "rainbow"){
      color <- rainbow(length(groups), start = rainbowStart, end = (rainbowStart + (max(1.1,length(groups)-1))/length(groups)) %% 1)   
    } else if (palette == "gray" | palette == "grey"){
      color <- shadesOfGrey(length(groups))
    } else  if (palette == "colorblind"){
      color <- colorblind(length(groups))
    } else if (palette == "R"){
      color <- seq_len(length(groups))
    } else if (palette == "ggplot2"){
      color <- ggplot_palette(length(groups))
    } else if (palette == "pastel"){
      color <- rainbow_hcl(length(groups), start = rainbowStart * 360, end = (360 * rainbowStart + 360*(length(groups)-1)/length(groups)))
    } else if (palette == "neon"){
      color <- neon(length(groups))
    } else if (palette == "pride"){
      if (length(groups) > 7){
        color <- rainbow(length(groups), start = rainbowStart, end = (rainbowStart + (max(1.1,length(groups)-1))/length(groups)) %% 1)   
      } else {
        pridecols <- c("#E50000","#FF8D00","#FFEE00","#028121","#004CFF","#760088")
        # Reorder:
        startcol <- round(1 + rainbowStart * 6)
        sequence <- startcol:(startcol+length(groups))%%6
        sequence[sequence==0] <- 6
        color <- pridecols[sequence]        
      }

    }  else stop(paste0("Palette '",palette,"' is not supported."))
  }
  
  # Default color:
  if (is.null(color))	color <- "background"  
  


  vertex.colors <- rep(color, length=nNodes)
  if (!is.null(groups)) 
  {
    vertex.colors <- rep("background", length=nNodes)
    for (i in 1:length(groups)) vertex.colors[groups[[i]]]=color[i] 
  } else vertex.colors <- rep(color, length=nNodes)
  if (length(color)==nNodes) vertex.colors <- color
  if (all(col2rgb(background,TRUE) == col2rgb("transparent",TRUE)))
  {
    vertex.colors[vertex.colors=="background"] <- "white"
  } else  vertex.colors[vertex.colors=="background"] <- background
  
  # Label color:
  if (length(lcolor) != nNodes){
    lcolor <- rep(lcolor,nNodes)
  }
  if (any(is.na(lcolor))){
    # if (!is.null(theme) && is.character(theme) && theme == "gray"){
    #   browser()
    #   lcolor[is.na(lcolor)] <- ifelse(vertex.colors == "background",
    #                                   ifelse(mean(col2rgb(background)/255) > 0.5,"black","white"),
    #                                   ifelse(colMeans(col2rgb(vertex.colors[is.na(lcolor)])) > 0.5,"black","white")
    #   )
    # } else {
      lcolor[is.na(lcolor)] <- ifelse(vertex.colors == "background",
                                      ifelse(mean(col2rgb(background)/255) > 0.5,"black","white"),
                                      ifelse(colMeans(col2rgb(vertex.colors[is.na(lcolor)])/255) > label.color.split,"black","white")
      )
    # }
  }
  
  # Dummy groups list:
  if (is.null(groups)) 
  {
    groups <- list(1:nNodes)
  }
  
  # Scores:
  if (!is.null(scores)) 
  {
    if (length(scores)!=nNodes)
    {
      warning ("Length of scores is not equal to nuber of items")
    } else
    {
      bcolor <- vertex.colors
      if (is.null(scores.range)) scores.range=c(min(scores),max(scores))
      scores[is.na(scores)]=scores.range[1]
      rgbmatrix=1-t(col2rgb(vertex.colors)/255)
      for (i in 1:nNodes) rgbmatrix[i,]=rgbmatrix[i,] * (scores[i]-scores.range[1] ) / (scores.range[2]-scores.range[1] )
      vertex.colors=rgb(1-rgbmatrix)
    }
  }
  
  if (diagCols)
  {
    if (diagCols & !is.null(scores)) stop("Multiple modes specified for vertex colors (diag and scores)")
    if (diagCols & weighted)
    {
      if (is.null(bcolor) & !all(vertex.colors=="white")) bcolor=vertex.colors
      if (cut==0) 
      {
        colV=(abs(diagWeights)-minimum)/(maximum-minimum)
      } else 
      {
        colV=(abs(diagWeights)-minimum)/(cut-minimum)
      }
      colV[colV>1]=1
      colV[colV<0]=0
      
      if (transparency) 
      {
        vertex.colors=rep("#00000000",nNodes)
        colV=colV^(2)
        neg=col2rgb(rgb(0.75,0,0))/255
        pos=col2rgb(rgb(0,0.6,0))/255
        
        # Set colors for edges over cutoff:
        vertex.colors[diagWeights< -1* minimum] <- rgb(neg[1],neg[2],neg[3],colV[diagWeights< -1*minimum])
        vertex.colors[diagWeights> minimum] <- rgb(pos[1],pos[2],pos[3],colV[diagWeights> minimum])
      } else 
      {
        vertex.colors=rep("white",nNodes)
        vertex.colors[diagWeights>minimum]=rgb(1-colV[diagWeights > minimum],1-(colV[diagWeights> minimum]*0.25),1-colV[diagWeights > minimum])
        vertex.colors[diagWeights< -1*minimum]=rgb(1-(colV[diagWeights< (-1)*minimum]*0.25),1-colV[diagWeights < (-1)*minimum],1-colV[diagWeights < (-1)*minimum])
      }
      if (cut!=0)
      {
        # Set colors for edges over cutoff:
        vertex.colors[diagWeights<= -1*cut] <- "red"
        vertex.colors[diagWeights>= cut] <- "darkgreen"
      }
    }
  }
  if (is.null(bcolor))
  {
    bcolor <- rep(ifelse(mean(col2rgb(background)/255)>0.5,"black","white"),nNodes)
  } else {
    bcolor <- rep(bcolor,length=nNodes)
  }
  
  if (any(vTrans<255) || length(vTrans) > 1)
  {
    if ( length(vTrans) > 1 && length(vTrans) != nNodes)
    {vTrans <- 255}

    # Transparance in vertex colors:
    num2hex <- function(x)
    {
      hex=unlist(strsplit("0123456789ABCDEF",split=""))
      return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
    }
    
    colHEX <- rgb(t(col2rgb(vertex.colors)/255))
    
    vertex.colors <- paste(sapply(strsplit(colHEX,split=""),function(x)paste(x[1:7],collapse="")),num2hex(vTrans),sep="")
  }
  
  
  # Vertex size:
  if (length(vsize)==1) vsize=rep(vsize,nNodes)
  if (length(vsize2)==1) vsize2=rep(vsize2,nNodes)
  if (!edgelist) Vsums=rowSums(abs(input))+colSums(abs(input))
  if (edgelist)
  {
    Vsums=numeric(0)
    for (i in 1:nNodes) Vsums[i]=sum(c(input[,1:2])==i)
  }
  if (length(vsize)==2 & nNodes>2 & length(unique(Vsums))>1) vsize=vsize[1] + (vsize[2]-vsize[1]) * (Vsums-min(Vsums))/(max(Vsums)-min(Vsums))
  if (length(vsize)==2 & nNodes>2 & length(unique(Vsums))==1) vsize=rep(mean(vsize),nNodes)
  
  if (length(vsize2)==2 & nNodes>2 & length(unique(Vsums))>1) vsize2=vsize2[1] + (vsize2[2]-vsize2[1]) * (Vsums-min(Vsums))/(max(Vsums)-min(Vsums))
  if (length(vsize2)==2 & nNodes>2 & length(unique(Vsums))==1) vsize2=rep(mean(vsize2),nNodes)
  
  # Vertex shapes:
  if (length(shape)==1) shape=rep(shape,nNodes)
  
  # means:
  if (length(means)==1) means <- rep(means,nNodes)
  if (length(SDs)==1) SDs <- rep(SDs, nNodes)
  
  
  
  #     
  #     pch1=numeric(0)
  #     pch2=numeric(0)
  #     
  #     for (i in 1:length(shape))
  #     {
  #       if (shape[i]=="circle")
  #       {
  #         pch1[i]=16
  #         pch2[i]=1
  #       }
  #       if (shape[i]=="square")
  #       {
  #         pch1[i]=15
  #         pch2[i]=0
  #       }
  #       if (shape[i]=="triangle")
  #       {
  #         pch1[i]=17
  #         pch2[i]=2
  #       }
  #       if (shape[i]=="diamond")
  #       {
  #         pch1[i]=18
  #         pch2[i]=5
  #       }
  #       if (!shape[i]%in%c("circle","square","triangle","diamond","rectangle")) stop(paste("Shape",shape[i],"is not supported"))
  #     }
  #     
  
  
  # Arrow sizes:
  if (length(asize)==1) asize=rep(asize,length(E$from))
  
  if (length(asize)!=length(E$from)) warning("Length of 'asize' is not equal to the number of edges")
  
  
  
  # Edge labels:
  # Make labels:
  
  if (!is.logical(edge.labels))
  {
    #       edge.labels=as.character(edge.labels)
    if (length(edge.labels)!=length(E$from))
    {
      warning("Number of edge labels did not correspond to number of edges, edge labes have been ommited")
      edge.labels <- FALSE
    }
    
    if (length(edge.labels) > 0 & is.character(edge.labels))
    {
      edge.labels[edge.labels=="NA"]=""
    }  
  } else
  {
    if (edge.labels)
    {
      edge.labels= as.character(round(E$weight,2))
    } else edge.labels <- rep('',length(E$from))
  }
  
  if (is.numeric(edge.labels)) edge.labels <- as.character(edge.labels)
  
  # Bars:
  length(bars) <- nNodes
  barSide <- rep(barSide,nNodes)
  barColor <- rep(barColor, nNodes)
  barLength <- rep(barLength, nNodes)
  barColor[barColor == 'border'] <- bcolor[barColor == 'border']
  
  
  # Compute loopRotation:
  if (DoNotPlot)
  {
    loopRotation[is.na(loopRotation)] <- 0
  } 
  else 
  {
    for (i in seq_len(nNodes))
    {
      if (is.na(loopRotation[i]))
      {
        centX <- mean(layout[,1])
        centY <- mean(layout[,2])
        for (g in 1:length(groups))
        {
          if (i%in%groups[[g]] & length(groups[[g]]) > 1)
          {
            centX <- mean(layout[groups[[g]],1])
            centY <- mean(layout[groups[[g]],2])
          }
        }
        loopRotation[i] <- atan2usr2in(layout[i,1]-centX,layout[i,2]-centY)
        if (shape[i]=="square")
        {
          loopRotation[i] <- c(0,0.5*pi,pi,1.5*pi)[which.min(abs(c(0,0.5*pi,pi,1.5*pi)-loopRotation[i]%%(2*pi)))]
        }
      } 
    } 
  } 
  
  
  # Node names:
  if (is.null(nodeNames)) nodeNames <- labels
  
  
  # Make labels:
  if (is.logical(labels))
  {
    if (labels)
    {
      labels=1:nNodes
    } else 
    {
      labels <- rep('',nNodes)
    }
  }
  
  border.width <- rep(border.width, nNodes)
  
  # Node argument setup:
  borders <- rep(borders,length=nNodes)
  label.font <- rep(label.font,length=nNodes)
  
  # Make negative dashed:
  if (negDashed){
    lty[] <- ifelse(E$weight < 0, 2, 1)
  }
  
  ########### SPLIT HERE ###########
  
  ### Fill cigraph object with stuff:
  ## Edgelist:
  cigraphObject$Edgelist$from <- E$from
  cigraphObject$Edgelist$to <- E$to
  cigraphObject$Edgelist$weight <- E$weight
  cigraphObject$Edgelist$directed <- directed
  cigraphObject$Edgelist$bidirectional <- bidirectional
  
  # Nodes:
  cigraphObject$graphAttributes$Nodes$border.color <- bcolor
  cigraphObject$graphAttributes$Nodes$borders <- borders
  cigraphObject$graphAttributes$Nodes$border.width <- border.width
  cigraphObject$graphAttributes$Nodes$label.cex <- label.cex
  cigraphObject$graphAttributes$Nodes$label.font <- label.font
  cigraphObject$graphAttributes$Nodes$label.color <- lcolor
  cigraphObject$graphAttributes$Nodes$labels <- labels
  cigraphObject$graphAttributes$Nodes$names <- nodeNames
  cigraphObject$graphAttributes$Nodes$loopRotation <- loopRotation
  cigraphObject$graphAttributes$Nodes$shape <- shape
  cigraphObject$graphAttributes$Nodes$color <- vertex.colors
  cigraphObject$graphAttributes$Nodes$width <- vsize
  cigraphObject$graphAttributes$Nodes$height <- vsize2
  cigraphObject$graphAttributes$Nodes$subplots <- subplots
  cigraphObject$graphAttributes$Nodes$images <- images
  # cigraphObject$graphAttributes$Nodes$tooltips <- tooltips
  # cigraphObject$graphAttributes$Nodes$SVGtooltips <- SVGtooltips
  cigraphObject$graphAttributes$Nodes$bars <- bars
  cigraphObject$graphAttributes$Nodes$barSide <- barSide
  cigraphObject$graphAttributes$Nodes$barColor <- barColor
  cigraphObject$graphAttributes$Nodes$barLength <- barLength
  cigraphObject$graphAttributes$Nodes$means <- means
  cigraphObject$graphAttributes$Nodes$SDs <- SDs
  cigraphObject$graphAttributes$Nodes$node.label.offset <- node.label.offset
  cigraphObject$graphAttributes$Nodes$node.label.position <- node.label.position

  
  # Pies:
  cigraphObject$graphAttributes$Nodes$pieColor <- pieColor
  cigraphObject$graphAttributes$Nodes$pieColor2 <- pieColor2
  cigraphObject$graphAttributes$Nodes$pieBorder <- pieBorder
  cigraphObject$graphAttributes$Nodes$pie <- pie
  cigraphObject$graphAttributes$Nodes$pieStart <- pieStart
  cigraphObject$graphAttributes$Nodes$pieDarken <- pieDarken
  
  
  
  # Edges:
  cigraphObject$graphAttributes$Edges$curve <- curve
  cigraphObject$graphAttributes$Edges$color <- edge.color
  cigraphObject$graphAttributes$Edges$labels <- edge.labels
  cigraphObject$graphAttributes$Edges$label.cex <- edge.label.cex
  cigraphObject$graphAttributes$Edges$label.bg <- edge.label.bg
  cigraphObject$graphAttributes$Edges$label.margin <- edge.label.margin
  cigraphObject$graphAttributes$Edges$label.font <- edge.label.font
  cigraphObject$graphAttributes$Edges$label.color <- ELcolor
  cigraphObject$graphAttributes$Edges$width <- edge.width
  cigraphObject$graphAttributes$Edges$lty <- lty
  cigraphObject$graphAttributes$Edges$fade <- fade
  cigraphObject$graphAttributes$Edges$edge.label.position <- edge.label.position
  cigraphObject$graphAttributes$Edges$residEdge <- residEdge
  cigraphObject$graphAttributes$Edges$CircleEdgeEnd <- CircleEdgeEnd
  cigraphObject$graphAttributes$Edges$asize <- asize
  if (mode == "sig") cigraphObject$graphAttributes$Edges$Pvals <- Pvals else Pvals <- NULL
  cigraphObject$graphAttributes$Edges$parallelEdge <- parallelEdge
  cigraphObject$graphAttributes$Edges$parallelAngle <- parallelAngle
  cigraphObject$graphAttributes$Edges$edgeConnectPoints <- edgeConnectPoints
  
  # Knots:
  cigraphObject$graphAttributes$Knots$knots <- knots
  cigraphObject$graphAttributes$Knots$knot.size <- knot.size
  cigraphObject$graphAttributes$Knots$knot.color <- knot.color
  cigraphObject$graphAttributes$Knots$knot.borders <- knot.borders
  cigraphObject$graphAttributes$Knots$knot.border.color <- knot.border.color
  cigraphObject$graphAttributes$Knots$knot.border.width <- knot.border.width
  
  # Graph:
  cigraphObject$graphAttributes$Graph$nNodes <- nNodes
  cigraphObject$graphAttributes$Graph$weighted <- weighted
  cigraphObject$graphAttributes$Graph$edgesort <- edgesort
  cigraphObject$graphAttributes$Graph$scores <- scores
  cigraphObject$graphAttributes$Graph$scores.range <- scores.range
  cigraphObject$graphAttributes$Graph$groups <- groups
  cigraphObject$graphAttributes$Graph$minimum <- minimum
  cigraphObject$graphAttributes$Graph$maximum <- maximum
  cigraphObject$graphAttributes$Graph$cut <- cut
  cigraphObject$graphAttributes$Graph$polygonList <- polygonList
  cigraphObject$graphAttributes$Graph$mode <- mode
  cigraphObject$graphAttributes$Graph$color <- color
  
  # Layout:
  cigraphObject$layout <- layout
  cigraphObject$layout.orig <- original.layout
  
  # Plot options:
  cigraphObject$plotOptions$filetype <- filetype
  cigraphObject$plotOptions$filename <- filename
  cigraphObject$plotOptions$background <- background
  cigraphObject$plotOptions$bg <- bg
  cigraphObject$plotOptions$normalize <- normalize
  cigraphObject$plotOptions$plot <- plot
  cigraphObject$plotOptions$mar <- mar
  cigraphObject$plotOptions$GLratio <- GLratio
  cigraphObject$plotOptions$legend <- legend
  cigraphObject$plotOptions$legend.cex <- legend.cex
  cigraphObject$plotOptions$pty <- pty
  cigraphObject$plotOptions$XKCD <- XKCD
  cigraphObject$plotOptions$residuals <- residuals
  cigraphObject$plotOptions$residScale <- residScale
  cigraphObject$plotOptions$arrows <- arrows
  cigraphObject$plotOptions$arrowAngle <- arrowAngle
  cigraphObject$plotOptions$open <- open
  cigraphObject$plotOptions$curvePivot <- curvePivot
  cigraphObject$plotOptions$curveShape <- curveShape
  cigraphObject$plotOptions$curveScale <- curveScale
  cigraphObject$plotOptions$curveScaleNodeCorrection <- curveScaleNodeCorrection
  cigraphObject$plotOptions$curvePivotShape <- curvePivotShape
  cigraphObject$plotOptions$label.scale <- label.scale
  cigraphObject$plotOptions$label.scale.equal <- label.scale.equal
  cigraphObject$plotOptions$label.fill.vertical <- label.fill.vertical
  cigraphObject$plotOptions$label.fill.horizontal <- label.fill.horizontal
  cigraphObject$plotOptions$label.prop <- label.prop
  cigraphObject$plotOptions$label.norm <- label.norm
  # cigraphObject$plotOptions$overlay <- overlay
  cigraphObject$plotOptions$details <- details
  cigraphObject$plotOptions$title <- title
  cigraphObject$plotOptions$title.cex <- title.cex
  cigraphObject$plotOptions$preExpression <- preExpression
  cigraphObject$plotOptions$postExpression <- postExpression
  cigraphObject$plotOptions$legend.mode <- legend.mode
  cigraphObject$plotOptions$srt <- srt
  cigraphObject$plotOptions$gray <- gray
  # cigraphObject$plotOptions$overlaySize <- overlaySize
  cigraphObject$plotOptions$plotELBG <- plotELBG
  cigraphObject$plotOptions$alpha <- alpha
  cigraphObject$plotOptions$width <- width
  cigraphObject$plotOptions$height <- height
  cigraphObject$plotOptions$aspect <- aspect
  cigraphObject$plotOptions$rescale <- rescale
  cigraphObject$plotOptions$barsAtSide <- barsAtSide
  cigraphObject$plotOptions$bgres <- bgres
  cigraphObject$plotOptions$bgcontrol <- bgcontrol
  cigraphObject$plotOptions$resolution <- res
  cigraphObject$plotOptions$subpars <- subpars
  cigraphObject$plotOptions$subplotbg <- subplotbg
  cigraphObject$plotOptions$usePCH <- usePCH
  cigraphObject$plotOptions$node.resolution <- node.resolution
  cigraphObject$plotOptions$noPar <- noPar
  cigraphObject$plotOptions$meanRange <- meanRange
  cigraphObject$plotOptions$drawPies <- drawPies
  cigraphObject$plotOptions$pieRadius <- pieRadius
  cigraphObject$plotOptions$pastel <- pastel
  cigraphObject$plotOptions$piePastel <- piePastel
  
  cigraphObject$plotOptions$rainbowStart <- rainbowStart
  cigraphObject$plotOptions$pieCIs <- pieCIs 
  
  
  
  if (!DoNotPlot)
  {
    plot(cigraphObject)
    invisible(cigraphObject)
  } else
  {
    return(cigraphObject)
  }
  
}
