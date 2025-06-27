getArgs <- function(args)
{
  if (length(args)>0)
  {
    iscigraph <- sapply(args,function(x)"cigraph"%in%class(x))
    argLists <- c(lapply(args[iscigraph],'[[','Arguments'),lapply(args[iscigraph],'[','layout'))
    args <- args[!iscigraph]
    newArgs <- lapply(argLists,getArgs)
    for (l in newArgs) args <- c(args,l[!names(l)%in%names(args)])
  }
  return(args)
}