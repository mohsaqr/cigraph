
print.cigraph <- function(x,...)
{
	out <- cbind(x$Edgelist$from,"\t",ifelse(x$Edgelist$directed,ifelse(x$Edgelist$bidir,"<->","-->"),"---"),"\t",x$Edgelist$to,"\t",round(x$Edgelist$weight,2),"\n")
	cat(c("From\t\tTo\tWeight\n"))
	apply(out,1,cat)
}

# plot.cigraph <- function(x,...) cigraph(x,...)

summary.cigraph <- function(object, ...)
{
	cat("Number of edges:\t",length(object$Edgelist$from),"\n","Number of directed edges:\t",sum(object$Edgelist$directed),"\n","Number of unique weights:\t",length(unique(object$Edgelist$weight)),"\n\n")
}