gfunc <- function(countData,colData,title) {
#### Rarefaction curve plotter ####  
  
  colData <- colData[names(countData),]
  
  # descending order each sample 
  DT <- data.table(apply(countData,2,sort,decreasing=T))
  
  # get log10(cummulative sum) of each sample
  DT <- log10(cumsum(DT))
  
  # set values larger than maximum for each column to NA
  DT <- data.table(apply(DT,2,function(x) {x[(which.max(x)+1):length(x)]<- NA;x}))
  
  # remove rows with all NA
  DT <- DT[rowSums(is.na(DT)) != ncol(DT), ]
  
  # add a count column to the data table
  DT$x <- seq(1,nrow(DT))
  
  # melt the data table for easy plotting 
  MDT <- melt(DT,id.vars="x")
  
  # create an empty ggplot object from the data table
  g <- ggplot(data=MDT,aes(x=x,y=value,colour=variable))
  
  # remove plot background and etc.
  g <- g + theme_classic_thin() %+replace% 
    theme(legend.position="none",axis.title=element_blank())
  
  # plot cumulative reads
  g <- g + geom_line(linewidth=1.5) + scale_colour_viridis(discrete=T)
  
  # add axis lables
  g <- g + ggtitle(title)
  #g <- g + ylab(expression("Log"[10]*" aligned sequenecs"))+xlab("OTU count")
  
  # print the plot
  g
}