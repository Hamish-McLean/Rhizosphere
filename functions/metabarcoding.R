cbPalette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", 
               "#332288", "#AA4499", "#44AA99", "#999933", 
               "#882255", "#661100", "#6699CC", "#888888")

cbPalette_small <-c("#000000", "#E69F00", "#56B4E9", "#009E73",
                    "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

em_splitter <- function(test,col="contrast",f=" - ") {
  DF <- do.call(rbind,
                strsplit(gsub(f,"\t",gsub(",","\t",test[[col]],perl=T)),"\t")
  )
  DT <- as.data.table(cbind(DF,test[,-1]))
  setnames(DT,names(DT)[1:ncol(DF)],paste0("V",names(DT)[1:ncol(DF)]))
}

triple.venn <- function(A,B,C,...) {
  ab <- sum(duplicated(c(A,B)))
  ac <- sum(duplicated(c(A,C)))
  bc <- sum(duplicated(c(B,C)))
  abc <- sum(duplicated(c(c(A,B)[duplicated(c(A,B))],C)))
  draw.triple.venn(length(A),length(B),length(C),ab,bc,ac,abc,...)
}

venn.out <- function(A,B,C) {
  ab <- c(A,B)[duplicated(c(A,B))]
  ac <- c(A,C)[duplicated(c(A,C))]
  bc <- c(B,C)[duplicated(c(B,C))]
  abc <- c(ab,C)[duplicated(c(ab,C))]
  list(AuB=ab,AuC=ac,BuC=bc,AuBuC=abc)
}



alpha_limit <- function(g,limits=c(0,4000),p=2,l=17) {
  g2 <- g + coord_cartesian(y=limits)
  g <- ggplotGrob(g)
  g2 <- ggplotGrob(g2)
  g[["grobs"]][[p]] <- g2[["grobs"]][[p]]
  g[["grobs"]][[l]] <- g2[["grobs"]][[l]]
  g
}

combineByTaxa <- 
  function (taxData,countData, rank = "species", confidence = 0.95, column_order = -8) {
    require(plyr)
    require(data.table)
    
    # reorder taxonomy table (remove rank column by default)
    taxData <- taxData[, column_order]
    
    # get rank level
    level <- which(c("kingdom", "phylum", "class", "order", 
                     "family", "genus", "species")==rank)
    
    # select taxonomy at given confidence
    taxData <- phyloTaxaTidy(taxData, confidence,level=level)
    
    # convert to data table
    TD <- data.table(taxData, keep.rownames = "OTU")
    
    # get list of OTUs in each entry at the given rank
    combinedOTUs <- ddply(TD, ~rank, summarize, OTUS = list(as.character(OTU)))
    
    # combine the OTUs into list format
    combinedOTUs <- combinedOTUs[lapply(TD[, 2], function(x) length(unlist(x))) > 
                                   1, ]
    
    list(countData = combCounts(combinedOTUs,ubiome_FUN$countData),
         taxData   = combTaxa(combinedOTUs,taxData))
    
  }



# adapted version of phyloTaxaTidy to mash names by conf
phyloTidy <- function(taxData,conf,setCount=T) {
  X <- taxData[,1:7]
  Y <- taxData[,9:15]
  Y[,1] <-1.00
  X[apply(Y,2,function(y)y<conf)] <- "unknown"
  n <- c("kingdom","phylum","class","order","family","genus","species")
  X <-apply(X,2,function(x) {
    x <- sub("\\(.*","",x)
    x <- sub("_SH[0-9]+\\..*","",x)
    x <- gsub("_"," ",x)
    #x <- sub("_"," ",x)
    x
  })
  
  td <- as.data.table(X,keep.rownames = "OTU")
  cols <- names(td)[-1]
  num <- td[,.N,by=cols]
  num2 <- invisible(apply(num,1,function(x) as.data.table(t(matrix(x[-8],nrow=7,ncol=as.numeric(x[8]))))))
  invisible(lapply(num2,setnames,cols))
  invisible(lapply(num2,function(x)x[,counts:=1:nrow(x)]))
  td2 <- do.call("rbind",num2)
  td[,(cols):=lapply(.SD,as.character),.SDcols=cols]
  setorderv(td, cols[-1])
  setorderv(td2,cols[-1])
  
  td$counts <- td2$counts
  td[,(cols):=lapply(.SD,as.factor),.SDcols=cols]
  taxData <- as.data.frame(td)
  rownames(taxData) <- taxData$OTU
  taxData <- taxData[,-1]
  if(setCount) {
    taxData$species <- as.factor(paste(taxData$species,taxData$counts,sep=" "))
    taxData <- taxData[,-8]
  }
  taxData
}

funq <- function(lfc,p){lfc[c(which(p>0.1),which(is.na(p)))] <- 0;lfc}

qf <- function(x)sub("Fungi;unknown.*","REMOVE",x)

# get the full taxonomy for each taxon_id
# get_full_taxa <- function(taxon_id,obj) {
#   t_var   <- c(obj$taxon_indexes()[[taxon_id]],obj$supertaxa()[[taxon_id]])
#   t_names <- sapply(t_var,function(x) {obj$taxa[[x]]$get_name()})
#   paste(rev(t_names),collapse = ";")
# }

get_full_taxa <- function(taxon_id,t1=t1,t2=t2) {
  t_names <- t2[c(taxon_id,unlist(t1[taxon_id]))]
  paste(rev(t_names),collapse = ";")
}

MCDESeq <- function (OBJ,formula,Rank,colData) {
  
  o <- OBJ$clone(deep=T)
  formula <- formula(formula)
  # set rank for collapsing
  #  Rank = "genus"
  # collapse OTU abundances at each taxonomic rank at and above Rank 
  o$data$tax_abund <- o %>% 
    taxa::filter_taxa(taxon_ranks == Rank,supertaxa=T) %>%
    calc_taxon_abund(data = "otu_table",cols = colData$sample_id)
  
  # get the otu table with taxon_ids  - rejig to get taxon_id as rownames
  otu_table <- as.data.table(o$data$otu_table) #metacoder:::get_taxmap_table(OBJ, "otu_table")
  
  # merge duplicate "species"
  numeric_cols <- which(sapply(otu_table, is.integer))
  otu_table <- as.data.frame(otu_table[, lapply(.SD, sum), by = taxon_id, .SDcols = numeric_cols])
  rownames(otu_table) <- otu_table[,1]
  otu_table <- otu_table[,c(-1)]
  
  # get Rank and above abundances - rejig to get taxon_id as rownames
  tax_abund <- as.data.frame(o$data$tax_abund) #metacoder:::get_taxmap_table(OBJ, "tax_abund")
  rownames(tax_abund) <- tax_abund[,1]
  tax_abund <- tax_abund[,c(-1)]
  
  # set character columns to factors
  colData <- as.data.frame(unclass(as.data.frame(colData)))
  
  # create new dds object from combined tables
  dds <- DESeqDataSetFromMatrix(rbind(otu_table,tax_abund),colData,~1)
  design(dds) <- formula
  dds <- DESeq(dds)
  dds@NAMES <- Rank
  return(dds)
}

MCDESres <- function (OBJ,contrast,alpha=0.1) {
  
  o <- OBJ$clone(deep=T)
  
  
  # set rank for collapsing
  #  Rank = "genus"
  # collapse OTU abundances at each taxonomic rank at and above Rank 
  o$data$tax_abund <- o %>% 
    taxa::filter_taxa(taxon_ranks == o$data$dds@NAMES,supertaxa=T) %>%
    calc_taxon_abund(data = "otu_table",cols = colData(o$data$dds)$sample_id)
  
  # get the otu table with taxon_ids  - rejig to get taxon_id as rownames
  otu_table <- as.data.table(o$data$otu_table) #metacoder:::get_taxmap_table(OBJ, "otu_table")
  # merge duplicate "species"
  numeric_cols <- which(sapply(otu_table, is.integer))
  otu_table <- as.data.frame(otu_table[, lapply(.SD, sum), by = taxon_id, .SDcols = numeric_cols])
  rownames(otu_table) <- otu_table[,1]
  otu_table <- otu_table[,c(-1,-2)]
  
  # get Rank and above abundances - rejig to get taxon_id as rownames
  tax_abund <- as.data.frame(o$data$tax_abund) #metacoder:::get_taxmap_table(OBJ, "tax_abund")
  rownames(tax_abund) <- tax_abund[,1]
  tax_abund <- tax_abund[,c(-1)]
  
  # set character columns to factors
  #colData <- as.data.frame(unclass(as.data.frame(colData))) # why this???
  
  # run results
  res <- results(o$data$dds,alpha=alpha,contrast=contrast)
  
  # make a data table from the resuults
  res_merge <- as.data.table(as.data.frame(res),keep.rownames="taxon_id")
  
  # add second log fold change column with fold changes for non sig taxa set to 0 (for colouring purposes on graphs)
  res_merge[,log2FoldChange2:=funq(log2FoldChange,padj)]
  
  # add full taxonomy to results (this is slow - good reason not to use S4 object model!)
  
  # order the results
  setkey(res_merge,taxon_id)
  
  # these are required due to the restrictive interace of taxa/metacoder 
  # reduces the time required to get the full taxonomy of each taxon id by about x10,000
  t1 <- o$supertaxa()
  t2 <- o$taxon_names()
  res_merge[,taxonomy:=sapply(1:nrow(res_merge),get_full_taxa,t1,t2)]
  
  # add results as diff_table (could be called anything - but for consistency with metacoder keep same names)
  #
  as_tibble(res_merge)
}


## Custom functions

gfunc <- function(countData,colData,title) {
  #### Rarefaction curve plotter ####  
  
  colData <- colData[names(countData),]
  
  # descending order each sample 
  DT <- data.table(apply(countData,2,sort,decreasing=T))
  
  # get cummulative sum of each sample
  DT <- cumsum(DT)    
  
  # log the count values                            
  DT <- log10(DT)
  
  # set values larger than maximum for each column to NA
  DT <- data.table(apply(DT,2,function(x) {
    i<-which.max(x)
    if(i<length(x)) {
      x[(i+1):length(x)]<- NA}
    x}))
  
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
  g <- g + geom_line(size=1.5) + scale_colour_viridis(discrete=T)
  
  # add axis lables
  g <- g + ggtitle(title)
  #g <- g + ylab(expression("Log"[10]*" aligned sequenecs"))+xlab("OTU count")
  
  # print the plot
  g
}

filter_otus <- function(countData,OTUFILTER=0.001){
  i <- sum(countData)
  y <- rowSums(countData)
  y<-y[order(y,decreasing = T)]
  keep <-  names(y[(cumsum(y)/i <=1-OTUFILTER)])
}

taxaCount <- function(taxData,conf=0.65,level="phylum",cutoff=1) {
  TD <- copy(taxData)
  setDT(TD)
  cols <- names(TD)[grep("_",names(TD))]
  
  TD[,(cols):=lapply(.SD,as.numeric),.SDcols=cols]
  TD[,(level):=taxaConfVec(TD[,-8],conf=conf,level=which(colnames(TD)==level))]
  TD <- TD[,.N,by=level]
  TD[,props:=N/sum(N)]
  
  txk <- TD[TD$props >= cutoff/100, ..level]
  TD_cut <- TD[txk,on=level]
  setorder(TD_cut)
  #others<-c("others",1-sum(TD[txk,on=level]$props))
  
  TD_cut <- rbind(TD_cut,
                  list("others",
                       (1-sum(TD[txk,on=level]$props))*sum(TD$N),
                       1-sum(TD[txk,on=level]$props)
                  ))
}

shift_legend <- function(p,position="center",...) {
  require(lemon)
  require(gtable)
  
  # convert plot to grob
  gpanels <- gtable::gtable_filter(ggplotGrob(p),"panel")
  
  # now we just need a simple call to reposition the legend
  lemon::reposition_legend(p, 
                           position, 
                           panel=gpanels$layout$name[sapply(gpanels$grobs,length)<4],
                           ...
  )
}

plotfun1 <- function(md1=md1,x="phylum",fill="Time.point") {
  g <- ggplot(md1,aes_string(x=x,y="value",fill=fill)) + theme_classic_thin()
  g <- g + geom_bar(stat="identity",colour="white",width = 0.8, position = position_dodge(width = 0.85))
  if(length(levels(as.factor(md1[[x]])))<=8) {
    g <- g + scale_fill_manual(values=cbPalette)
  }
  g <- g  + xlab("")+ ylab("Frequency (%)")
  g <- g + scale_y_continuous(expand = c(0, 0),labels = scaleFUN,limits=NULL)
  g <- g + guides(fill=guide_legend(title = fill))
  g <- g + theme_blank()
  g <- g + theme(legend.position = "bottom",legend.direction = "horizontal",legend.justification = "left",
                 axis.text.x = element_text(angle = 45, hjust = 1,size=14),
                 plot.margin=unit(c(0.2,0,0.1,.5),"cm"),
                 axis.line.y = element_line(colour = "black",size=1),
                 axis.ticks.x=element_blank(),
                 text=element_text(size=14),
                 axis.title.y=element_text(size=(14-2)))
  g
}

plotfun2 <- function(md1=md1,x="phylum") {
  g <- ggplot(md1,aes_string(x=x,y="value"))
  g <- g + geom_bar(stat="identity",colour="black")
  g <- g  + xlab("")+ ylab("")
  scaleFUN<-function(x) sprintf("%.0f", x)
  g <- g + scale_y_continuous(expand = c(0, 0),labels = scaleFUN,limits=NULL)
  g <- g + guides(fill=guide_legend(ncol=1))
 # g <- g+ facet_wrap(~Treatment,nrow=6) + theme_blank()
  g <- g + theme_blank()+ theme(
    strip.background = element_rect(colour="white", fill="white"),
    strip.text = element_text(size=16),
    axis.text.x = element_text(angle = 45, hjust = 1,size=14),
    plot.margin=unit(c(0.2,0,0.2,1.5),"cm"),
    panel.spacing.x = unit(2, "lines"),
    axis.line.y = element_line(colour = "black",size=1),
    axis.ticks.x=element_blank(),
    text=element_text(size=12),
    plot.title = element_text(hjust = -0.11),
    axis.title.y=element_text(size=(14-2)))
  g
}

getSummedTaxa <- function(...) {
  DT <- melt(sumTaxaAdvanced(...))
  setDT(DT)
  DT[,(design):=tstrsplit(variable," : ",fixed=TRUE)]
  DT[,variable:=NULL]
  setnames(DT,c("taxon","value",design))
}


# updated version of sumTaxaAdvanced to correct bug in confidence levels
sumTaxaAdvanced <- 
function (obj, taxon = "phylum", conf = 0.65, design = "all", 
          proportional = T, cutoff = 1, topn = 0, others = T, ordered = F) 
{
  obj$taxData <- taxa_sum <- sumTaxa(obj, taxon = taxon, design = design, 
                                     conf = conf)
  taxa_sum$taxon[grep("\\(", taxa_sum$taxon)] <- taxa_sum$taxon[sub("\\(.*", 
                                                                    " incertae sedis", taxa_sum$taxon)]
  if (!topn) {
    obj[[3]]$MLUflop <- 1
    tx <- sumTaxa(obj, taxon = taxon, "MLUflop",conf=conf)
    tx[, -1] <- prop.table(as.matrix(tx[, -1]), 2) * 100
    txk <- tx[tx[, 2] >= cutoff, 1]
  }
  else {
    taxa_sum[, ncol(taxa_sum) + 1] <- 0
    taxa_sum <- taxa_sum[order(rowSums(taxa_sum[, -1]), decreasing = T), 
    ]
    taxa_sum <- taxa_sum[, -ncol(taxa_sum)]
    txk <- taxa_sum[1:topn, 1]
  }
  if (proportional) {
    taxa_sum[, -1] <- prop.table(as.matrix(taxa_sum[, -1]), 
                                 2) * 100
  }
  taxa_cut <- taxa_sum[taxa_sum[, 1] %in% txk, ]
  taxa_cut <- taxa_cut[order(taxa_cut[, 1], decreasing = F), 
  ]
  if (others) {
    taxa_cut <- rbind(taxa_cut, setNames(data.frame(x = "others", 
                                                    t(colSums(as.data.frame(taxa_sum[!taxa_sum[, 1] %in% 
                                                                                       txk, -1])))), names(taxa_cut)))
  }
  taxa_cut <- na.omit(taxa_cut)
  taxa_cut[, 1] <- as.factor(taxa_cut[, 1])
  if (ordered) {
    taxa_cut <- taxa_cut[order(rowSums(as.data.frame(taxa_cut[, 
                                                              -1])), decreasing = T), ]
  }
  return(taxa_cut)
}

