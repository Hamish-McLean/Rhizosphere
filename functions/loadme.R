theme_blank <- 
function (base_size = 16, base_family = "") 
{
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
}

theme_classic_thin <-
function (base_size = 16, base_family = "") 
{
  theme_blank(base_size = base_size, base_family = base_family) %+replace% 
    theme(axis.line.x = element_line(size = 0.3, colour = "black"), 
          axis.line.y = element_line(size = 0.3, colour = "black"), 
          axis.text = element_text(colour = "black"))
}

theme_facet_blank <- 
  function (base_size = 11, base_family = "", angle = -90, t = 2, 
            r = 0, b = 0, l = 0, hjust = 0, vjust = 0) 
  {
    theme_classic_thin(base_size = base_size, base_family = base_family) %+replace% 
      theme(panel.border = element_rect(colour = "black", fill = NA, 
                                        size = 0.5), axis.text.x = element_text(angle = angle, 
                                                                                margin = margin(t = t, r = r, b = b, l = l), hjust = hjust, 
                                                                                vjust = vjust))
  }

em_splitter <- function(test,col="contrast",f=" - ") {
  DF <- do.call(rbind,
                strsplit(gsub(f,"\t",gsub(",","\t",test[[col]],perl=T)),"\t")
  )
  DT <- as.data.table(cbind(DF,test[,-1]))
  setnames(DT,names(DT)[1:ncol(DF)],paste0("V",names(DT)[1:ncol(DF)]))
}

cbPalette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", 
               "#332288", "#AA4499", "#44AA99", "#999933", 
               "#882255", "#661100", "#6699CC", "#888888")

cbPalette_small <-c("#000000", "#E69F00", "#56B4E9", "#009E73",
                    "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
as.number <- 
  function (f, convert = as.numeric) 
  {
    convert(levels(f))[f]
  }

scaleFUN <- function(x) sprintf("%.2f", x)

check.model <- function(model){
  plot(resplot(model))
  print(summary(model))
  cat("\n")
  print(overdisp_fun(model))
}

resplot <- function(m){
  r <- residuals(m,type="pearson")
  f <- fitted(m)
  d <- data.frame(residuals=r,fitted=f)
  ggplot(d,aes(x=fitted,y=residuals))+geom_point(size=1)+geom_point(size=3.5,alpha=0.5)
  
}

filter_otus <- function(countData,OTUFILTER){
  i <- sum(countData)
  y <- rowSums(countData)
  y<-y[order(y,decreasing = T)]
  keep <-  names(y[(cumsum(y)/i <=1-OTUFILTER)])
}

# checks a glm/glmm for over dispersion
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,rdf=rdf,ratio=prat,p=pval)
}

# corrects Chisq analysis of deviance for over/under dispersed glm
qanova <- function(model,...) {
  rdf <- df.residual(model)
  
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  od_par <- Pearson.chisq/rdf
  a<-anova(model)
  a$`Pr(>Chi)` <- pchisq(a[[2]]/od_par,a[[1]],lower.tail = F)
  attributes(a)$heading <- sub("poisson",paste0("poisson (",round(od_par,3),")"),attributes(a)$heading)
  a
}

qemm <- function(model,formula,...) {
  # calculate over dispersion parameter
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  od_par <- Pearson.chisq/rdf
  
  # get the marginal means object
  emm <- summary(emmeans(model,formula,...) )
  emm$SE <- emm$SE*sqrt(od_par) 
  emm
}

qcontrast <- function(res,od_par,...) {
  # get emmeans contrasts
  contr <- summary(contrast(res,...))
  #  contr <- summary(pairs(emmeans(fm,~Treatment|Pesticide)))
  
  # multiply standard errors by od parameter
  contr$SE <- contr$SE*sqrt(od_par) #contr[,4]*sqrt(overdisp_fun(fm)[[3]])
  
  # get the z ratio 
  mean.no <- which(names(contr)=="SE")-1
  contr$z.ratio <- contr[,mean.no]/contr$SE
  
  # calculate unadjusted p values
  contr$p.value <- 2*pnorm(-abs(contr$z.ratio))
  
  # get tukey adjusted p values
  contr$p.adj <- if(is.null(attr(contr,"by.vars"))){
    if(attr(contr,"adjust")=="dunnettx"|attr(contr,"adjust")=="sidak"){
      require(MHTdiscrete)
      Sidak.p.adjust(contr$p.value)
    }
    else {
     ptukey(abs(contr$z.ratio)*sqrt(2),nmeans=nrow(contr),df=Inf,lower=F)
     # ptukey(abs(contr$z.ratio),nmeans=nrow(summary(res)),df=Inf,lower=F)
    }
  } else {
    by.vars <- attr(contr,"by.vars")
    unlist(lapply(abs(reshape2::dcast(contr,contrast~get(by.vars),value.var = "z.ratio")[,-1])*sqrt(2),ptukey,nmeans=nrow(summary(res)),df=Inf,lower=F))
  }
  contr
}

qpairs <- function(res,od_par,...) {
  
  # get emmeans contrasts
  contr <- summary(pairs(res,...))
  #  contr <- summary(pairs(emmeans(fm,~Treatment|Pesticide)))
  
  # multiply standard errors by od parameter
  contr$SE <- contr$SE*sqrt(od_par) #contr[,4]*sqrt(overdisp_fun(fm)[[3]])
  
  # get the z ratio 
  contr$z.ratio <- contr[[(ncol(contr)-4)]]/contr$SE
  
  # calculate unadjusted p values
  contr$p.value <- 2*pnorm(-abs(contr$z.ratio))
  
  # get tukey adjusted p values
  contr$p.adj <- if(!is.null(attr(contr,"by.vars"))){
   # ptukey(abs(contr$z.ratio),nmeans=nrow(summary(res)),df=Inf,lower=F)
    ptukey(abs(contr$z.ratio)*sqrt(2),nmeans=nrow(summary(res)),df=Inf,lower=F)
  } else {
    p.adjust(contr$p.value,method="BH")
    #by.vars <- attr(contr,"by.vars")
    #unlist(lapply(abs(reshape2::dcast(contr,contrast~get(by.vars),value.var = "z.ratio")[,-1]),ptukey,nmeans=2,df=Inf,lower=F))
  }
  contr
}

#z.vals <- abs(reshape2::dcast(contr,contrast~get(by.vars),value.var = "z.ratio")[,-1])
#ptukey(z.vals,length(z.vals),df=Inf,lower =F)

# deviance_ratio <- function(model,specs){
# 
#   #pair labels as contrasts                      
#   contrasts <- combn(model$xlevels[[specs]],2)
#   
#   ##add self pairs
#   self=rbind(model$xlevels[[specs]],model$xlevels[[specs]])
#   
#   contrasts=cbind(contrasts,self)
#   
# }


t_func <- function (object, ..., model.names = NULL) 
{
  mCall <- match.call(expand.dots = TRUE)
  dots <- list(...)
  .sapply <- function(L, FUN, ...) unlist(lapply(L, FUN, ...))
  modp <- as.logical(vapply(dots, FUN = is, "glmmTMB", 
                            FUN.VALUE = NA))
  if (any(modp)) {
    mods <- c(list(object), dots[modp])
    nobs.vec <- vapply(mods, nobs, 1L)
    vapply(mods[-1], glmmTMB:::CompareFixef, mod1 = mods[[1]], FUN.VALUE = TRUE)
    if (var(nobs.vec) > 0) 
      stop("models were not all fitted to the same size of dataset")
    if (is.null(mNms <- model.names)) 
      mNms <- vapply(as.list(mCall)[c(FALSE, TRUE, modp)], 
                     deparse1, "")
    if (any(duplicated(mNms))) {
      warning("failed to find unique model names, assigning generic names")
      mNms <- paste0("MODEL", seq_along(mNms))
    }
    if (length(mNms) != length(mods)) 
      stop("model names vector and model list have different lengths")
    names(mods) <- sub("@env$", "", mNms)
    llks <- lapply(mods, logLik)
    ii <- order(Df <- vapply(llks, attr, FUN.VALUE = numeric(1), 
                             "df"))
    mods <- mods[ii]
    llks <- llks[ii]
    Df <- Df[ii]
    calls <- lapply(mods, getCall)
    data <- lapply(calls, `[[`, "data")
    if (!all(vapply(data, identical, NA, data[[1]]))) 
      #  stop("all models must be fit to the same data object")
      header <- paste("Data:", glmmTMB:::abbrDeparse(data[[1]]))
    subset <- lapply(calls, `[[`, "subset")
    if (!all(vapply(subset, identical, NA, subset[[1]]))) 
      #stop("all models must use the same subset")
      if (!is.null(subset[[1]])) 
        header <- c(header, paste("Subset:", glmmTMB:::abbrDeparse(subset[[1]])))
    llk <- unlist(llks)
    chisq <- 2 * pmax(0, c(NA, diff(llk)))
    dfChisq <- c(NA, diff(Df))
    val <- data.frame(Df = Df, AIC = .sapply(llks, AIC), 
                      BIC = .sapply(llks, BIC), logLik = llk, deviance = -2 * 
                        llk, Chisq = chisq, `Chi Df` = dfChisq, 
                      `Pr(>Chisq)` = pchisq(chisq, dfChisq, lower.tail = FALSE), 
                      row.names = names(mods), check.names = FALSE)
    class(val) <- c("anova", class(val))
    forms <- lapply(lapply(calls, `[[`, "formula"), 
                    deparse)
    ziforms <- lapply(lapply(calls, `[[`, "ziformula"), 
                      deparse)
    dispforms <- lapply(lapply(calls, `[[`, "dispformula"), 
                        deparse)
    # structure(val, heading = c(header, "Models:", paste(paste(paste(rep(names(mods), 
    #                                                                      times = lengths(forms)), unlist(forms), sep = ": "), 
    #                                                            unlist(ziforms), sep = ", zi="), unlist(dispforms), 
    #                                                     sep = ", disp=")))
    val[,c(1,5,7,6,8)]
  }
  else stop("no single-model anova() method for glmmTMB")
}


deviance_ratio <-  function(model,fact,type="pairwise",ref=1,adjust="BH",print=T,...) {
  
  FUN <- function(model,fact,type="pairwise",ref=1) {
    
    modelClass <- class(model)
    
    frame <- if("glm"%in%modelClass) {
      "model"
    } else if ("glmmTMB"%in%modelClass){
      "frame"
    } else {
      stop("\nUnsupported model")
    }
    
    factorLevels <- levels(model[[frame]][[fact]])
    
    #contrasts <- combn(model$xlevels[[fact]],2)
    
    X<-do.call(rbind,lapply((1:(length(factorLevels)-1)),function(i) {
      do.call(rbind,lapply((i+1):length(factorLevels),function(j) {
        x<-1:length(factorLevels)
        x[[i]] <- i
        x[[j]] <- i
        c(factorLevels[x],factorLevels[i],factorLevels[j])
      }))
    }))
    
    contrast <- paste(X[,(ncol(X)-1)],X[,ncol(X)],sep="-")
    
    X <- X[,-(ncol(X)-1):-ncol(X)]
    
    ratios <- apply(as.data.frame(X),1,function(v){
      m <- model
      m[[frame]][[fact]] <- factor(m[[frame]][[fact]],
                                   levels=levels(m[[frame]][[fact]]),
                                   labels=v)
      names(m[[frame]]) <-  gsub("(.*\\()([a-zA-Z0-9]*)(\\).*)","\\2",names(m[[frame]]))
      
      m <- update(m,.~.,data=m[[frame]])
      a <- if("glmmTMB"%in%modelClass)t_func(m,model)
      else(anova(m,model))
      
      dev <- a[2,4] 
      df <- a[2,3]
      
      rdf <- df.residual(model)
      rp <- residuals(model,type="pearson")
      Pearson.chisq <- sum(rp^2)
      od_par <- Pearson.chisq/rdf
      
      r <- pchisq(dev/od_par,df,lower.tail = F)
      c(dev,df,r)
    }) 
    
    colnames(ratios) <- contrast
    rownames(ratios) <- c("Dev","Df","P.value")
    
    
    t(ratios)
  }
  X <- FUN(model,fact,type,ref) 
  X <- cbind(X,P.adj=p.adjust(X[,"P.value"], method=adjust,...))
  if(print) {
    stats::printCoefmat(
      X,signif.stars = T,signif.legend = T,
      tst.ind = 1,P.values=T,has.Pvalue = T,
      zap.ind=F,cs.ind=3
    )
  }
  invisible(X)
  # class(X) <- "lm"
  #X
}




shift_legend <- function(p,position="center",...) {
  require(lemon)
  
  # convert plot to grob
  gpanels <- gtable::gtable_filter(ggplotGrob(p),"panel")
  
  # now we just need a simple call to reposition the legend
  lemon::reposition_legend(p, 
                           position, 
                           panel=gpanels$layout$name[sapply(gpanels$grobs,length)<4],
                           ...
  )
}


rplot <- function(m) {
  r <- resid(m,type="pearson")
  yh <- fitted(m)
  show.r <- "Resid vs fitted"
  d <- tibble(predicted = yh, resid = r, label = show.r)
  g <- ggplot(d, aes(x = predicted, y = resid, label = label))
  g <- g + labs(x = "Fitted Value", y = "Residual", title = "Residual vs Fitted Values")
  g <- g + geom_point(shape = 1, colour = "blue")
  g <- g + geom_hline(yintercept = 0, colour = "red")
  #g <- g + geom_label(na.rm = T)
  g
}

logit     <- function(x) { log(x/(1-x)) }
inv.logit <- function(x) {exp(x)/(1+exp(x))}

angular <- function() {
  ## link
  linkfun <- function(y) asin(sqrt(y))
  ## inverse link
  linkinv <- function(eta) sin(eta)^2
  mu.eta <- function(eta) { 2*sin(eta)}
  valideta <- function(eta) T
  link <- "asin.sqrt"
  structure(list(linkfun=linkfun,linkinv=linkinv,mu.eta=mu.eta,valideta=valideta,name=link),class="link-glm")
}

string_to_seed <- function(s) sum(utf8ToInt(s))


# ang <- angular()
# 
# tran <- make.tran("asin.sqrt", 100)
# 
# invAngular <- function(x){
#   x <- sin(x)^2
#   x[x>1] <- 1
#   x[x<0] <- 0
#   x*100
# }
# 
# tran$linkinv <- invAngular
# 
# tran_list <- list(list(linkfun <-function(y)y,
#                        linkinv <-function(eta)eta,
#                        mu.eta  <-function(eta)eta,
#                        valideta<-function(eta)TRUE,
#                        name="untransformed"),
#                   tran)


SE_fun <- function(model,parNums){
  
  covmat <- vcov(model)[parNums,parNums]
  sqrt(c(1,-1)%*%covmat%*%c(1,-1))
  
}

# get standard error
standard_error <- function(x) sd(x,na.rm=T) / sqrt(length(x[!is.na(x)])) 


# Emmeans pairs output in distance matrix format
pairs_dist <- function(res,div=" - ",over=NULL,...){
  X<-summary(pairs(res))
  setDT(X)
  X[,c("row","col"):=tstrsplit(contrast,div)]
  Y<-X[,c("p.value","row","col",over),with=F]
  setDF(Y)
  Y<-rbind(Y,list(NA,Y[nrow(Y),3],NA))
  Y<-rbind(Y,list(NA,NA,Y[1,2]))

  Y<-spread(Y,row, p.value, fill=1)
  Y<- Y[-nrow(Y),]
  Y <- Y[,-ncol(Y)]
#  Y <- Y[,-ncol(Y)]
  #Y[nrow(Y),] <- list(names(Y)[ncol(Y)],)

  Y %>%
    column_to_rownames(var="col") %>%
    as.dist
}


pairs_dist <- function(res,div=" - ",over=NULL,p.value="p.value",...){
  X<-if(class(res)[[1]]=="emmGrid"){
    summary(pairs(res))
    } else res
  setDT(X)
  X[,c("row","col"):=tstrsplit(contrast,div)]
  if(!is.null(over)){
    X<-split(X,by=over)
    } else X <- list(X)
  Y<- lapply(X,function(X) X[,c(p.value,"row","col",over),with=F])
  lapply(Y,setDF)  
  lapply(Y,function(Y) {
    Y<-Y[,1:3]
    ord<-unique(c(Y$row,Y$col))
    Y<-rbind(Y,list(NA,Y[nrow(Y),3],NA))
    Y<-rbind(Y,list(NA,NA,Y[1,2]))
    Y<-spread(Y,row, p.value, fill=1)
    Y<- Y[-nrow(Y),]
    Y <- Y[,-ncol(Y)]
    Y<-Y[match(ord,Y$col),match(ord,names(Y))]
    
    #Y <- Y[,-ncol(Y)]
    #Y[nrow(Y),] <- list(names(Y)[ncol(Y)],)
    rownames(Y)<-ord
   # Y<-Y[,-1]
    as.dist(Y)
    # 
 #    Y %>% 
    #  # pivot_wider(names_from=row,values_from=p.value) %>% 
  #     column_to_rownames(var="col") %>%
    #   as.dist
  })
}




plot_lm <-
  
  function(
    x,
    plots=c(1L:4L),
    id.n=3,
    labels.id = names(residuals(x)),
    caption=NULL,
    arrange=T,
    ...
  ){
    suppressPackageStartupMessages(require(tidyverse))
    suppressPackageStartupMessages(require(gridExtra))
    dropInf <- function(x, h) {
      if (any(isInf <- h >= 1)) {
        warning(gettextf("not plotting observations with leverage one:\\n  %s",
                         paste(which(isInf), collapse = ", ")), call. = FALSE,
                domain = NA)
        x[isInf] <- NaN
      }
      x
    }
    if (!inherits(x, "lm")) stop("use only with '(g)lm' objects")
    show <- rep(F,6)
    show[plots] <- T
    r <- residuals(x)
    yh <- predict(x)
    n <- length(r)
    s <-  sqrt(deviance(x)/df.residual(x))
    hii <- (infl <- influence(x, do.coef = FALSE))$hat
    ylab23 <- "Standardized residuals"
    r.w <- r
    rs <- dropInf(r.w/(s * sqrt(1 - hii)), hii)
    
    l.fit <- "Fitted values"
    
    if (is.null(id.n))
      id.n <- 0
    else {
      id.n <- as.integer(id.n)
      if (id.n < 0L || id.n > n)
        stop(gettextf("'id.n' must be in {1,..,%d}", n),
             domain = NA)
    }
    if (id.n > 0L) {
      if (is.null(labels.id))
        labels.id <- paste(1L:n)
      show.r <- show.rs <- labels.id
      show.r[sort.list(abs(r), decreasing = TRUE)[(id.n+1):length(show.rs)]] <- NA
      show.rs[sort.list(abs(rs), decreasing = TRUE)[(id.n+1):length(show.rs)]] <- NA
    }
    
    gl <- list()
    
    # plot histogram
    if (show[1]) {
      h <- hist(r, plot = FALSE)
      xfit <- seq(min(r),max(r),length=80)
      yfit <- dnorm(xfit, mean = mean(r), sd = sd(r)) * diff(h$mids[1:2]) * n
      
      d2 <- tibble(x = xfit, y = yfit)
      d <- tibble(x = r)
      g <- ggplot(d, aes(x = x))
      g <- g + labs(x="Residuals",title="Residual Histogram")
      g <- g + geom_histogram(bins = 6, color = "black", fill = "#ADD8E6")
      g <- g + geom_line(data = d2, aes(x = x, y = y), color = "#0000A0", size = 1.2)
      gl$histogram <- g
    }
    # end plot
    
    # plot normal q-q
    if (show[2]) {
      y_q <- quantile(r[!is.na(r)], c(0.25, 0.75))
      x_q <- qnorm(c(0.25, 0.75))
      slope <- diff(y_q)/diff(x_q)
      int <- y_q[1L] - slope * x_q[1L]
      d <- merge(ggplot_build(ggplot(tibble(x = r), aes(sample = x)) + stat_qq())$data[[1]],tibble(y = r,label=show.rs))
      g <- ggplot(d, aes(x = x, y=y,label=label))
      g <- g + labs(x= "Theoretical Quantiles", y="Sample Quantiles", title = "Normal Q-Q Plot")
      g <- g + geom_point(colour="blue")
      g <- g + geom_abline(slope = slope, intercept = int, color = "red")
      g <- g + geom_label(na.rm=T)
      gl$norm_qq <- g
    }
    # end plot
    
    # plot residual fit
    if(show[3]) {
      d <- tibble(predicted=yh,resid=r,label=show.r)
      g <- ggplot(d, aes(x = predicted, y = resid,label=label))
      g <- g + labs(x="Fitted Value", y="Residual", title="Residual vs Fitted Values")
      g <- g + geom_point(shape = 1,colour = "blue")
      g <- g + geom_hline(yintercept = 0,colour = "red")
      g <- g + geom_label(na.rm=T)
      gl$resid_fit <- g
    }
    # end plot
    
    # plot Scale-Location
    if(show[4]) {
      sqrtabsr <- sqrt(abs(rs))
      ylim <- c(0, max(sqrtabsr, na.rm = TRUE))
      yl <- as.expression(substitute(sqrt(abs(YL)), list(YL = as.name(ylab23))))
      yhn0 <-  yh
      d <- tibble(x=yhn0,y=sqrtabsr,label=show.rs)
      g <- ggplot(d, aes(x,y,label=label))
      g <- g + labs(x=l.fit,y=yl,title="Scale-Location Plot")
      g <- g + geom_point() + geom_smooth(...)
      g <- g + geom_label(na.rm=T)
      gl$scale_loc <- g
    }
    # end plot
    
    # plot density
    if(show[5]) {
      d <- tibble(x = r)
      g <- ggplot(d, aes(x = x))
      g <- g + labs(x="Residuals",title="Residual Density Plot")
      g <- g + geom_density()
      gl$density <- g
    }
    if(arrange) {
      return(grid.arrange(grobs=gl,se=F))
    } else {
      return(gl)
    }
  }


read.delim <- function(...)data.table::data.table(utils::read.delim(...))

cb <- function(){
  X<-list(dat=data.table::data.table(utils::read.delim("clipboard")))
  invisible(mapply(assign,names(X),X,MoreArgs=list(envir = globalenv())))
}

# henderson tilton analysis - this is a fudge
# a better method is reqiured for zeros (i.e. not ht analysis)
henderson_tilton <- function(control_pre,control_x,treat_pre,treat_x) {

  i <- control_pre/control_x
  j <- treat_x/treat_pre
  i[abs(i)==Inf] <- NA
  j[abs(j)==Inf] <- NA
  x<-(1 - (i)*(j))
  #x[is.na(x)]<-0
  x
}

# Constants

months <- data.frame(month=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),
                     postition = 1:12)

seasons <- data.frame(month=months$month,season=c(rep("Winter",2),rep("Spring",3),rep("Summer",3),rep("Autumn",3),"Winter"))

# get_grouping <- function(a,b,p) {
#   library(data.table)
#   library(combinat)
#   #imports
#   #from itertools import combinations
# 
#   # Remove elements from a list(x) if the element contains both a and b
#   remove_levels <-function(a,b,X){
#   
#     Y<-lapply(X, function(x){
#       as.data.table(lapply(x,function(i){
#         return(i[!(a%in%i&b%in%i)])
#       })) 
#     
#     })
#     Y <- lapply(Y,function(dt){
#       dt[,which(unlist(lapply(dt, function(x)!all(is.na(x))))),with=F]
#     })
#     
#     Y[sapply(Y, function(x) dim(x)[1]) > 0]
#     
# # for i in range(len(x) - 1, -1, -1) # loop over x backwards
#   #   if(a in x[i] and b in x[i]):
#   #   # print(x[i])
#   #   x.remove(x[i])
#   }
# 
#   # Generate all combinations from input list
#   get_combinations <- function(r) {
#     as.list(lapply(2:length(r), function(i){
#       as.data.table(combn(r,i)  )
#     }))
#   }
# 
#   # for i in range(2,len(r)+1):
#   #     x.append(list(combinations(r,i)))
#   #   #for j in combinations(r,i):
#   #   #   x.append(j)
#   #   return(x)
#   # }
# 
#   #X = itertools.combinations(range(4),2)
# 
#   x=get_combinations(1:5)
# 
#   x <- remove_levels(1,2,x)
#   x <- remove_levels(1,3,x)
#   x<- remove_levels(1,4,x)
#   x <- remove_levels(2,4,x)
# 
#   y=data.table()
#   while( x) {
#     y.append(x.pop())
#     for (i in get_combinations(y[-1])) {
#       if(i in x) {
#         x.remove(i)
#       }
#     }
#   }
# 
# }