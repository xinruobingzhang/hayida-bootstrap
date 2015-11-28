library(rms)
library(MASS)
##全模型部分自动二分类数据
rm( list = ls ( all = TRUE))
#求bootstrap的三个指标AUC,R2,Slope
calcum.boot<-function(fit){
  stats <- fit$stats
  Dxy<-stats["Dxy"]
  shrink <- 1  #此时将Slope的值默认为1
  R2 <- stats["R2"] #R2的值
  AUC<-(1+Dxy)/2 #AUC的值
  z<-c(AUC,R2,shrink)
  names(z)<-c("AUC","R2","Slope")
  z
}

#将list列表组装成form公式
pack.form<-function(listX){
  #  xnam<-paste0("X",listX)
  form<-as.formula(paste("Y~",paste(listX,collapse="+")))
  form
}

#将form公式分离成list表
division.form<-function(form){
  team<-terms(form)
  chrlist<-attr(team,"term.labels")
  chrlist
}

#单变量分析函数
univar.analy<-function(form,pvalue=0.1,mydata){
  listchr<-division.form(form)
  nlen<-length(listchr)
  lisre<-c()
  for(i in 1:nlen){
    X<-listchr[i]
    mylogit <- lrm(as.formula(paste("Y~",X)), x=TRUE, y=TRUE, data = mydata)
    stats<-mylogit$stats#得到一系列值
    chip<-stats["P"]
    if(chip<=pvalue){ #判断P值是否小于给定的阈值
      lisre<-append(lisre,listchr[i])#注意其返回值
      #      print(X)
    }
  }
  lisre 
}

#抽样函数
boot.sample<-function(mydata){
  #  index<-list(1:264)
  #  data1<-t(mydata)
  #  data2<-as.data.frame(data1)
  #  list<-sample(data2,replace=T)
  
  #  data3<-t(list)
  #  as.data.frame(data3)
  mydata[sample(nrow(mydata),replace=T),] #有放回的进行抽样，得到抽样后的数据
  # write.table(sample(nrow(mydata),replace=T),file='suiji.csv')
}

#公式计算slope的值，主要利用rms包里面的函数进行计算的。
val.probnew <- function(p, y, logit, group, weights=rep(1,length(y)),
                        normwt=FALSE, smooth=TRUE, logistic.cal=TRUE,
                        lim=c(0,1), m, g, cuts, emax.lim=c(0,1),
                        legendloc=lim[1] + c(.55 * diff(lim), .27 * diff(lim)),
                        statloc=c(0,.99), riskdist="calibrated", cex=.7, mkh=.02,
                        connect.group=FALSE, connect.smooth=TRUE, 
                        g.group=4, evaluate=100, nmin=0)
{
  
  if(missing(p)) p <- plogis(logit)	else logit <- qlogis(p)
  if(length(p) != length(y)) stop("lengths of p or logit and y do not agree")
  names(p) <- names(y) <- names(logit) <- NULL
  
  Spi <- function(p, y) {
    z <- sum((y - p)*(1 - 2*p)) /
      sqrt(sum((1 - 2 * p) * (1 - 2 * p) * p * (1-p)))
    P <- 2 * (1 - pnorm(abs(z)))
    c(Z=z, P=P)
  }
  
  if(!missing(group)) {
    if(length(group)==1 && is.logical(group) && group)
      group <- rep('', length(y))
    if(!is.factor(group)) group <- 
        if(is.logical(group) || is.character(group)) 
          as.factor(group) else cut2(group, g=g.group)
    names(group) <- NULL
    nma <- !(is.na(p + y + weights) | is.na(group))
    ng  <- length(levels(group))
  }
  else {
    nma <- !is.na(p + y + weights)
    ng <- 0
  }
  
  logit <- logit[nma]
  y     <- y[nma]
  p     <- p[nma]
  if(ng > 0) {
    group   <- group[nma]
    weights <- weights[nma]
    return(val.probg(p, y, group, evaluate, weights, normwt, nmin) )
  }
  
  if(length(unique(p)) == 1) {
    P     <- mean(y)
    Intc  <- qlogis(P)
    n     <- length(y)
    D     <- -1 / n
    L01   <- -2 * sum(y * logit - logb(1 + exp(logit)), na.rm=TRUE)
    L.cal <- -2 * sum(y * Intc - logb(1 + exp(Intc)), na.rm=TRUE)
    U.chisq <- L01 - L.cal
    U.p <- 1 - pchisq(U.chisq, 1)
    U <- (U.chisq - 1) / n
    Q <- D - U
    spi <- unname(Spi(p, y))
    stats <- c(0, .5, 0, D, 0, 1, U, U.chisq, U.p, Q,
               mean((y - p[1]) ^ 2), Intc, 0,
               rep(abs(p[1] - P), 2), spi)
    names(stats) <- c("Dxy","C (ROC)", 
                      "R2","D","D:Chi-sq","D:p","U","U:Chi-sq","U:p","Q",
                      "Brier","Intercept","Slope","Emax","Eavg",
                      "S:z", "S:p")
    return(stats)
  }
  
  i <- ! is.infinite(logit)
  nm <- sum(!i)
  if(nm > 0) warning(paste(nm,
                           "observations deleted from logistic calibration due to probs. of 0 or 1"))
  f.fixed <- lrm.fit(logit[i], y[i], initial=c(0., 1.), maxit=1L)
  f.recal <- lrm.fit(logit[i], y[i])
  stats <- f.fixed$stats
  n <- stats["Obs"]
  predprob <- seq(emax.lim[1], emax.lim[2], by=.0005)
  lt <- f.recal$coef[1] + f.recal$coef[2] * qlogis(predprob)
  calp <- plogis(lt)
  emax <- max(abs(predprob - calp))
  
  Sm <- lowess(p, y, iter=0)
  cal.smooth <- approx(Sm, xout=p, ties=mean)$y
  eavg <- mean(abs(p - cal.smooth))
  lr <- stats["Model L.R."]
  p.lr <- stats["P"]
  D <- (lr - 1) / n
  L01 <- -2 * sum(y * logit - logb(1 + exp(logit)), na.rm=TRUE)
  U.chisq <- L01 - f.fixed$deviance[2]
  p.U <- 1 - pchisq(U.chisq, 2)
  U <- (U.chisq - 2)/n
  Q <- D - U
  Dxy <- stats["Dxy"]
  C   <- stats["C"]
  R2  <- stats["R2"]
  B   <- mean((p - y) ^ 2)
  spi <- unname(Spi(p, y))
  stats <- c(Dxy, C, R2, D, lr, p.lr, U, U.chisq, p.U, Q, B, f.recal$coef,
             emax, spi)
  names(stats) <- c("Dxy","C (ROC)", 
                    "R2","D","D:Chi-sq","D:p","U","U:Chi-sq","U:p","Q",
                    "Brier","Intercept","Slope","Emax","S:z","S:p")
  stats <- c(stats, c(Eavg=eavg))
  stats#返回很多指标，可以取自己想要的指标
}

#利用上面公式得到slope的值
calcum.slope<-function(fit,mydata){
  p=predict(fit,mydata) #模型的预测值
  pred=exp(p)/(1+exp(p))#模型的预测概率
  slo<-val.probnew(p=pred,mydata$Y,logistic.cal = FALSE)
  slo["Slope"] #得到slope的值
}

#计算test中的三个指标的值
calcum.test<-function(fit,mydata,form){
  team<-terms(form)
  chrlist<-attr(team,"term.labels")
  x<-as.matrix(mydata[,chrlist])
  #  x<-as.matrix(fit$x)
  coee=fit$coefficients
  k<-fit$non.slopes
  null.model <- length(fit$coefficients)==k
  penalty.matrix=fit$penalty.matrix
  penalty.matrix
  refit <- if(null.model) lrm.fit(y=y) else lrm.fit( x=x, y=mydata$Y,penalty.matrix=penalty.matrix,initial=coee,maxit=1)
  kr <- refit$non.slopes
  stats <- refit$stats
  lr <- stats["Model L.R."]
  Dxy <- stats["Dxy"]
  shrink <- if(null.model) 1 else refit$coefficients[kr + 1]
  slope<-calcum.slope(fit,mydata)
  R2 <- stats["R2"]
  AUC<-(1+Dxy)/2 
  z<-c(AUC,R2,slope) #得到AUC R2 Slope三个因素的值
  names(z)<-c("AUC","R2","Slope")
  z
}

#calcum apparent three factor
apparent.factor<-function(form,mydata){
  fit<-lrm(form,x=TRUE, y=TRUE, data = mydata)
  calcum.boot(fit)
}
stepAICC <-
  function(object, scope, scale = 0,
           direction = c("both", "backward", "forward"),
           trace = 1, keep = NULL, steps = 1000, use.start = FALSE, k = 2, ...)
  {
    mydeviance <- function(x, ...)
    {
      dev <- deviance(x)
      if(!is.null(dev)) dev else extractAIC(x, k=0)[2L]
    }
    
    cut.string <- function(string)
    {
      if(length(string) > 1L)
        string[-1L] <- paste("\n", string[-1L], sep = "")
      string
    }
    
    re.arrange <- function(keep)
    {
      namr <- names(k1 <- keep[[1L]])
      namc <- names(keep)
      nc <- length(keep)
      nr <- length(k1)
      array(unlist(keep, recursive = FALSE), c(nr, nc), list(namr, namc))
    }
    
    step.results <- function(models, fit, object, usingCp=FALSE)
    {
      change <- sapply(models, "[[", "change")
      rd <- sapply(models, "[[", "deviance")
      dd <- c(NA, abs(diff(rd)))
      rdf <- sapply(models, "[[", "df.resid")
      ddf <- c(NA, abs(diff(rdf)))
      AIC <- sapply(models, "[[", "AIC")
      heading <- c("Stepwise Model Path \nAnalysis of Deviance Table",
                   "\nInitial Model:", deparse(formula(object)),
                   "\nFinal Model:", deparse(formula(fit)),
                   "\n")
      aod <-
        if(usingCp)
          data.frame(Step = change, Df = ddf, Deviance = dd,
                     "Resid. Df" = rdf, "Resid. Dev" = rd,
                     Cp = AIC, check.names = FALSE)
      else data.frame(Step = change, Df = ddf, Deviance = dd,
                      "Resid. Df" = rdf, "Resid. Dev" = rd,
                      AIC = AIC, check.names = FALSE)
      attr(aod, "heading") <- heading
      class(aod) <- c("Anova", "data.frame")
      fit$anova <- aod
      fit
      #      print(aod)
    }
    
    Terms <- terms(object)
    object$formula <- Terms
    if(inherits(object, "lme")) object$call$fixed <- Terms
    else if(inherits(object, "gls")) object$call$model <- Terms
    else object$call$formula <- Terms
    if(use.start) warning("'use.start' cannot be used with R's version of 'glm'")
    md <- missing(direction)
    direction <- match.arg(direction)
    backward <- direction == "both" | direction == "backward"
    forward <- direction == "both" | direction == "forward"
    if(missing(scope)) {
      fdrop <- numeric()
      fadd <- attr(Terms, "factors")
      if(md) forward <- FALSE
    } else {
      if(is.list(scope)) {
        fdrop <- if(!is.null(fdrop <- scope$lower))
          attr(terms(update.formula(object, fdrop)), "factors")
        else numeric()
        fadd <- if(!is.null(fadd <- scope$upper))
          attr(terms(update.formula(object, fadd)), "factors")
      } else {
        fadd <- if(!is.null(fadd <- scope))
          attr(terms(update.formula(object, scope)), "factors")
        fdrop <- numeric()
      }
    }
    models <- vector("list", steps)
    if(!is.null(keep)) keep.list <- vector("list", steps)
    n <- nobs(object, use.fallback = TRUE)  # might be NA
    fit <- object
    bAIC <- extractAIC(fit, scale, k = k, ...)
    edf <- bAIC[1L]
    bAIC <- bAIC[2L]
    if(is.na(bAIC))
      stop("AIC is not defined for this model, so 'stepAIC' cannot proceed")
    if(bAIC == -Inf)
      stop("AIC is -infinity for this model, so 'stepAIC' cannot proceed")
    nm <- 1
    Terms <- terms(fit)
    if(trace) {
      cat("Start:  AIC=", format(round(bAIC, 2)), "\n",
          cut.string(deparse(formula(fit))), "\n\n", sep='')
      utils::flush.console()
    }
    models[[nm]] <- list(deviance = mydeviance(fit), df.resid = n - edf,
                         change = "", AIC = bAIC)
    if(!is.null(keep)) keep.list[[nm]] <- keep(fit, bAIC)
    usingCp <- FALSE
    while(steps > 0) {
      steps <- steps - 1
      AIC <- bAIC
      ffac <- attr(Terms, "factors")
      ## don't drop strata terms
      if(!is.null(sp <- attr(Terms, "specials")) &&
         !is.null(st <- sp$strata)) ffac <- ffac[-st,]
      scope <- factor.scope(ffac, list(add = fadd, drop = fdrop))
      aod <- NULL
      change <- NULL
      if(backward && length(scope$drop)) {
        aod <- dropterm(fit, scope$drop, scale = scale,
                        trace = max(0, trace - 1), k = k, ...)
        rn <- row.names(aod)
        row.names(aod) <- c(rn[1L], paste("-", rn[-1L], sep=" "))
        ## drop all zero df terms first.
        if(any(aod$Df == 0, na.rm=TRUE)) {
          zdf <- aod$Df == 0 & !is.na(aod$Df)
          nc <- match(c("Cp", "AIC"), names(aod))
          nc <- nc[!is.na(nc)][1L]
          ch <- abs(aod[zdf, nc] - aod[1, nc]) > 0.01
          if(any(is.finite(ch) & ch)) {
            warning("0 df terms are changing AIC")
            zdf <- zdf[!ch]
          }
          ## drop zero df terms first: one at time since they
          ## may mask each other
          if(length(zdf) > 0L)
            change <- rev(rownames(aod)[zdf])[1L]
        }
      }
      if(is.null(change)) {
        if(forward && length(scope$add)) {
          aodf <- addterm(fit, scope$add, scale = scale,
                          trace = max(0, trace - 1), k = k, ...)
          rn <- row.names(aodf)
          row.names(aodf) <- c(rn[1L], paste("+", rn[-1L], sep=" "))
          aod <-
            if(is.null(aod)) aodf
          else rbind(aod, aodf[-1, , drop=FALSE])
        }
        attr(aod, "heading") <- NULL
        if(is.null(aod) || ncol(aod) == 0) break
        ## need to remove any terms with zero df from consideration
        nzdf <- if(!is.null(aod$Df)) aod$Df != 0 | is.na(aod$Df)
        aod <- aod[nzdf, ]
        if(is.null(aod) || ncol(aod) == 0) break
        nc <- match(c("Cp", "AIC"), names(aod))
        nc <- nc[!is.na(nc)][1L]
        o <- order(aod[, nc])
        if(trace) {
          print(aod[o,  ])
          utils::flush.console()
        }
        if(o[1L] == 1) break
        change <- rownames(aod)[o[1L]]
      }
      usingCp <- match("Cp", names(aod), 0) > 0
      ## may need to look for a 'data' argument in parent
      fit <- update(fit, paste("~ .", change), evaluate = FALSE)
      fit <- eval.parent(fit)
      nnew <- nobs(fit, use.fallback = TRUE)
      if(all(is.finite(c(n, nnew))) && nnew != n)
        stop("number of rows in use has changed: remove missing values?")
      Terms <- terms(fit)
      bAIC <- extractAIC(fit, scale, k = k, ...)
      edf <- bAIC[1L]
      bAIC <- bAIC[2L]
      if(trace) {
        cat("\nStep:  AIC=", format(round(bAIC, 2)), "\n",
            cut.string(deparse(formula(fit))), "\n\n", sep='')
        utils::flush.console()
      }
      ## add a tolerance as dropping 0-df terms might increase AIC slightly
      if(bAIC >= AIC + 1e-7) break
      nm <- nm + 1
      models[[nm]] <-
        list(deviance = mydeviance(fit), df.resid = n - edf,
             change = change, AIC = bAIC)
      if(!is.null(keep)) keep.list[[nm]] <- keep(fit, bAIC)
    }
    if(!is.null(keep)) fit$keep <- re.arrange(keep.list[seq(nm)])
    step.results(models = models[seq(nm)], fit, object, usingCp)
  }
#逐步回归后退法
cal.stepwise<-function(fit){
  s2<-stepAICC(fit,direction='both',trace=0)
  coe<-s2$coefficients
  l<-attr(coe,"names")
  len<-length(l)
  l[2:len]
}

#组装模型
pack.model<-function(form,coee){
  chrlist<-division.form(form)
  n<-length(coee)
  coee.new<-coee[2:n]
  
  xmode<-paste(coee.new,chrlist,collapse="+",sep="")
  #  str(xmode)
  xmode<-paste(coee[1],'+',xmode,sep="")
  #  print(xmode)
  formnew<-as.formula(paste("Y~",xmode))
  formnew
}

#写日志文件函数
write.opt<-function(filename,name,opt,boot,test){
  s1<-paste("BOOT剩余变量：",paste(name,collapse=","))
  optname<-list('AUC:',"R2:","slope:")
  s2<-paste("boot校正量：",paste(optname,opt,collapse=" ",sep=" "))
  bootname<-list('AUC:',"R2:","slope:")
  s3<-paste("boot训练量：",paste(optname,boot,collapse=" ",sep=" "))
  testname<-list('AUC:',"R2:","slope:")
  s4<-paste("boot测试量：",paste(optname,test,collapse=" ",sep=" "))
  write.table(s1,file=filename,append=TRUE,col.names=FALSE,row.names=FALSE)
  write.table(s3,file=filename,append=TRUE,col.names=FALSE,row.names=FALSE)
  write.table(s4,file=filename,append=TRUE,col.names=FALSE,row.names=FALSE)
  write.table(s2,file=filename,append=TRUE,col.names=FALSE,row.names=FALSE)
}

#用来写boot和test的日志文件
write.boot<-function(filename,boot,opt){
  boot.name<-list('boot.auc:','boot.r2','boot.slope')
  test.name<-list('test.auc','test.r2','test,slope')
  s1<-paste("Bootstrap performance:",paste(boot.name,boot,collapse=" ",sep=" "))
  s2<-paste("Test performance:",paste(test.name,opt,collapse=" ",sep=" "))
  write.table(s1,file=filename,append=TRUE,col.names=FALSE,row.names=FALSE)
  write.table(s2,file=filename,append=TRUE,col.names=FALSE,row.names=FALSE)
}

#写Apparent日志文件
wite.apparent<-function(filename,appar){
  app.name<-list('apparent.auc:','apparent.r2','apparent.slope')
  s1<-paste("Apparent performance:",paste(app.name,appar,collapse=" ",sep=" "))
  write.table(s1,file=filename,append=TRUE,col.names=FALSE,row.names=FALSE)
}

#写expected日志文件
wite.expected<-function(filename,appar){
  app.name<-list('expected.auc:','expected.r2','expected.slope')
  s1<-paste("expected optimism:",paste(app.name,appar,collapse=" ",sep=" "))
  write.table(s1,file=filename,append=TRUE,col.names=FALSE,row.names=FALSE)
}

#写Optimism-corrected日志文件
write.correct<-function(filename,appar){
  app.name<-list('correct.auc:','correct.r2','correct.slope')
  s1<-paste("Optimism-corrected:",paste(app.name,appar,collapse=" ",sep=" "))
  write.table(s1,file=filename,append=TRUE,col.names=FALSE,row.names=FALSE)
}

#写出文件头
writer.bianliangtou<-function(appar){
  s1<-paste(appar,collapse=",")
  s1
}

#写单变量日志文件
write.dan<-function(filename,appar){
  s1<-paste("单变量分析后剩余变量：",paste(appar,collapse=" "))
  write.table(s1,file=filename,append=TRUE,col.names=FALSE,row.names=FALSE)
}

#写集中文件
write.jizhong1<-function(filename,appar){
  write.table(appar,file=filename,col.names = TRUE,row.names = FALSE,sep = ",",eol = "\n")
}

#写集中日志文件
write.jizhong2<-function(filename,appar){
  write.table(appar,file=filename,append=TRUE,col.names = FALSE,row.names = FALSE,sep = ",",eol = "\n")
}

#计算出相应的标准差
cal.sd<-function(fileread){
  mydata<-read.table(file=fileread,header = TRUE,sep=',')
  boot.auc<-mydata$Bootstrap.Auc
  boot.r2<-mydata$Bootstrap.R2
  boot.slope<-mydata$Bootstrap.slope
  test.auc<-mydata$Test.AUC
  test.r2<-mydata$Test.R2
  test.slope<-mydata$Test.slope
  opt.auc<-mydata$Expected.optimism.AUC
  opt.r2<-mydata$Expected.optimism.R2
  opt.slope<-mydata$Expected.optimism.slope
  true.auc<-mydata$Optimism.corrected.AUC
  true.r2<-mydata$Optimism.corrected.R2
  true.slope<-mydata$Optimism.corrected.Slope
  boot.aucsd<-sd(boot.auc)
  boot.r2sd<-sd(boot.r2)
  boot.slopesd<-sd(boot.slope)
  test.aucsd<-sd(test.auc)
  test.r2sd<-sd(test.r2)
  test.slopesd<-sd(test.slope)
  opt.aucsd<-sd(opt.auc)
  opt.r2sd<-sd(opt.r2)
  opt.slopesd<-sd(opt.slope)
  true.aucsd<-sd(true.auc)
  true.r2sd<-sd(true.r2)
  true.slopesd<-sd(true.slope)
  boot.sd<-c(boot.aucsd,boot.r2sd,boot.slopesd)
  names(boot.sd)<-c("AUC","R2","Slope")
  test.sd<-c(test.aucsd,test.r2sd,test.slopesd)
  names(test.sd)<-c("AUC","R2","Slope")
  opt.sd<-c(opt.aucsd,opt.r2sd,opt.slopesd)
  names(opt.sd)<c("AUC","R2","Slope")
  true.sd<-c(true.aucsd,true.r2sd,true.slopesd)
  names(true.sd)<c("AUC","R2","Slope")
  z<-c(boot.sd,test.sd,opt.sd,true.sd)
  z
}
#calcum apparent data
calcum.apparent<-function(form,mydata,file){
  fit<-glm(form,data=mydata,family = binomial())
  name.list<-cal.stepwise(fit)
  name.list<-division.form(form)#使用全模型
  form.v1<-pack.form(name.list)
  fitnew<-lrm(form.v1,x=TRUE, y=TRUE, data = mydata)
  print(fitnew)
  par<-apparent.factor(form.v1,mydata)
  write.csv(fitnew$coefficients,file)
  z<-c(par,name.list)#是否输出变量
  print(fitnew$coefficients)
  z
}

#bootstrap过程，其中有单因素分析
boot.optimate<-function(formal,mydata,var.univ,mult.in,mult.out,filezong,bootnum,app.z){
  form<-formal
  auc.boot<-0.0
  R2.boot<-0.0
  sloop.boot<-0.0
  auc.test<-0.0
  R2.test<-0.0
  sloop.test<-0.0
  z.boot<-c(auc.boot,R2.boot,sloop.boot)
  names(z.boot)<-c("AUC","R2","Slope")
  z.test<-c(auc.test,R2.test,sloop.test)
  names(z.test)<-c("AUC","R2","Slope")
  num<-bootnum
  #boot过程
  for(i in 1:bootnum){
    a<-runif(10000,min=1,max=1000)
    set.seed(13*i*a[i])
    bootdata<-boot.sample(mydata)
    val.list<-univar.analy(formal,var.univ,bootdata)#单变量分析过程
    val.retu<-val.list
    #   write.dan(filedan,val.list)#写单变量日志文件
    vlen<-length(val.list)
    if(vlen==0){
      boot.par<-c(0,0,0)
      test.par<-c(0,0,0)
      names(boot.par)<-c("AUC","R2","Slope")
      names(test.par)<-c("AUC","R2","Slope")
      opti.par<-boot.par-test.par
      #     write.opt(filewriter,name.list,opti.par,boot.par,test.par)#写文档
      z.boot<-z.boot+boot.par
      z.test<-z.test+test.par
      num<-num-1
      
    }
    else{
      formt<-pack.form(val.list)#组装模型
      fit<-lrm(formt,x=TRUE, y=TRUE, data = bootdata)
      name.list<-cal.stepwise(fit,mult.in,mult.out)
      name.retu<-name.list
      nlen<-length(name.list)
      if(nlen==0){
        boot.par<-c(0,0,0)
        test.par<-c(0,0,0)
        names(boot.par)<-c("AUC","R2","Slope")
        names(test.par)<-c("AUC","R2","Slope")
        opti.par<-boot.par-test.par
        #       write.opt(filewriter,name.list,opti.par,boot.par,test.par)#写文档
        z.boot<-z.boot+boot.par
        z.test<-z.test+test.par
        num<-num-1
      }
      else{
        form.v1<-pack.form(name.list)
        fitnew<-lrm(form.v1,x=TRUE, y=TRUE, data = bootdata)
        boot.par<-calcum.boot(fitnew)#boot的三个指标
        test.par<-calcum.test(fitnew,mydata,form.v1)
        opti.par<-boot.par-test.par
        #     write.table(name.list,file=filewriter,append=TRUE,col.names=FALSE,row.names=FALSE)
        #     write.table(opti.par,file=filewriter,append=TRUE,col.names=FALSE,row.names=FALSE)
        #       write.opt(filewriter,name.list,opti.par,boot.par,test.par)#写文档
        z.boot<-z.boot+boot.par
        z.test<-z.test+test.par
      }
    }
    s1<-as.character(writer.bianliangtou(val.retu))
    s2<-as.character(writer.bianliangtou(name.retu))
    if(i==1){ 
      true.z<-app.z-opti.par
      stemp<-list(i,s2,boot.par[1],test.par[1]
                  ,opti.par[1],true.z[1],boot.par[2],test.par[2],opti.par[2],true.z[2]
                  ,boot.par[3],test.par[3],opti.par[3],true.z[3])
      #print(stemp)
      sframe<-as.data.frame(stemp)
      colnames(sframe)<-c('bootstrap 次数','多因素逐步回归选择剩余变量','Bootstrap Auc','Test AUC'
                          ,'Expected optimism AUC','Optimism corrected AUC'
                          ,'Bootstrap R2','Test R2','Expected optimism R2','Optimism corrected R2'
                          ,'Bootstrap slope','Test slope','Expected optimism slope','Optimism corrected Slope')
      write.jizhong1(filezong,sframe)
    }
    else{
      true.z<-app.z-opti.par
      stemp<-list(i,s2,boot.par[1],test.par[1]
                  ,opti.par[1],true.z[1],boot.par[2],test.par[2],opti.par[2],true.z[2]
                  ,boot.par[3],test.par[3],opti.par[3],true.z[3])
      sframe<-as.data.frame(stemp)
      write.jizhong2(filezong,sframe) 
    }
  }
  
  new.boot<-z.boot/num
  new.test<-z.test/num
  new.opt<-new.boot-new.test
  #  new.opt #
  z<-c(new.boot,new.test,new.opt)
  z
}

#有单因素分析
boot.LR<-function(form,fileread,fileframe,filemodel,filezong,var.univ,mult.in,mult.out,bnum){
  mydata<-read.table(file=fileread,header=TRUE,sep=",")
  # form=Y~X1+X2+X3+X4+X7+X10+X11+X12+X13+X14+X15
  app.z<-calcum.apparent(form,mydata,var.univ,mult.in,mult.out,filemodel)
  nlen<-length(app.z)
  form.v1<-pack.form(app.z[4:nlen])
  len<-length(app.z)
  auc<-as.numeric(app.z[1])
  r2<-as.numeric(app.z[2])
  slope<-as.numeric(app.z[3])#需要转化为数值型
  app.z<-c(auc,r2,slope)#原始数据
  z<-boot.optimate(form,mydata,var.univ,mult.in,mult.out,filezong,bnum,app.z)
  opt.z<-z[7:9]
  boot.z<-z[1:3]
  test.z<-z[4:6]
  true.z<-app.z-opt.z
  name<-c('predict factor','Apparent performance','Bootstrap performance'," Test performance","Excepted optimism","Optimism-corrected")
  cs<-cal.sd(filezong)
  boot.sd<-cs[1:3]
  test.sd<-cs[4:6]
  opt.sd<-cs[7:9]
  true.sd<-cs[10:12]
  frame<-data.frame(app.z,boot.z,boot.sd,test.z,test.sd,opt.z,opt.sd,true.z,true.sd)
  write.table(frame,fileframe,sep=",",col.names = NA,qmethod = "double")
}

#bootstrap过程，无单因素分析
boot.optimateUnvali<-function(form,mydata,file,filezong,bootnum,app.z){
  auc.boot<-0.0
  R2.boot<-0.0
  sloop.boot<-0.0
  auc.test<-0.0
  R2.test<-0.0
  sloop.test<-0.0
  z.boot<-c(auc.boot,R2.boot,sloop.boot)
  names(z.boot)<-c("AUC","R2","Slope")
  z.test<-c(auc.test,R2.test,sloop.test)
  names(z.test)<-c("AUC","R2","Slope")
  #boot过程
  num<-bootnum
  for(i in 1:bootnum){
    a<-runif(10000,min=1,max=10000)
    set.seed(13*a[i])
    bootdata<-boot.sample(mydata)
    fit<-glm(form,data=bootdata,family = binomial())
    name.list<-cal.stepwise(fit)
    name.list<-division.form(form)
    name.retu<-name.list
    nlen<-length(name.list)
    if(nlen==0){
      boot.par<-c(0,0,0)
      test.par<-c(0,0,0)
      names(boot.par)<-c("AUC","R2","Slope")
      names(test.par)<-c("AUC","R2","Slope")
      opti.par<-boot.par-test.par
      z.boot<-z.boot+boot.par
      z.test<-z.test+test.par
      num<-num-1
    }
    else{
      form.v1<-pack.form(name.list)
      fitnew<-lrm(form.v1,x=TRUE, y=TRUE, data = bootdata)
      boot.par<-calcum.boot(fitnew)#boot的三个指标
      test.par<-calcum.test(fitnew,mydata,form.v1)
      opti.par<-boot.par-test.par
      z.boot<-z.boot+boot.par
      z.test<-z.test+test.par
    }
    s2<-as.character(writer.bianliangtou(name.retu))
    if(i==1){ 
      true.z<-app.z-opti.par
      stemp<-list(i,s2,boot.par[1],test.par[1]
                  ,opti.par[1],true.z[1],boot.par[2],test.par[2],opti.par[2],true.z[2]
                  ,boot.par[3],test.par[3],opti.par[3],true.z[3])
      #print(stemp)
      sframe<-as.data.frame(stemp)
      colnames(sframe)<-c('bootstrap 次数','多因素逐步回归选择剩余变量','Bootstrap Auc','Test AUC'
                          ,'Expected optimism AUC','Optimism corrected AUC'
                          ,'Bootstrap R2','Test R2','Expected optimism R2','Optimism corrected R2'
                          ,'Bootstrap slope','Test slope','Expected optimism slope','Optimism corrected Slope')
      write.jizhong1(filezong,sframe)
    }
    else{
      true.z<-app.z-opti.par
      stemp<-list(i,s2,boot.par[1],test.par[1]
                  ,opti.par[1],true.z[1],boot.par[2],test.par[2],opti.par[2],true.z[2]
                  ,boot.par[3],test.par[3],opti.par[3],true.z[3])
      sframe<-as.data.frame(stemp)
      write.jizhong2(filezong,sframe) 
    }
  }
  new.boot<-z.boot/num
  new.test<-z.test/num
  new.opt<-new.boot-new.test
  #  write.boot(file,new.boot,new.test)
  #  wite.expected(file,new.opt)
  #  new.opt #
  z<-c(new.boot,new.test,new.opt)
  z
}
############################
#无单因素分析
boot.LRUnvali<-function(form,fileread,fileframe,filemodel,filezong,bnum){
  mydata<-read.table(file=fileread,header=TRUE,sep=",")
  app.z<-calcum.apparent(form,mydata,filemodel)#计算出原始值  使用全模型
  nlen<-length(app.z)
  form.v1<-pack.form(app.z[4:nlen])
  len<-length(app.z)
  auc<-as.numeric(app.z[1])
  r2<-as.numeric(app.z[2])
  slope<-as.numeric(app.z[3])#需要转化为数值型
  app.z<-c(auc,r2,slope)#原始数据
  z<-boot.optimateUnvali(form,mydata,filemodel,filezong,bnum,app.z)
  opt.z<-z[7:9]
  boot.z<-z[1:3]
  test.z<-z[4:6]
  true.z<-app.z-opt.z
  name<-c('predict factor','Apparent performance','Bootstrap performance'," Test performance","Excepted optimism","Optimism-corrected")
  cs<-cal.sd(filezong)
  boot.sd<-cs[1:3]
  test.sd<-cs[4:6]
  opt.sd<-cs[7:9]
  true.sd<-cs[10:12]
  frame<-data.frame(app.z,boot.z,boot.sd,test.z,test.sd,opt.z,opt.sd,true.z,true.sd)
  write.table(frame,fileframe,sep=",",col.names = NA,qmethod = "double")
}



#哈医大bootstrap整个过程的函数，其中form是调用的公式，fileread文件是存放数据的文件，如本模型中数据.csv，其中要注意格式
#fileFrame是输出的frame文件的名称，记录实验的结果
#fileopt是输出opti文件，记录的每次bootstarp剩余变量的个数，以及相应的AUC ,R2 Slope的值
#filerecord是输出record文件，记录一些记录文件
#filedan是输出danbianl文件，记录每次bootstarap中单因素分析后剩余变量
#filezong是输出zonghe 文件，是将整个综合起来的文件
#var.univ是单因素分析的参数
#mult.in mult.out是逐步回归后退法的参数
#bnum是bootstrap抽样次数
#是否需要单因素分析，若TRUE则需要单因素分析，若FALSE则不需要单因素分析
hayida<-function(form,fileread,fileframe,filemodel,filezong,bnum){
  boot.LRUnvali(form,fileread,fileframe,filemodel,filezong,bnum)
}

#画出每组实验的slope图形
draw.slope<-function(form,filepng,fileread,var.univ=0.2,mult.in=0.05,mult.out=0.05,num.ovfit=1){
  mydata<-read.table(file=fileread,header=TRUE,sep=",")
  #  form<-Y~X1+X2+X3+X4+X7+X10+X11+X12+X13+X14+X15
  val.list<-univar.analy(form,var.univ,mydata)
  form<-pack.form(val.list)
  # fit<-lrm(form,x=TRUE, y=TRUE, data = mydata)
  fit<-glm(form,data=bootdata,family = binomial())
  name.list<-cal.stepwise(fit,mult.in,mult.out)
  form.v1<-pack.form(name.list)
  fitnew<-lrm(form.v1,x=TRUE, y=TRUE, data = mydata)
  p<-predict(fitnew,mydata)
  p<-p/num.ovfit
  pred=exp(p)/(1+exp(p))
  jpeg(file=filepng,bg='white',width=600,height = 600,quality = 600)
  val.prob(p=pred,mydata$Y,logistic.cal = FALSE,statloc = c(1,3))
  dev.off()
}

#计算出对比模型model2的slope值
cal.model2slope<-function(fileread){
  mydata<-read.table(file='预测数据.csv',header=TRUE,sep=",")
  pred<-mydata$预测
  v.true<-mydata$真实
  y<-exp(v.true/(1-v.true))
  x<-exp(pred/(1-pred))
  fit<-lm(y~x)
  slope<-fit$coefficients
  slope
}

#画出对比模型model2的slope图形
draw.model2slope<-function(fileread,filepng){
  mydata<-read.table(file=fileread,header=TRUE,sep=",")
  pred<-mydata$预测
  jpeg(file=filepng,bg='white',width=600,height = 600,quality = 600)
  val.prob(p=pred,mydata$Y,logistic.cal = FALSE,statloc = c(1,3))
  dev.off()
}

#定义函数接口
hayida.bootstrap<-function(form,fileread,fileout='data',bnum){
  form=as.formula(form)
  if(!file.exists(fileread)){
    print("请输入数据文件路径")
  }
  else{
    if(length(list.dirs(fileout))==0){
      dir.create(fileout)#如果该文件夹不存在，则创建相应的文件夹
    } 
    fileframe=paste0(fileout,'/final result.csv')
    filemodel=paste0(fileout,'/model.csv')
    filezong=paste0(fileout,'/bootstrap.csv')
    hayida(form,fileread,fileframe,filemodel,filezong,bnum)
    myres<-read.table(file=fileframe,header=TRUE,sep=",")
    true.z<-myres$true.z
    opt.slope<-true.z[3]
    filepng=paste0(fileout,'/slope.png')
    #   draw.slope(form,filepng,fileread,var.univ,mult.in,mult.out,opt.slope)#画出相应的slope图形
  }
}

#提取对应的阀值点
cal.auc<-function(mydata){
  data<-mydata[order(mydata$P,decreasing=F),]#以概率P从小到大进行排序
  len.1<-length(which(data$Y==1))
  
  len.0<-length(data$Y)-len.1
  
  xstep<-1.0/len.0
  ystep<-1.0/len.1
  #  print(data)
  max<--10.0
  max.index<-0
  x<-1.0
  y<-1.0
  label<-data$Y#取出标签值
  xvalue<-data$X
  temp<-0
  #这里有问题存在
  for(i in 1:length(data$Y)){
    if(label[i]!=1){
      x<-x-xstep
      y<-y
      if((y-x)>0)
        temp<-y-x
      else
        temp<-x-y
      len<-temp
      if(len>max){
        max.index<-i
        max<-len
        
      }
    }
    else{
      x<-x
      y<-y-ystep
      if((y-x)>0)
        temp<-y-x
      else
        temp<-x-y
      len<-temp
      if(len>max){
        max.index=i
        max=len
      }
    }
  }
  return (xvalue[max.index])
}
#对单列数据进行二分类
cal.binary<-function(mydata,thre,xname){
  len<-nrow(mydata)#求出数据框的行数
  for(i in 1:len){
    if(mydata[i,xname]>thre){
      if(xname=='X10'||xname=='X13'){
        mydata[i,xname]<-0
      }
      else{
        mydata[i,xname]<-1
      }
      
    }
    
    else{
      if(xname=='X10'||xname=='X13'){
        mydata[i,xname]<-1
      }
      else{
        mydata[i,xname]<-0
      }
    }
  }
  mydata
}
#将数据进行二分类
cal.binaryAll<-function(mydata,form){
  #只对这四个变量进行二分类
  list.x<-division.form(form)#相应元素
  len<-length(list.x)
  for(i in 1:len){
    x<-list.x[i]
    #判断是否需要二分类
    if(length(unique(mydata[,x]))!=2){
      data1<-mydata[,c(x,'Y')]#提取对应的两列
      names(data1)[1]<-"X"#将名称改掉
      fit<-lrm(Y~X,x=TRUE, y=TRUE, data = data1)
      pred<-predict(fit,data1)#预测概率
      data1$P=pred#添加新列
      thre<-cal.auc(data1)
      mydata<-cal.binary(mydata,thre,x)
    }
  }
  mydata
}
#统计有多少不一样的数据
cal.statDiff<-function(mydata,data){
  m<-264
  n<-13
  for(j in 1:n){
    count<-0
    for(i in 1:m){
      if(mydata[i,j]!=data[i,j]){
        count<-count+1
      }
    }
    print(count)
  }
}
#如何对数据进行二分类

form=Y~X1+X2+X3+X4+X7+X10+X11+X12+X13+X14+X15
mydata<-read.table(file='混合数据.csv',header=TRUE,sep=',')
bootdata<-boot.sample(mydata)
hayida.bootstrap(form,'混合数据.csv','data/v35',bnum=500)
