
## 1. CDFs -- LDA

require(aucdiff)
require(mvtnorm)
require(parallel)
start <- Sys.time()
set.seed(0)
alpha <- .05
mc.reps <- 1e2
p.f <- 6
p.r <- round(p.f*2/3)
Sigma.f <- diag(p.f)
Sigma.r <- Sigma.f[1:p.r,1:p.r]
beta.star.r <- runif(p.r)
auc.true <- .75
beta.star.r <- beta.star.r / sqrt(t(beta.star.r)%*%Sigma.r%*%beta.star.r) * qnorm(auc.true)
beta.star.f <- c(beta.star.r,rep(0,p.f-p.r))
mu.0.f <- rep(0,p.f)
mu.1.f <- c(Sigma.f%*%beta.star.f)
start <- Sys.time()
Ns <- c(50,150)
settings <- cbind(N=Ns)
by.setting <- mclapply(1:nrow(settings), mc.cores=detectCores()-2, FUN=function(i) {
    N <- settings[i,'N']
    m <- round(2/3*N)
    n <- N-m
    pi.0 <- m/(m+n); pi.1 <- 1-pi.0
    ## empirical
    obs <- replicate(mc.reps, {
        x.0.f <- rmvnorm(m,mu.0.f,Sigma.f)
        x.1.f <- rmvnorm(n,mu.1.f,Sigma.f)
        gamma.hat <- c(coefs.lda(x.0.f,x.1.f,params=list()))
        beta.hat <- coefs.lda(x.0.f[,1:p.r],x.1.f[,1:p.r],params=list())
        beta.hat <- c(beta.hat,rep(0,p.f-p.r))
        auc.hat(x.0.f%*%beta.hat,x.1.f%*%beta.hat) - auc.hat(x.0.f%*%gamma.hat,x.1.f%*%gamma.hat)  
    })
    ## oracle
    knots <- sort(obs)
    params.f   <- params <-  list(mu.0=mu.0.f,mu.1=mu.1.f,Sigma.0=Sigma.f,Sigma.1=Sigma.f,pi.0=pi.0,p.f=p.f,p.r=p.r)
    fit.oracle <- aucdiff(w=NULL,d=NULL,N=m+n,beta=beta.star.f,p.reduced=p.r,infl.fn=infl.lda,vcov.fn=var.normal.lda,hessian.fn=auc.hessian.normal,eps.finite.diff=NULL,dgm.params=params.f)
    ecdf.oracle <- paucdiff(knots,fit.oracle)
    ## parametric and nonparametric
    fits.estimated <- replicate(mc.reps, simplify=FALSE, expr={
        x.0.f <- rmvnorm(m,mu.0.f,Sigma.f)
        x.1.f <- rmvnorm(n,mu.1.f,Sigma.f)
        beta.hat <- coefs.lda(x.0.f[,1:p.r],x.1.f[,1:p.r],params=list(Sigma=Sigma.f[1:p.r,1:p.r]))
        beta.hat <- c(beta.hat,rep(0,p.f-p.r))
        gamma.hat <- c(coefs.lda(x.0.f,x.1.f,params=list(Sigma=Sigma.f)))
        w <- rbind(x.0.f,x.1.f)
        d <- c(rep(0,nrow(x.0.f)),rep(1,nrow(x.1.f)))
        Sigma.hat.f <- pi.0*var(x.0.f)+pi.1*var(x.1.f)
        params.f.hat <- list(mu.0=colMeans(x.0.f),mu.1=colMeans(x.1.f),Sigma.0=Sigma.hat.f,Sigma.1=Sigma.hat.f,pi.0=pi.0,pi.1=pi.1,p.r=p.r,p.f=p.f)
        fit.para <- aucdiff(w=NULL,d=NULL,N=m+n,beta=beta.hat,p.reduced=p.r,infl.fn=infl.lda,vcov.fn=var.normal.lda,hessian.fn=auc.hessian.normal,eps.finite.diff=NULL,dgm.params=params.f.hat)
        fit.nonpara <- aucdiff(w=w,d=d,beta=beta.hat,p.reduced=p.r,infl.fn=infl.lda,vcov.fn=NULL,hessian.fn=NULL,eps.finite.diff=.7,dgm.params=params.f.hat)
        list(para=fit.para,nonpara=fit.nonpara)
    })
    fits.para <- lapply(fits.estimated,function(lst)lst$para)
    ecdf.para  <- sapply(fits.para, function(fit.para) paucdiff(knots,fit.para))
    fits.nonpara <- lapply(fits.estimated,function(lst)lst$nonpara)
    ecdf.nonpara  <- sapply(fits.nonpara, function(fit.nonpara) paucdiff(knots,fit.nonpara))
    list(obs=obs,oracle=ecdf.oracle,para=rowMeans(ecdf.para),nonpara=rowMeans(ecdf.nonpara))
})
op <- par(mfrow=c(1,2))
for(i in 1:length(by.setting)) {
    with(by.setting[[i]], {
        knots <- sort(obs)
        plot(ecdf(obs),main=paste0('M+N=',settings[i,'N']),xlab='',ylab='')
        lines(knots,oracle,lty=1)
        lines(knots,para,lty=2)
        lines(knots,nonpara,lty=3)
        legend('bottomright',lty=1:4,legend=c('oracle','para','nonpara'))
    })
}
par(op)


extrafont::loadfonts()
png('../figs/cdf_logit.png', width = 1024, height = 768, pointsize=15, family='CM Roman')
op <- par(mfrow=c(1,2))
for(i in 1:length(by.setting)) {
    with(by.setting[[i]], {
        knots <- sort(obs)
        plot(ecdf(obs),main=paste0('M+N=',settings[i,'N']),xlab='',ylab='')
        lines(knots,oracle,lty=1)
        lines(knots,para,lty=2)
        lines(knots,nonpara,lty=3)
    })
}
par(op)
dev.off()





## 2. CDFs -- logit

require(aucdiff)
require(mvtnorm)
require(parallel)
start <- Sys.time()
set.seed(1)
infl.logit <- with(list(h=plogis), {
    h.1 <- function(x)h(x)*(1-h(x))
    function(w,d)infl.glm(w,d,params=list(link=plogis,link.deriv=function(x)h(x)*(1-h(x)),link.deriv2=function(x)h.1(x)*(1-2*h(x)),link.name='logit'))
    })
mc.reps <- 5e2
alpha <- .05
p.f <- 6
p.r <- round(p.f*2/3)
beta.star.r <- runif(p.r)
beta.star.f <- c(beta.star.r,rep(0,p.f-p.r))
start <- Sys.time()
Ns <- c(50,150)
settings <- cbind(N=Ns)
by.setting <- sapply(1:nrow(settings), simplify=FALSE, FUN=function(i) {
    N <- settings[i,'N']
    obs <- replicate(mc.reps, {
        w <- x.f <- matrix(rnorm(p.f*N),ncol=p.f)
        risk <- plogis(x.f%*%beta.star.f)
        d <- rbinom(N,1,risk)
        x.r <- x.f[,1:p.r]
        x.0.f <- x.f[d==0,]; x.1.f<- x.f[d==1,]    
        x.0.r <- x.r[d==0,]; x.1.r <- x.r[d==1,]    
        gamma.hat <- coef(glm(d~x.f-1,family=binomial(link='logit')))
        beta.hat <- coef(glm(d~x.r-1,family=binomial(link='logit')))
        beta.hat <- c(beta.hat,rep(0,p.f-p.r))
        delta.hat <- auc.hat(x.0.f%*%beta.hat,x.1.f%*%beta.hat) - auc.hat(x.0.f%*%gamma.hat,x.1.f%*%gamma.hat)  
    })
    fits.nonpara <- replicate(mc.reps, simplify=FALSE, expr={
        w <- x.f <- matrix(rnorm(p.f*N),ncol=p.f)
        risk <- plogis(x.f%*%beta.star.f)
        d <- rbinom(N,1,risk)
        x.r <- x.f[,1:p.r]
        x.0.f <- x.f[d==0,]; x.1.f<- x.f[d==1,]    
        x.0.r <- x.r[d==0,]; x.1.r <- x.r[d==1,]    
        gamma.hat <- coef(glm(d~x.f-1,family=binomial(link='logit')))
        beta.hat <- coef(glm(d~x.r-1,family=binomial(link='logit')))
        beta.hat <- c(beta.hat,rep(0,p.f-p.r))
        fit.nonpara <-  aucdiff(w=w,d=d,beta=beta.hat,p.reduced=p.r,infl.fn=infl.logit,vcov.fn=NULL,hessian.fn=NULL,eps.finite.diff=exp(2)/sqrt(N),dgm.params=NULL)
    })
    knots <- sort(obs)
    cdf.hats  <- sapply(fits.nonpara, function(fit.nonpara) paucdiff(knots,fit.nonpara))
    list(obs=obs,cdf.hats=cdf.hats)
})
print(Sys.time() - start)
op <- par(mfrow=c(1,2))
for(i in 1:length(by.setting)) {
    with(by.setting[[i]], {
        plot(ecdf(obs),main=paste0('N=',settings[i,'N']),xlab='',ylab='')
        nonpara.mean <- rowMeans(cdf.hats)
        knots <- sort(obs)
        nonpara.CI <- apply(cdf.hats,1,quantile,probs=c(alpha/2,1-alpha/2))
        lines(knots,nonpara.mean,lty=3)        
        lines(knots,nonpara.CI[1,],lty=3)        
        lines(knots,nonpara.CI[2,],lty=3)        
    })
}
par(op)

extrafont::loadfonts()
png('../figs/cdf_logit.png', width = 1024, height = 768, pointsize=15, family='CM Roman')
op <- par(mfrow=c(1,2))
for(i in 1:length(by.setting)) {
    with(by.setting[[i]], {
        plot(ecdf(obs),main=paste0('M+N=',settings[i,'N']),xlab='',ylab='')
        nonpara.mean <- rowMeans(cdf.hats)
        knots <- sort(obs)
        nonpara.CI <- apply(cdf.hats,1,quantile,probs=c(alpha/2,1-alpha/2))
        lines(knots,nonpara.mean,lty=3)        
        lines(knots,nonpara.CI[1,],lty=3)        
        lines(knots,nonpara.CI[2,],lty=3)        
    })
}
par(op)
dev.off()








## 3. power -- logit

require(aucdiff)
require(mvtnorm)
set.seed(0)
start <- Sys.time()
p.f <- 6
p.r <- round(p.f*2/3)
h <- plogis
h.1 <- function(x)h(x)*(1-h(x))
infl.logit <- function(w,d)infl.glm(w,d,params=list(link=plogis,link.deriv=function(x)h(x)*(1-h(x)),link.deriv2=function(x)h.1(x)*(1-2*h(x)),link.name='logit'))
beta.star.r <- rnorm(p.r)
alpha <- .05
start <- Sys.time()
Ns <- c(50,150)
offsets <- rev(seq(-1,0,len=20))
offsets <- seq(-1,1,len=21)
settings <- lapply(Ns, function(N)cbind(N=N,offset=offsets*10/sqrt(N)))
settings <- do.call(rbind,settings)
by.setting <- parallel::mclapply(1:nrow(settings), mc.cores=5, FUN=function(i) {
    N <- settings[i,'N']
    offset <- settings[i,'offset']
    print(N)
    print(offset)
    beta.star.f <- c(beta.star.r,rep(offset,p.f-p.r))
    obs <- replicate(1e3, {
        w <- x.f <- matrix(rnorm(p.f*N),ncol=p.f)
        risk <- h(x.f%*%beta.star.f)
        d <- rbinom(N,1,risk)
        x.r <- x.f[,1:p.r]
        x.0.f <- x.f[d==0,]; x.1.f<- x.f[d==1,]    
        x.0.r <- x.r[d==0,]; x.1.r <- x.r[d==1,]    
        gamma.hat <- coef(glm(d~x.f-1,family=binomial(link='logit')))
        beta.hat <- coef(glm(d~x.r-1,family=binomial(link='logit')))
        beta.hat <- c(beta.hat,rep(0,p.f-p.r))
        delta.hat <- auc.hat(x.0.f%*%beta.hat,x.1.f%*%beta.hat) - auc.hat(x.0.f%*%gamma.hat,x.1.f%*%gamma.hat)  
        reject.nonpara <- tryCatch(
        {
            fit.nonpara <-  aucdiff(w=w,d=d,beta=beta.hat,p.reduced=p.r,infl.fn=infl.logit,vcov.fn=NULL,hessian.fn=NULL,eps.finite.diff=exp(2)/sqrt(N),dgm.params=NULL)
            CI <- confint(fit.nonpara,1-alpha)
            prod(CI - delta.hat) > 0
        }, error=function(e)NA,warning=function(w)NA)
        reject.delong <- suppressMessages(pROC::roc.test(d,c(w%*%beta.hat),c(w%*%gamma.hat),method='delong',smooth=FALSE)$p.value < alpha)
        c(non.para=reject.nonpara,delong=reject.delong)
    })
    rowMeans(obs,na.rm=TRUE)
})
print(Sys.time() - start)
beepr::beep()
by.setting <- simplify2array(by.setting)
cbind(settings,t(by.setting))


extrafont::loadfonts()
png(paste0('figs/power_logit.png'), width = 1024, height = 768, pointsize=15, family='CM Roman')
op <- par(mfrow=c(1,2))
for(i in 1:2) {
    df <- split(as.data.frame(t(by.setting)),settings[,'N'])[[i]]
    plot(0,type='n',xlim=range(offsets),ylim=range(df),ylab='rejection rate',xlab='offset from null',lty=4,main=paste0('N=',Ns[i]))
    for(j in 1:ncol(df)) lines(offsets,df[,j],lty=j)
    ## legend('bottomright',lty=1:ncol(df),legend=colnames(df))
    abline(h=alpha,lty=3)
}
dev.off()





## 4. power -- LDA
require(aucdiff)
require(mvtnorm)
require(parallel)
start <- Sys.time()
set.seed(4)
alpha <- .05
p.f <- 6
p.r <- round(p.f*2/3)
Sigma.f <- diag(p.f)
Sigma.r <- Sigma.f[1:p.r,1:p.r]
beta.star.r <- rnorm(p.r)
auc.true <- .75
beta.star.r <- beta.star.r / sqrt(t(beta.star.r)%*%Sigma.r%*%beta.star.r) * qnorm(auc.true)
pnorm(sqrt(t(beta.star.r)%*%Sigma.r%*%beta.star.r))
alpha <- .05
start <- Sys.time()
Ns <- c(50,150)
offsets <- seq(-1,1,len=21)
auc.by.offset <- sapply(offsets,function(offset) {
    beta.star.f <- c(beta.star.r,rep(offset,p.f-p.r))
    mu.0.f <- rep(0,p.f)
    mu.1.f <- c(Sigma.f%*%beta.star.f)
    auc.normal(beta.star.f,list(mu.0=mu.0.f,mu.1=mu.1.f,Sigma.0=Sigma.f,Sigma.1=Sigma.f))-auc.normal(beta.star.r,list(mu.0=mu.0.f[1:p.r],mu.1=mu.1.f[1:p.r],Sigma.0=Sigma.r,Sigma.1=Sigma.r))
    pnorm(sqrt(1/2*t(beta.star.f)%*%Sigma.f%*%beta.star.f)) -     pnorm(sqrt(1/2*t(beta.star.r)%*%Sigma.r%*%beta.star.r))
})
plot(offsets,auc.by.offset)
settings <- lapply(Ns, function(N)cbind(N=N,offset=offsets*10/sqrt(N)))
settings <- do.call(rbind,settings)
by.setting <- parallel::mclapply(1:nrow(settings), mc.cores=6, FUN=function(i) {
    N <- settings[i,'N']
    offset <- settings[i,'offset']
    print(N)
    print(offset)
    m <- round(2/3*N)
    n <- N-m
    pi.0 <- m/(m+n); pi.1 <- 1-pi.0
    beta.star.f <- c(beta.star.r,rep(offset,p.f-p.r))
    mu.0.f <- rep(0,p.f)
    mu.1.f <- c(Sigma.f%*%beta.star.f)
    obs <- replicate(2e3, {
        x.0.f <- rmvnorm(m,mu.0.f,Sigma.f)
        x.1.f <- rmvnorm(n,mu.1.f,Sigma.f)
        w <- rbind(x.0.f,x.1.f)
        d <- c(rep(0,nrow(x.0.f)),rep(1,nrow(x.1.f)))
        beta.hat <- coefs.lda(x.0.f[,1:p.r],x.1.f[,1:p.r],params=list())
        beta.hat <- as.numeric(c(beta.hat,rep(0,p.f-p.r))) 
        gamma.hat <- as.numeric(coefs.lda(x.0.f,x.1.f,params=list()))
        offset.hat <- auc.hat(x.0.f%*%beta.hat,x.1.f%*%beta.hat) - auc.hat(x.0.f%*%gamma.hat,x.1.f%*%gamma.hat)
        params.f <- list(mu.0=mu.0.f,mu.1=mu.1.f,Sigma.0=Sigma.f,Sigma.1=Sigma.f,pi.0=pi.0,p.r=p.r,p.f=p.f)
        Sigma.hat <- pi.0*cov(x.0.f)+pi.1*cov(x.1.f)
        params.hat.f <- list(mu.0=colMeans(x.0.f),mu.1=colMeans(x.1.f),Sigma.0=Sigma.hat,Sigma.1=Sigma.hat,m=m,n=n,pi.0=pi.0,p.r=p.r,p.f=p.f)
        fit.oracle <- aucdiff(w=NULL,d=NULL,N=m+n,beta=beta.star.f,p.reduced=p.r,infl.fn=infl.lda,vcov.fn=var.normal.lda,hessian.fn=auc.hessian.normal,eps.finite.diff=NULL,dgm.params=params.f)
        fit.para <- aucdiff(w=NULL,d=NULL,N=m+n,beta=gamma.hat,p.reduced=p.r,infl.fn=infl.lda,vcov.fn=var.normal.lda,hessian.fn=auc.hessian.normal,eps.finite.diff=NULL,dgm.params=params.hat.f)
        fit.nonpara <- aucdiff(w=w,d=d,beta=beta.hat,p.reduced=p.r,infl.fn=infl.lda,vcov.fn=NULL,hessian.fn=NULL,eps.finite.diff=exp(2)/sqrt(N),dgm.params=params.hat.f)
        reject.proposed <- sapply(list(oracle=fit.oracle,para=fit.para,nonpara=fit.nonpara), function(fit) {
            CI <- confint(fit,1-alpha)
            prod(CI - offset.hat) > 0
        })
        S <- m*cov(x.0.f)+n*cov(x.1.f)
        Sigma.hat <- (1/m+1/n)*S/(m+n)
        mu.1.bar.f <- colMeans(x.1.f); mu.0.bar.f <- colMeans(x.0.f)
        mu.1.bar.r <- mu.1.bar.f[1:p.r];     mu.0.bar.r <- mu.0.bar.f[1:p.r]
        maha.f <- t(mu.1.bar.f-mu.0.bar.f)%*%solve(Sigma.hat)%*%(mu.1.bar.f-mu.0.bar.f)
        maha.r <- t(mu.1.bar.r-mu.0.bar.r)%*%solve(Sigma.hat[1:p.r,1:p.r])%*%(mu.1.bar.r-mu.0.bar.r)
        f.stat <- (maha.f-maha.r) / (m+n+maha.r) * (m+n-p.f+1)/(p.f-p.r) 
        CI <- qf(c(alpha/2,1-alpha/2),p.f-p.r,n-p.f+1)
        reject.ftest <- prod(CI - f.stat) > 0
        reject.delong <- suppressMessages(pROC::roc.test(d,c(w%*%beta.hat),c(w%*%gamma.hat),method='delong',smooth=FALSE)$p.value < alpha)
        c(reject.proposed,ftest=reject.ftest,delong=reject.delong)
    })
    rowMeans(obs,na.rm=TRUE)## !!
})
by.setting <- simplify2array(by.setting)
print(Sys.time() - start)
beepr::beep()
cbind(settings,t(by.setting))
op <- par(mfrow=c(1,2))
for(i in 1:2) {
    df <- split(as.data.frame(t(by.setting)),settings[,'N'])[[i]]
    plot(0,type='n',xlim=range(offsets),ylim=range(df),ylab='rejection rate',xlab='distance from null',lty=4,main=paste0('N=',Ns[i]))
    for(j in 1:ncol(df)) lines(offsets,df[,j],lty=j)
    legend('bottomright',lty=1:ncol(df),legend=colnames(df))
    abline(h=alpha,lty=3)
}
par(op)


extrafont::loadfonts()
png(paste0('../figs/power_lda.png'), width = 1024, height = 768, pointsize=15, family='CM Roman')
op <- par(mfrow=c(1,2))
for(i in 1:2) {
    df <- split(as.data.frame(t(by.setting)),settings[,'N'])[[i]]
    plot(0,type='n',xlim=range(offsets),ylim=range(df),ylab='rejection rate',xlab='distance from null',lty=4,main=paste0('M+N=',Ns[i]))
    for(j in 1:ncol(df)) lines(offsets,df[,j],lty=j)
    abline(h=alpha,lty=3)
}
par(op)
dev.off()
