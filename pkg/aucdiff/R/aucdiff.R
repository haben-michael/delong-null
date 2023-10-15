## ## ## parametric formulas for homoskedastic guassian covariates with lda coefs
## ## coefs.lda <- function(x.0,x.1,params=NULL) {
## ##     ## x.0 <- x.0[,-1]; x.1 <- x.1[,-1]
## ##     ## browser()
## ##     n.0 <- nrow(x.0); n.1 <- nrow(x.1); n <- n.0+n.1
## ##     mu.0.hat <- colMeans(x.0);  mu.1.hat <- colMeans(x.1)
## ##     Sigma <- params$Sigma
## ##     if(is.null(Sigma)) Sigma <- (  with(list(x.scaled=scale(x.0,scale=FALSE)), t(x.scaled)%*%x.scaled) + with(list(x.scaled=scale(x.1,scale=FALSE)), t(x.scaled)%*%x.scaled)  ) / n
## ##     pi.0 <- params$pi.0 ## for heteroskedasticity
## ##     ## browser()
## ##     Sigma.0 <- params$Sigma.0; Sigma.1 <- params$Sigma.1
## ##     if(is.null(pi.0))pi.0 <- 1/2 
## ##     if(is.null(Sigma.0))Sigma.0 <- Sigma
## ##     if(is.null(Sigma.1))Sigma.1 <- Sigma
## ##     Sigma <- pi.0*Sigma.0 + (1-pi.0)*Sigma.1
## ##     Sigma.inv <- solve(Sigma)
## ##     ## if(max(abs(Sigma-ss))>.003)browser()
## ##     beta.hat <- Sigma.inv%*%(mu.1.hat-mu.0.hat)
## ##     return(beta.hat)
## ## }
## hajek.lda <- function(x,y,beta,params) {
##     mu.0 <- params$mu.0; mu.1 <- params$mu.1; Sigma <- params$Sigma
##     qform <- as.numeric(t(beta)%*%Sigma%*%beta)
##     x <- t(x); y <- t(y) ## expects x,y in model matrix format
##     mean(pnorm(t(beta)%*%(y-mu.0)/sqrt(qform))) + mean(pnorm(t(beta)%*%(mu.1-x)/sqrt(qform))) - 2*pnorm(t(beta)%*%(mu.1-mu.0)/sqrt(2*qform))
## }
## hajek.lda.deriv <- function(x,y,beta,params) {
##     mu.0 <- params$mu.0; mu.1 <- params$mu.1; Sigma <- params$Sigma
##     qform <- as.numeric(t(beta)%*%Sigma%*%beta)
##     f.prime.pre <- function(u) {
##         colMeans(as.numeric(dnorm(t(beta)%*%u/sqrt(qform))) * (t(u/sqrt(qform)) - t(t(beta)%*%u/qform^(3/2))%*%t(beta)%*%Sigma))
##     }
##     f.prime.pre(t(y)-mu.0) + f.prime.pre(-t(x)+mu.1) -  2*f.prime.pre(matrix((mu.1-mu.0)/sqrt(2),ncol=1))
## }
## hajek.lda.deriv.var <- function(beta,params) {
##     u.1 <- function(beta,w) {
##         qform <- as.numeric(t(beta)%*%Sigma%*%beta)
##         as.numeric(-dnorm(t(beta)%*%w/sqrt(qform))/sqrt(qform))*(as.numeric(t(beta)%*%w/qform)*Sigma%*%beta-w)
##     }
##     Sigma <- params$Sigma
##     mu <- with(params, mu.1-mu.0)
##     qform <- as.numeric(t(beta)%*%Sigma%*%beta)
##     moment.1 <- u.1(beta,mu/sqrt(2))
##     moment.2 <- 1/sqrt(2*pi)*as.numeric(dnorm(sqrt(2/3)*t(beta)%*%mu/sqrt(qform))/sqrt(3)/qform) * (
##         mu%*%t(mu) + Sigma - as.numeric(t(beta)%*%mu/qform) * (mu%*%t(beta)%*%Sigma + Sigma%*%beta%*%t(mu)) + as.numeric((t(beta)%*%mu/qform)^2-1/qform)*Sigma%*%beta%*%t(beta)%*%Sigma
##     )
##     moment.2 - moment.1%*%t(moment.1)
## }
## vcov.full.hajek.lda <- function(beta,params) {
##     mu <- with(params, mu.1-mu.0); Sigma <- params$Sigma
##     u.1 <- function(beta,w) {
##         qform <- as.numeric(t(beta)%*%Sigma%*%beta)
##         as.numeric(-dnorm(t(beta)%*%w/sqrt(qform))/sqrt(qform))*(as.numeric(t(beta)%*%w/qform)*Sigma%*%beta-w)
##     }
##     q <- as.numeric(t(beta)%*%Sigma%*%beta)
##     A.1 <- mu - as.numeric(t(beta)%*%mu/2)*Sigma%*%beta/q
##     A.2 <- t(mu)/sqrt(2) - as.numeric(t(beta)%*%mu*2^(-3/2))*t(beta)%*%Sigma/q
##     A <- as.numeric(dnorm(t(beta)%*%mu/sqrt(2*q))) * ( A.1%*%A.2 + Sigma/sqrt(2) - (Sigma%*%beta)%*%t(Sigma%*%beta)/q/2^(3/2))
##     mean.try <- 1/sqrt(q)*(diag(p) - 1/q*Sigma%*%beta%*%t(beta)) %*%A
##     cov.try <- solve(Sigma)%*%(mean.try - u.1(beta,mu/sqrt(2))%*%t(mu)) 
## }
## vcov.lda <- function(beta,p.reduced,params) {
##     hajek.lda.deriv.var <- function(beta,params) {
##         u.1 <- function(beta,w) {
##             qform <- as.numeric(t(beta)%*%Sigma%*%beta)
##             as.numeric(-dnorm(t(beta)%*%w/sqrt(qform))/sqrt(qform))*(as.numeric(t(beta)%*%w/qform)*Sigma%*%beta-w)
##         }
##         Sigma <- params$Sigma
##         mu <- with(params, mu.1-mu.0)
##         qform <- as.numeric(t(beta)%*%Sigma%*%beta)
##         moment.1 <- u.1(beta,mu/sqrt(2))
##         moment.2 <- 1/sqrt(2*pi)*as.numeric(dnorm(sqrt(2/3)*t(beta)%*%mu/sqrt(qform))/sqrt(3)/qform) * (
##             mu%*%t(mu) + Sigma - as.numeric(t(beta)%*%mu/qform) * (mu%*%t(beta)%*%Sigma + Sigma%*%beta%*%t(mu)) + as.numeric((t(beta)%*%mu/qform)^2-1/qform)*Sigma%*%beta%*%t(beta)%*%Sigma
##         )
##         moment.2 - moment.1%*%t(moment.1)
##     }
##     vcov.full.hajek.lda <- function(beta,params) {
##         mu <- with(params, mu.1-mu.0); Sigma <- params$Sigma
##         u.1 <- function(beta,w) {
##             qform <- as.numeric(t(beta)%*%Sigma%*%beta)
##             as.numeric(-dnorm(t(beta)%*%w/sqrt(qform))/sqrt(qform))*(as.numeric(t(beta)%*%w/qform)*Sigma%*%beta-w)
##         }
##         q <- as.numeric(t(beta)%*%Sigma%*%beta)
##         A.1 <- mu - as.numeric(t(beta)%*%mu/2)*Sigma%*%beta/q
##         A.2 <- t(mu)/sqrt(2) - as.numeric(t(beta)%*%mu*2^(-3/2))*t(beta)%*%Sigma/q
##         A <- as.numeric(dnorm(t(beta)%*%mu/sqrt(2*q))) * ( A.1%*%A.2 + Sigma/sqrt(2) - (Sigma%*%beta)%*%t(Sigma%*%beta)/q/2^(3/2))
##         mean.try <- 1/sqrt(q)*(diag(p) - 1/q*Sigma%*%beta%*%t(beta)) %*%A
##         cov.try <- solve(Sigma)%*%(mean.try - u.1(beta,mu/sqrt(2))%*%t(mu)) 
##     }
##     p.red <- p.reduced
##     cov.beta.hajek <- matrix(0,p,p)
##     cov.beta.hajek[1:p.red,1:p.red] <- vcov.full.hajek.lda(beta,params)[1:p.red,1:p.red] #* (1/m+1/n)
##     cov.gamma.hajek <- vcov.full.hajek.lda(beta,params)# * (1/m+1/n)
##     var.gamma <- solve(Sigma) #* (1/m+1/n)
##     var.beta <- matrix(0,p,p)
##     var.beta[1:p.red,1:p.red] <- solve(Sigma[1:p.red,1:p.red]) #* (1/m+1/n)
##     cov.beta.gamma <- solve(Sigma[1:p.red,1:p.red]) %*% Sigma[1:p.red,] %*% solve(Sigma) #* (1/m+1/n)
##     cov.beta.gamma <- rbind(cov.beta.gamma,matrix(0,p.red,p))
##     vcov.beta.gamma <- rbind(cbind(var.beta,cov.beta.gamma),cbind(t(cov.beta.gamma),var.gamma))
##     var.hajek <- hajek.lda.deriv.var(beta,params)#*(1/m+1/n)
##     vcov.beta.gamma.hajek <- rbind(cbind(vcov.beta.gamma,rbind(cov.beta.hajek,cov.gamma.hajek)), cbind(t(cov.beta.hajek),t(cov.gamma.hajek),var.hajek)) 
## }







## ############### mrc with normal-laplace model



#' Estimate the AUC
#' @export
auc.hat <- function(x,y,ties=FALSE) {
    stopifnot(ties==FALSE)
    w <- c(x,y)
    m <- length(x); n <- length(y)
    d <- rep(0:1,c(m,n))
    idx <- order(w)
    d <- d[idx]
    ## w <- w[idx]
    sum(cumsum(1-d)[d==1]) / (m*n)
}


## ## coef.mrc <- function(x,y,plot=FALSE) {
## ##     ## browser()
## ##     knots <- sort(as.numeric( - outer(y[,1],x[,1],`-`) / outer(y[,2],x[,2],`-`)))
## ##     midpoints <- (knots[-1]+knots[-length(knots)])/2
## ##     auc.hats <- sapply(midpoints, function(beta2)auc.hat(x%*%c(1,beta2),y%*%c(1,beta2)))
## ##     beta.hat <-     midpoints[which.max(auc.hats)]
## ##     if(plot){
## ##         f <- stepfun(knots[-c(1,length(knots))],auc.hats)
## ##         plot(f,do.points=FALSE)
## ##         abline(v=beta.hat)
## ##     }
## ##     return(beta.hat)
## ## }
## coef.mrc <- function(x,y,plot=FALSE) {
##     idx <- expand.grid(seq(nrow(x)),seq(nrow(y)))
##     paired <- cbind(x[idx[,1],1] - y[idx[,2],1],x[idx[,1],2] - y[idx[,2],2])
##     increasing <- paired[paired[,2]<0,]
##     decreasing <- paired[paired[,2]>0,]
##     knots.inc <- as.matrix(-increasing[,1] / increasing[,2])
##     knots.dec <- as.matrix(-decreasing[,1] / decreasing[,2])
##     knots <- rbind(cbind(knots.inc,1),cbind(knots.dec,-1))
##     knots <- knots[order(knots[,1]),]
##     vals <- cumsum(knots[,2])
##     beta.hat <-     knots[which.max(vals),1]
##     if(plot){
##         plot(knots[,1],(vals-min(vals))/diff(range(vals)),col=2,type='l')
##         abline(v=beta.hat)
##     }
##     return(beta.hat)
## }


## coef.mrc2 <- function(x,y,epsilon=1e-4) {
##     idx <- expand.grid(1:nrow(x),1:nrow(y))
##     uvw <- x[idx[,1],] - y[idx[,2],]

##     idx <- expand.grid(rep(list(1:nrow(uvw)),2))
##     idx <- idx[idx[,1] < idx[,2],]
##     uvw.pairs <- cbind(uvw[idx[,1],],uvw[idx[,2],])
##     uvw.1 <- uvw.pairs[,1:3]; uvw.2 <- uvw.pairs[,4:6]
##     betas <- cbind( uvw.1[,1]*uvw.2[,3] - uvw.2[,1]*uvw.1[,3], uvw.1[,2]*uvw.2[,1] - uvw.1[,1]*uvw.2[,2]) / (uvw.2[,2]*uvw.1[,3] - uvw.1[,2]*uvw.2[,3])
##     betas[is.infinite(betas)] <- NA
##     vals.beta <- apply(betas, 1, function(beta) mean(uvw%*%c(1,beta[1],beta[2])<0))

##     ## (cont) look in nbhds of maximizing knots for maximum
##     max.idx <- which(vals.beta==max(vals.beta,na.rm=TRUE))
##     origins <- betas[max.idx,,drop=FALSE]
##     ray.1 <- uvw.1[max.idx,3:2,drop=FALSE]
##     ray.2 <- uvw.2[max.idx,3:2,drop=FALSE]
##     ray.1[,1] <- -ray.1[,1]; ray.2[,1] <- -ray.2[,1]
##     ray.1 <- t(apply(ray.1,1,function(x)x/sqrt(sum(x^2))))
##     ray.2 <- t(apply(ray.2,1,function(x)x/sqrt(sum(x^2))))
##     max.knot.vals <- mapply(split(origins,1:nrow(origins)), split(ray.1,1:nrow(ray.1)), split(ray.2,1:nrow(ray.2)),FUN=function(origin.i,ray.1i,ray.2i) {
##         ## browser()
##         test.vals <- apply(expand.grid(c(-1,1),c(-1,1)), 1, function(signs) {
##             test.beta <- origin.i + epsilon*(signs[1]*ray.1i+signs[2]*ray.2i)
##             c(beta2=test.beta[1],beta2=test.beta[2], val=mean(uvw%*%c(1,test.beta)<0))
##         })#,simplify=FALSE)
##         test.vals <- cbind(test.vals,c(origin.i,val=mean(uvw%*%c(1,origin.i)<0)))
##         return(test.vals[,which.max(test.vals['val',])])
##     })
##     beta.hat <- max.knot.vals[1:2,which.max(max.knot.vals['val',])]
##     return(beta.hat)
## }


## coef.mrc.shin <- function(x,y,lim=1e5) {
##     n <- length(y)
##     d <- c(rep(0,nrow(x)),rep(1,nrow(y)))
##     x <- rbind(x,y)
##     y <- d
##     idx <- expand.grid(seq(nrow(x)),seq(nrow(x)))
##     idx <- idx[idx[,1]<idx[,2],]
##     x.ij <- x[idx[,1],] - x[idx[,2],]
##     y.ij <- y[idx[,1]] - y[idx[,2]]
##     obj <- c(rep(0,p-1), as.numeric(y.ij > 0)) / (n*(n-1))
##     M.ij <- rep(lim,nrow(x.ij))
##     A.1  <- A.2 <- cbind(x.ij[,-1], diag(-M.ij))
##     rhs.1 <- -M.ij-x.ij[,1]
##     sense.1 <- rep('>',nrow(x.ij))
##     rhs.2 <- -x.ij[,1]
##     sense.2 <- rep('<',nrow(x.ij))
##     vtype <- c(rep('C',p-1),rep('B',nrow(x.ij)))
##     model <- within(list(), {modelsense <- 'max'; obj <- obj; A <- rbind(A.1,A.2); rhs <- c(rhs.1,rhs.2); sense <- c(sense.1,sense.2); vtype <- vtype}) # clean
##     params <- list(OutputFlag=0)
##     result <- gurobi::gurobi(model,params)
##     return(result$x[1:(p-1)])
## }



## pnl <- function(q,mean,sd,rate.1,rate.2=rate.1) {
##     mills <- function(x)ifelse(is.infinite(x),0,(1-pnorm(x))/dnorm(x))
##     z <- (q-mean)/sd
##     pnorm(z) - dnorm(z)*(mills(rate.1*sd-z)-rate.1/rate.2*mills(rate.2*sd+z))/(1+rate.1/rate.2)
## }

## ## todo: switch a,b to params arguments or vice versa in hajek functions
## auc.nl <- function(beta,a,b)1-pnl(0,beta[2]+a*beta[1],sqrt(beta[2]^2+beta[3]^2),1/b/beta[1])

## auc.grad.nl <- function(beta,a,b) {
##     mills <- function(x)ifelse(is.infinite(x),0,(1-pnorm(x))/dnorm(x))
##     mills.prime <- function(x) x*(1-pnorm(x))/dnorm(x)-1
##     sigma <- sqrt(beta[2]^2+beta[3]^2)
##     z <- -(beta[2]+beta[1]*a)/sigma
##     alpha.sigma <- sigma/beta[1]/b
##     z.grad <- -c(a/sigma,(sigma-(beta[2]+beta[1]*a)*beta[2]/sigma)/sigma^2,-beta[3]*(beta[2]+beta[1]*a)/sigma^3)
##     alpha.sigma.grad <- 1/b*c(-sigma^2/beta[1],beta[2],beta[3])/beta[1]/sigma
##     out <- -dnorm(z)*z.grad - z*dnorm(z)*z.grad*1/2*(mills(alpha.sigma-z)-mills(alpha.sigma+z)) + 1/2*dnorm(z)*(mills.prime(alpha.sigma-z)*(alpha.sigma.grad-z.grad)-mills.prime(alpha.sigma+z)*(alpha.sigma.grad+z.grad))
##     return(out)
## }

## auc.hessian.nl <- function(beta,a,b) {
##     sigma <- sqrt(beta[2]^2+beta[3]^2)
##     sigma.grad <- c(0,beta[2:3])/sigma
##     z <- -(beta[2]+beta[1]*a)/sigma
##     z.grad <- -c(a/sigma,1/sigma-beta[2]*(beta[2]+beta[1]*a)/sigma^3,-beta[3]*(beta[2]+beta[1]*a)/sigma^3)
##     z.hessian <- -rbind(-a/sigma^2*sigma.grad,
##                         -1/sigma^2*sigma.grad + 3*beta[2]*(beta[2]+beta[1]*a)/sigma^4*sigma.grad - c(a*beta[2],2*beta[2]+a*beta[1],0)/sigma^3,
##                         3*beta[3]*(beta[2]+beta[1]*a)/sigma^4*sigma.grad-c(a*beta[3],beta[3],beta[2]+beta[1]*a)/sigma^3)
##     mills <- function(x)ifelse(is.infinite(x),0,(1-pnorm(x))/dnorm(x))
##     mills.prime <- function(x) x*(1-pnorm(x))/dnorm(x)-1
##     mills.prime.prime <- function(x)mills(x)*(1+x^2)-x
##     alpha.sigma <- sigma/beta[1]/b
##     alpha.sigma.grad <- 1/b*c(-sigma^2/beta[1],beta[2],beta[3])/beta[1]/sigma
##     alpha.sigma.hessian <- 1/b*rbind(-sigma*c(-2/beta[1]^3,0,0)-1/beta[1]^2*sigma.grad, 1/sigma*c(-beta[2]/beta[1]^2,1/beta[1],0)-beta[2]/beta[1]/sigma^2*sigma.grad,1/sigma*c(-beta[3]/beta[1]^2,0,1/beta[1])-beta[3]/beta[1]/sigma^2*sigma.grad  )
##     mills.term <- mills(alpha.sigma-z)-mills(alpha.sigma+z)
##     mills.term.grad <- mills.prime(alpha.sigma-z)*(alpha.sigma.grad-z.grad)-mills.prime(alpha.sigma+z)*(alpha.sigma.grad+z.grad)
##     mills.term.hessian <- mills.prime.prime(alpha.sigma-z)*outer(alpha.sigma.grad-z.grad,alpha.sigma.grad-z.grad) - mills.prime.prime(alpha.sigma+z)*outer(alpha.sigma.grad+z.grad,alpha.sigma.grad+z.grad) + mills.prime(alpha.sigma-z)*(alpha.sigma.hessian-z.hessian) - mills.prime(alpha.sigma+z)*(alpha.sigma.hessian+z.hessian)
##     out <- z*dnorm(z)*outer(z.grad,z.grad) - dnorm(z)*z.hessian - 1/2*mills.term*(dnorm(z)*outer(z.grad,z.grad) - z^2*dnorm(z)*outer(z.grad,z.grad) + z*dnorm(z)*z.hessian) - 1/2*z*dnorm(z)*outer(z.grad,mills.term.grad) - 1/2*z*dnorm(z)*outer(mills.term.grad,z.grad) + 1/2*dnorm(z)*mills.term.hessian
##     return(out)
##     }


## hajek.nl <- function(x,y,beta,params) {
##     a <- params$a; b <- params$b
##     auc.beta <- auc.nl(beta,a,b)
##     F <- function(q,beta,a,b) pnl(q, -beta[2],sd=sqrt(1/2*(beta[2]^2+beta[3]^2)),1/beta[1]/b,Inf)
##     G <- function(q,beta,a,b) pnl(q, a*beta[1],sd=sqrt(1/2*(beta[2]^2+beta[3]^2)),1/beta[1]/b,Inf)
##     mean(F(y%*%beta,beta,a,b)) + mean(1-G(x%*%beta,beta,a,b)) - 2*auc.beta
## }



## hajek.grad.nl <- function(x,y,beta,params) {
##     ## browser()
##     a <- params$a; b <- params$b
##     F.grad <- function(y,beta,a,b) {
##         sigma <- sqrt(beta[2]^2+beta[3]^2) 
##         z  <- sqrt(2)*(y%*%beta+beta[2])/sigma
##         z.grad <- sqrt(2)/sigma*c(y[1],1+y[2]-beta[2]*(y%*%beta+beta[2])/sigma^2, y[3] - beta[3]*(y%*%beta+beta[2])/sigma^2)
##         mills <- function(x)ifelse(is.infinite(x),0,(1-pnorm(x))/dnorm(x))
##         mills.prime <- function(x) x*(1-pnorm(x))/dnorm(x)-1
##         mills.arg <- sqrt(sigma^2)/(sqrt(2)*beta[1]*b) - z
##         mills.arg.grad <-  c(-sqrt(sigma^2)/sqrt(2)/b/beta[1]^2,beta[2]/b/beta[1]/sqrt(2*sigma^2),beta[3]/b/beta[1]/sqrt(2*sigma^2)) - z.grad  
##         dnorm(z)*z.grad + z*dnorm(z)*z.grad*mills(mills.arg) - dnorm(z)*mills.prime(mills.arg)*mills.arg.grad
##     }
##     ## todo: combine redundant code in F.grad,G.grad
##     G.grad <- function(x,beta,a,b) {
##         sigma <- sqrt(beta[2]^2+beta[3]^2) 
##         z  <- sqrt(2)*(x%*%beta-a*beta[1])/sigma
##         z.grad <- sqrt(2)/sigma*c(x[1]-a,x[2]-beta[2]*(x%*%beta-a*beta[1])/sigma^2, x[3] - beta[3]*(x%*%beta-a*beta[1])/sigma^2)
##         mills <- function(x)ifelse(is.infinite(x),0,(1-pnorm(x))/dnorm(x))
##         mills.prime <- function(x) x*(1-pnorm(x))/dnorm(x)-1
##         mills.arg <- sqrt(sigma^2)/(sqrt(2)*beta[1]*b) - z
##         mills.arg.grad <-  c(-sqrt(sigma^2)/sqrt(2)/b/beta[1]^2,beta[2]/b/beta[1]/sqrt(2*sigma^2),beta[3]/b/beta[1]/sqrt(2*sigma^2)) - z.grad  
##         dnorm(z)*z.grad + z*dnorm(z)*z.grad*mills(mills.arg) - dnorm(z)*mills.prime(mills.arg)*mills.arg.grad
##     }
##     ## todo: vectorize F.grad, G.grad
##     F.y <- apply(y,1,function(y.i)F.grad(y.i,beta,a,b))
##     G.x <- apply(x,1,function(x.i)G.grad(x.i,beta,a,b))
##     auc.grad.beta <- auc.grad.nl(beta,a,b)
##     rowMeans(F.y) - rowMeans(G.x) - 2*auc.grad.beta
## }




## parametric routines for gaussian covariates with lda coefs
## todo: x,y used here, x,d used in delong-alt. rename hajek to hoeffding?

#' @export
hajek.normal <- function(x,y,beta,params) {
    mu.0 <- params$mu.0; mu.1 <- params$mu.1
    Sigma.0 <- params$Sigma.0; Sigma.1 <- params$Sigma.1
    if(is.null(Sigma.0)||is.null(Sigma.1)) Sigma.0 <- Sigma.1 <- params$Sigma
    Sigma <- (Sigma.0+Sigma.1)/2
    qform <- function(Sigma)as.numeric(t(beta)%*%Sigma%*%beta)
    x <- t(x); y <- t(y) ## expects x,y in model matrix format
    mean(pnorm(t(beta)%*%(y-mu.0)/sqrt(qform(Sigma.0)))) + mean(pnorm(t(beta)%*%(mu.1-x)/sqrt(qform(Sigma.1)))) - 2*pnorm(t(beta)%*%(mu.1-mu.0)/sqrt(2*qform(Sigma)))
}

#' @export
hajek.grad.normal <- function(x,y,beta,params) {
    mu.0 <- c(params$mu.0); mu.1 <- c(params$mu.1)
    Sigma.0 <- params$Sigma.0; Sigma.1 <- params$Sigma.1
    if(is.null(Sigma.0)||is.null(Sigma.1)) Sigma.0 <- Sigma.1 <- params$Sigma
    f.prime.pre <- function(u,Sigma) {
        qform <- as.numeric(t(beta)%*%Sigma%*%beta)
        colMeans(as.numeric(dnorm(t(beta)%*%u/sqrt(qform))) * (t(u/sqrt(qform)) - t(t(beta)%*%u/qform^(3/2))%*%t(beta)%*%Sigma))
    }
    ## browser()
    f.prime.pre(t(y)-mu.0,Sigma.0) + f.prime.pre(-t(x)+mu.1,Sigma.1) -  2*f.prime.pre(matrix((mu.1-mu.0),ncol=1),Sigma.0+Sigma.1)
}

u.1 <- function(beta,Sigma,w) {## f.prime.pre is a vectorized version; use that
    qform <- as.numeric(t(beta)%*%Sigma%*%beta)
    as.numeric(-dnorm(t(beta)%*%w/sqrt(qform))/sqrt(qform))*(as.numeric(t(beta)%*%w/qform)*Sigma%*%beta-w)
}

J.1 <- function(beta,mu,Sigma) {
    quad <- as.numeric(t(beta)%*%Sigma%*%beta)
    A <- mu-Sigma%*%beta%*%t(beta)%*%mu/(1+quad)
    B <- t(mu)/sqrt(1+quad)-(t(beta)%*%mu*(1+quad)^(-3/2))%*%(t(beta)%*%Sigma)
    C <- Sigma/sqrt(1+quad) - Sigma%*%beta%*%t(beta)%*%Sigma*(1+quad)^(-3/2)
    as.numeric(dnorm(t(beta)%*%mu/sqrt(1+quad))) * (A%*%B+C)
}
J.2 <- function(beta,mu,Sigma.0,Sigma.1) {
    p <- length(beta) ## patched without debugging
    qform <- as.numeric(t(beta)%*%Sigma.1%*%beta)
    (sqrt(qform)*diag(p)-Sigma.1%*%beta%*%t(beta)/sqrt(qform))%*%J.1(beta,mu/sqrt(qform),Sigma.0/qform)-u.1(beta,Sigma.0+Sigma.1,mu/sqrt(2))%*%t(mu)
}
J.3 <- function(beta,mu,Sigma.0,Sigma.1) {
    q.1 <- as.numeric(t(beta)%*%Sigma.1%*%beta)
    J <- J.1(beta,mu/sqrt(q.1/2),2*Sigma.0/q.1)
    (J/2 + 1/2/q.1^2*as.numeric(t(beta)%*%J%*%beta)*Sigma.1%*%beta%*%t(beta)%*%Sigma.1 - 1/2/q.1*(Sigma.1%*%beta%*%t(beta)%*%J + J%*%beta%*%t(beta)%*%Sigma.1)) / sqrt(2*pi)    
}
J.4 <- function(beta,mu,Sigma,p.r) {
    Sigma.f <- Sigma
    Sigma.r <- Sigma.f[1:p.r,1:p.r]
    mu.f <- mu; mu.r <- mu.f[1:p.r]
    beta.f <- beta; beta.r <- beta.f[1:p.r]
    q <- as.numeric(t(beta.f)%*%Sigma.f%*%beta.f)
    A <- 1/2/q*Sigma.r%*%beta.r%*%beta.f%*%Sigma.f * as.numeric(t(beta.f)%*%J.1(beta.f,sqrt(2/q)*mu.f,2/q*Sigma.f)%*%beta.f)
    B <- 1/2*Sigma.r%*%beta.r%*%beta.f%*%J.1(beta.f,sqrt(2/q)*mu.f,2/q*Sigma.f)
    C <- 1/2*J.1(beta.r,sqrt(2/q)*mu.r,2/q*Sigma.r)%*%beta.r%*%t(beta.f)%*%Sigma.f
    D <- q/2*J.1(beta.f,sqrt(2/q)*mu.f,2/q*Sigma.f)[1:p.r,]
    (A-B-C+D)/q/sqrt(2*pi)
}


## update 8/30: m,n no longer extern, must be included in params
## update 10/1: m,n no longer used, use pi.0 instead. return value
## must be normalized by m+n

#' @export
var.hajek.normal <- function(beta,params) {  
    ## browser()
    ## N <- params$N
    pi.0 <- params$pi.0; pi.1 <- 1-pi.0
    ## m <- params$m; n <- params$n
    mu.0 <- params$mu.0; mu.1 <- params$mu.1
    Sigma.0 <- params$Sigma.0; Sigma.1 <- params$Sigma.1
    if(is.null(Sigma.0)||is.null(Sigma.1)) Sigma.0 <- Sigma.1 <- params$Sigma
    ## mu <- mu.1 - mu.0
    moment.1 <- u.1(beta,Sigma.0+Sigma.1,mu.1-mu.0)
    ## 1/n*J.3(beta,mu.1-mu.0,Sigma.1,Sigma.0)+1/m*J.3(beta,mu.1-mu.0,Sigma.0,Sigma.1) - moment.1%*%t(moment.1)*(1/m+1/n)
    1/pi.1*J.3(beta,mu.1-mu.0,Sigma.1,Sigma.0)+1/pi.0*J.3(beta,mu.1-mu.0,Sigma.0,Sigma.1) - moment.1%*%t(moment.1)*(1/pi.0+1/pi.1)

    ## by.status <- lapply(list(control=Sigma.0,case=Sigma.1), function(Sigma) {
    ##     q <- as.numeric(t(beta)%*%Sigma%*%beta)    
    ##     J <- J.1(beta,mu/sqrt(q/2),2*Sigma/q)
    ##     moment.2 <- (J/2 + 1/2/q^2*as.numeric(t(beta)%*%J%*%beta)*Sigma%*%beta%*%t(beta)%*%Sigma - 1/2/q*(Sigma%*%beta%*%t(beta)%*%J + J%*%beta%*%t(beta)%*%Sigma)) / sqrt(2*pi)
    ##     moment.1 <- u.1(beta,Sigma,mu/sqrt(2))
    ##     moment.2 - moment.1%*%t(moment.1)
    ## })
    ## ## return(by.status[['control']]/N/pi.0 + by.status[['case']]/N/pi.1)
    ## return(by.status[['control']]/m + by.status[['case']]/n)
}

var.hajek.grad.normal <- var.hajek.normal #deprecate

## ## deprecate in favor of auc.hessian.normal--second derivative of auc(beta), but only valid when evaluated at beta==the lda coefs
## auc.scores.deriv2.lda <- function(coefs,Sigma.diff) {
##       quad <- as.numeric(t(coefs)%*%Sigma.diff%*%coefs)
##       -dnorm(sqrt(quad)/2)/sqrt(quad)*(Sigma.diff - (Sigma.diff%*%coefs)%*%t(Sigma.diff%*%coefs)/quad)/2
##       }

#' @export
auc.normal <- function(beta,params) {
    mu <- with(params,mu.1-mu.0)
    Sigma <- with(params,Sigma.0+Sigma.1)
    q <- as.numeric(t(beta)%*%Sigma%*%beta)
    ## with(with(params,list(mu=mu.1-mu.0,Sigma=Sigma.0+Sigma.1)),
    pnorm(t(beta)%*%mu / sqrt(q))
    }

#' @export
auc.grad.normal <- function(beta,params) {
    mu <- with(params,mu.1-mu.0)
    Sigma <- with(params,Sigma.0+Sigma.1)
    q <- as.numeric(t(beta)%*%Sigma%*%beta)
    q.grad  <- 2*Sigma%*%beta
    ## with(with(params,list(mu=mu.1-mu.0,Sigma=Sigma.0+Sigma.1)),
    c(dnorm(t(beta)%*%mu/sqrt(q)))*(mu/sqrt(q)-1/2*c(t(beta)%*%mu)/q^(3/2)*q.grad)
    }

#' @export
auc.hessian.normal <- function(beta,params) {
    mu <- with(params,mu.1-mu.0)
    Sigma <- with(params,Sigma.0+Sigma.1)
    q <- as.numeric(t(beta)%*%Sigma%*%beta)
    q.grad  <- 2*Sigma%*%beta
    ## with(with(params,list(mu=mu.1-mu.0,Sigma=Sigma.0+Sigma.1)),
    -c(t(beta)%*%mu/sqrt(q))*auc.grad.normal(beta,params)%*%t(auc.grad.normal(beta,params))/c(dnorm(t(beta)%*%mu/sqrt(q))) + c(dnorm(t(beta)%*%mu/sqrt(q))) * ( -1/2/q^(3/2)* mu%*%t(q.grad) - 1/2*q.grad%*%t(mu/q^(3/2)-3/2*c(t(beta)%*%mu/q^(5/2))*q.grad) - 1/2*c(t(beta)%*%mu/q^(3/2))*2*Sigma )
}


## routines specific to lda betahat with nested models

#' Compute LDA coefficients
#'
#' @param x.0 control observations
#' @param x.1 case observations
#' @param params parameters of the normal model, which are estimated
#'     if NULL value is passed
#' @return coefficient estimates
#' @export
coefs.lda <- function(x.0,x.1,params=NULL) {
    ## x.0 <- x.0[,-1]; x.1 <- x.1[,-1]
    n.0 <- nrow(x.0); n.1 <- nrow(x.1); n <- n.0+n.1
    mu.0.hat <- colMeans(x.0);  mu.1.hat <- colMeans(x.1)
    Sigma <- params$Sigma
    if(is.null(Sigma)) Sigma <- (  with(list(x.scaled=scale(x.0,scale=FALSE)), t(x.scaled)%*%x.scaled) + with(list(x.scaled=scale(x.1,scale=FALSE)), t(x.scaled)%*%x.scaled)  ) / n
    pi.0 <- params$pi.0 ## for heteroskedasticity
    ## browser()
    Sigma.0 <- params$Sigma.0; Sigma.1 <- params$Sigma.1
    if(is.null(pi.0))pi.0 <- 1/2 
    if(is.null(Sigma.0))Sigma.0 <- Sigma
    if(is.null(Sigma.1))Sigma.1 <- Sigma
    ## browser()
    Sigma <- pi.0*Sigma.0 + (1-pi.0)*Sigma.1
    Sigma.inv <- solve(Sigma)
    ## if(max(abs(Sigma-ss))>.003)browser()
    beta.hat <- Sigma.inv%*%(mu.1.hat-mu.0.hat)
    return(beta.hat)
}

## ## based on delong_alt/utils.R. expects x in model matrix format. making
## ## hetero/homoskedasticity explicit. adding empirical estimates to
## ## allow for params=NULL.
## infl.lda <- function(w,d,params=NULL,var.equal=TRUE,terms.only=TRUE) {
##     ## browser()
##     g <- d # clean up
##     n <- length(d) # clean up redundancies here
##     n.1 <- sum(d); n.0 <- n-n.1
##     mu.0 <- params$mu.0
##     mu.1 <- params$mu.1
##     pi.0 <- params$pi.0
##     x.0 <- w[d==0,]; x.1 <- w[d==1,]
##     if(is.null(mu.0))mu.0 <- colMeans(x.0)
##     if(is.null(mu.1))mu.1 <- colMeans(x.1)
##     if(is.null(pi.0))pi.0 <- n.0/n#nrow(x.0)/(nrow(x.0)+nrow(x.1))
##     pi.1 <- 1-pi.0
##     if(var.equal) {
##         Sigma <- params$Sigma
##         if(is.null(Sigma))Sigma <- (var(x.0)*n.0+var(x.1)*n.1)/n
##     } else {
##         Sigma.0 <- params$Sigma.0; Sigma.1 <- params$Sigma.1
##         if(is.null(Sigma.0))Sigma.0 <- var(x.0)
##         if(is.null(Sigma.1))Sigma.1 <- var(x.1)
##         Sigma <- pi.0*Sigma.0 + pi.1*Sigma.1
##     }
##     w <- t(w)#[-1,]
##     ## x.0 <- t(x[,g==0]); x.1 <- t(x[,g==1])
##     d <- sort(d)
##     ## mu.0.hat <- rowMeans(x.0);  mu.1.hat <- rowMeans(x.1)
##     ## pi.1.hat <- n.1/(n.0+n.1)
##     ## infl.0 <-   g/pi.1/(1-pi.1)-1/(1-pi.1) - t(mu.1)%*%solve(Sigma)%*%(x-mu.1) * n/n.1 * g  + t(mu.0)%*%solve(Sigma)%*%(x-mu.0) * n/n.0 * (1-g)
##     infl <- t(d*t(solve(Sigma)%*%(w-mu.1)*n/n.1)  -  (1-d)*t(solve(Sigma)%*%(w-mu.0)*n/n.0))
##     ## infl <- rbind(infl.0,infl)
##     if(terms.only) return(infl) else return(rowMeans(infl))
## }


## expects w in column major format. differs from the delong_alt
## version as it allows for unknown Sigmas.
## 10-3: w in model matrix format


#' Influence function for LDA estimates
#'
#' @param w covariates
#' @param d class labels, 0 or 1
#' @param params parameters of the normal model, which are estimated
#'     if NULL value is passed
#' @param terms.only boolean, whether to return the terms of the
#'     influence function directly or a class-weighted average
#' @return depending on "terms.only", either the terms of the
#'     influence function directly or a class-weighted average
#' @export
infl.lda <- function(w,d,params=NULL,terms.only=TRUE) {
    ## g <- d # clean up
    ## browser()
    n <- length(d) # clean up redundancies here
    p <- ncol(w)
    n.1 <- sum(d); n.0 <- n-n.1
    mu.0 <- c(params$mu.0)
    mu.1 <- c(params$mu.1)
    pi.0 <- params$pi.0
    x.0 <- w[d==0,]; x.1 <- w[d==1,]
    Sigma.0 <- params$Sigma.0; Sigma.1 <- params$Sigma.1
    if(is.null(mu.0))mu.0 <- c(colMeans(x.0))
    if(is.null(mu.1))mu.1 <- c(colMeans(x.1))
    if(is.null(pi.0))pi.0 <- n.0/n#nrow(x.0)/(nrow(x.0)+nrow(x.1))
    ## browser()
    pi.1 <- 1-pi.0
    if(is.null(Sigma.0))Sigma.0 <- var(x.0)
    if(is.null(Sigma.1))Sigma.1 <- var(x.1)
    Sigma <- pi.0*Sigma.0 + pi.1*Sigma.1
    x.1 <- t(x.1); x.0 <- t(x.0)
    terms.case <-  solve(Sigma)%*%x.1 - t(as.numeric(t(x.1-mu.1)%*%solve(Sigma)%*%(mu.1-mu.0)) * t(    pi.1*solve(Sigma)%*%(x.1-mu.1) ) )       
    terms.control <- -solve(Sigma)%*%x.0 - t(as.numeric(t(x.0-mu.0)%*%solve(Sigma)%*%(mu.1-mu.0)) * t(    pi.0*solve(Sigma)%*%(x.0-mu.0) ) )
        ## terms <- cbind((m+n)/m * terms.control, (m+n)/n * terms.case)
    terms <- matrix(nrow=p,ncol=n.0+n.1)
        ## if(dim(    terms[,d==0])[2]==36)browser()
        ## print(dim(    terms[,d==0]))
    terms[,d==0] <- (n.0+n.1)/n.0 * terms.control
    terms[,d==1] <- (n.0+n.1)/n.1 * terms.case
    if(terms.only) return(terms) else return(rowMeans(terms))
}


#' @export
var.beta.gamma.lda <- function(params) {
    p.r <- params$p.r; p.f <- params$p.f
    pi.0 <- params$pi.0; pi.1 <- 1-pi.0
    params.f <- params    
    params.r <- with(params.f, list(mu.0=mu.0[1:p.r],mu.1=mu.1[1:p.r],Sigma.0=Sigma.0[1:p.r,1:p.r],Sigma.1=Sigma.1[1:p.r,1:p.r],pi.0=pi.0))#,n.0=n.0,n.1=n.1))
    vars <- lapply(list(full=params.f,reduced=params.r), function(params)
        with(params, {
            pi.1 <- 1-pi.0
            mu <- mu.1-mu.0
            Sigma <- pi.0*Sigma.0 + pi.1*Sigma.1
            beta.star <- solve(Sigma)%*%mu
            var.control <- solve(Sigma)%*%Sigma.0%*%solve(Sigma) + pi.0^2*(solve(Sigma)%*%Sigma.0%*%beta.star%*%t(beta.star)%*%Sigma.0%*%solve(Sigma)+as.numeric(t(beta.star)%*%Sigma.0%*%beta.star)*solve(Sigma)%*%Sigma.0%*%solve(Sigma))
            var.case <- solve(Sigma)%*%Sigma.1%*%solve(Sigma) + pi.1^2*(solve(Sigma)%*%Sigma.1%*%beta.star%*%t(beta.star)%*%Sigma.1%*%solve(Sigma)+as.numeric(t(beta.star)%*%Sigma.1%*%beta.star)*solve(Sigma)%*%Sigma.1%*%solve(Sigma))
            ## var.control/n.0 + var.case/n.1
            var.control/pi.0 + var.case/pi.1
        }))
    vars$reduced <- cbind(rbind(vars$reduced,matrix(0,p.f-p.r,p.r)),matrix(0,p.f,p.f-p.r))
    ## vars[['reduced']]
    Sigma.f <- with(params.f, pi.0*Sigma.0 + (1-pi.0)*Sigma.1)
    Sigma.r <- Sigma.f[1:p.r,1:p.r]
    Sigma.0.f <- params.f$Sigma.0;     Sigma.1.f <- params.f$Sigma.1
    beta.star.f <- with(params.f, solve(Sigma.f)%*%(mu.1-mu.0))
    ## n.0 <- params.f$n.0; n.1 <- params.f$n.1
    covs <- lapply(list(control=list(Sigma.d.f=Sigma.0.f,pi=pi.0),case=list(Sigma.d.f=Sigma.1.f,pi=pi.1)), function(params) with(params, {
        A.d <- pi^2*Sigma.d.f%*%beta.star.f%*%t(beta.star.f)
        solve(Sigma.f)%*%( (1+sum(diag(A.d)))*diag(p.f) + A.d)%*%Sigma.d.f[,1:p.r]%*%solve(Sigma.r)        
    }))
    ## cov <- covs$control/n.0+covs$case/n.1
    cov <- covs$control/pi.0+covs$case/pi.1
    ## rbind(cbind(vars$reduced, t(cov)), cbind(cov,vars$full))
    cov <- rbind(t(cov),matrix(0,p.f-p.r,p.f))
    rbind(cbind(vars$reduced, cov), cbind(t(cov),vars$full))
}


#' @export
cov.beta.gamma.hajek.normal.lda  <- function(params) {
    f <- function(beta,mu,Sigma){
        q <- function(beta)c(t(beta)%*%Sigma%*%beta)
        pnorm( sum(mu*beta)/sqrt(1+q(beta)) )
    }
    f.grad <- function(beta,mu,Sigma){
        q <- function(beta)c(t(beta)%*%Sigma%*%beta)
        ## pnorm( sum(mu*beta)/sqrt(1+q(beta)) )
        c( dnorm(sum(mu*beta)/sqrt(1+q(beta))) * ( mu/sqrt(1+q(beta))  -sum(mu*beta)/(1+q(beta))^(3/2) * Sigma%*%beta ) )
    }
    f.hessian <- function(beta,mu,Sigma){
        ## browser()
        q <- function(beta)c(t(beta)%*%Sigma%*%beta)
        aa <- mu/sqrt(1+q(beta))  -sum(mu*beta)/(1+q(beta))^(3/2) * Sigma%*%beta 
        dnorm(sum(mu*beta)/sqrt(1+q(beta))) * sum(mu*beta)/sqrt(1+q(beta)) * aa%*%t(aa) - dnorm(sum(mu*beta)/sqrt(1+q(beta)))*(1+q(beta))^(-3/2)*( -mu%*%t(beta)%*%Sigma - Sigma%*%beta%*%t(mu) - sum(beta*mu)*Sigma + 3/(1+q(beta))*sum(beta*mu)*Sigma%*%beta%*%t(beta)%*%Sigma)
    }
    ## browser()
    p <- p.f <- params$p.f; p.r <- params$p.r
    pi.0 <- params$pi.0; pi.1 <- 1-pi.0
    mu.0 <- params$mu.0; mu.1 <- params$mu.1
    mu <- mu.1-mu.0
    Sigma.0 <- params$Sigma.0;         Sigma.1 <- params$Sigma.1
    Sigma <- pi.0*Sigma.0 + pi.1*Sigma.1
    beta.star <- solve(Sigma)%*%mu
    q <- c(t(beta.star)%*%Sigma%*%beta.star)
    q.0 <- c(t(beta.star)%*%Sigma.0%*%beta.star)
    q.1 <- c(t(beta.star)%*%Sigma.1%*%beta.star)
    by.class <- lapply(list(control=list(q.d = q.1, Sigma.d = Sigma.0, Sigma.not.d = Sigma.1, pi.d = pi.0), case=list(q.d = q.0, Sigma.d = Sigma.1, Sigma.not.d = Sigma.0, pi.d = pi.1)), function(params) 
        with(params, {
            linear.term <- J.1(beta.star,mu/sqrt(q.d),Sigma.d/q.d) - outer(mu/sqrt(q.d),f.grad(beta.star,mu/sqrt(q.d),Sigma.d/q.d))
            quad.term.inner <- f.hessian(beta.star,mu/sqrt(q.d),Sigma.d/q.d) - (Sigma.d/q.d+mu%*%t(mu)/q.d)%*%beta.star%*%t(f.grad(beta.star,mu/sqrt(q.d),Sigma.d/q.d))
            quad.term <- q.d*pi.d*solve(Sigma)%*%quad.term.inner%*%(beta.star%*%t(beta.star)%*%Sigma.not.d/q.d-diag(p))  
            - (  (sqrt(q.d)*solve(Sigma) + sqrt(q.d)*sum(mu*beta.star)*pi.d*solve(Sigma) + sqrt(q.d)*pi.d*solve(Sigma)%*%mu%*%t(beta.star)) %*% linear.term %*% (beta.star%*%t(beta.star)%*%Sigma.not.d/q.d-diag(p)) - quad.term  ) 
        }))
    cov.f <- with(by.class, control/pi.0 + case/pi.1)
    cov.r <- solve(Sigma[1:p.r,1:p.r])%*%diag(p.f)[1:p.r,]%*%Sigma%*%cov.f
    cov.r <- rbind(cov.r,matrix(0,p.f-p.r,p.f))
    ## return(list(reduced=cov.r,full=cov.f))
    rbind(cov.r,cov.f)
}

## cov.beta.gamma.hajek.normal.lda <- function(beta,params) {
##     p <- length(beta)
##     Sigma.0 <- params$Sigma.0; Sigma.1 <- params$Sigma.1
##     if(is.null(Sigma.0)||is.null(Sigma.1)) Sigma.0 <- Sigma.1 <- params$Sigma
##     Sigma <- (Sigma.0+Sigma.1)/2
##     cov.beta.hajek <- matrix(0,p,p)
##     ## cov.beta.hajek[1:p.red,1:p.red] <- cov.coefs.hajek.grad.normal(beta,params)[1:p.red,1:p.red]#vcov.full.hajek.lda(beta,params)[1:p.red,1:p.red] * (1/m+1/n)
##     cov.beta.hajek[1:p.red,1:p.red] <- solve(Sigma[1:p.red,1:p.red])%*%Sigma[1:p.red,1:p.red]%*%(cov.gamma.hajek.normal(beta,params)[1:p.red,1:p.red])
## }

## this is for the IF that assumes Sigma.0,Sigma.1 are available
## var.beta.gamma.lda <- function(beta,p.reduced,params) { # m,n should be params
##     p <- length(beta)
##     m <- params$m; n <- params$n
##     Sigma.0 <- params$Sigma.0; Sigma.1 <- params$Sigma.1
##     if(is.null(Sigma.0)||is.null(Sigma.1)) Sigma.0 <- Sigma.1 <- params$Sigma
##     Sigma <- (Sigma.0+Sigma.1)/2
##     var.gamma <- solve(Sigma)%*%(Sigma.0/m+Sigma.1/n)%*%solve(Sigma)
##     var.beta <- matrix(0,p,p)
##     var.beta[1:p.reduced,1:p.reduced] <- solve(Sigma[1:p.reduced,1:p.reduced])%*%((Sigma.0/m+Sigma.1/n)[1:p.reduced,1:p.reduced])%*%solve(Sigma[1:p.reduced,1:p.reduced])
##     cov.beta.gamma <- solve(Sigma[1:p.reduced,1:p.reduced]) %*% ((Sigma.0/m+Sigma.1/n)[1:p.reduced,]  ) %*% solve(Sigma)
##     cov.beta.gamma <- rbind(cov.beta.gamma,matrix(0,p-p.reduced,p))
##     vcov.beta.gamma <- rbind(cbind(var.beta,cov.beta.gamma),cbind(t(cov.beta.gamma),var.gamma))
## }


## version assuming known class Sigmas (think this is a wrong description--safe to delte this code?)x
## var.beta.gamma.lda <- function(p.reduced,params) {
##     p.r <- p.reduced; p.f <- length(params$mu.1)
##     params.f <- params    
##     params.r <- with(params.f, list(mu.0=mu.0.f[1:p.r],mu.1=mu.1[1:p.r],Sigma.0=Sigma.0[1:p.r,1:p.r],Sigma.1=Sigma.1[1:p.r,1:p.r],pi.0=pi.0,n.0=n.0,n.1=n.1))
##     vars <- lapply(list(full=params.f,reduced=params.r), function(params)
##         with(params, {
##             pi.1 <- 1-pi.0
##             mu <- mu.1-mu.0
##             Sigma <- pi.0*Sigma.0 + pi.1*Sigma.1
##             beta.star <- solve(Sigma)%*%mu
##             var.control <- solve(Sigma)%*%Sigma.0%*%solve(Sigma) + pi.0^2*(solve(Sigma)%*%Sigma.0%*%beta.star%*%t(beta.star)%*%Sigma.0%*%solve(Sigma)+as.numeric(t(beta.star)%*%Sigma.0%*%beta.star)*solve(Sigma)%*%Sigma.0%*%solve(Sigma))
##             var.case <- solve(Sigma)%*%Sigma.1%*%solve(Sigma) + pi.1^2*(solve(Sigma)%*%Sigma.1%*%beta.star%*%t(beta.star)%*%Sigma.1%*%solve(Sigma)+as.numeric(t(beta.star)%*%Sigma.1%*%beta.star)*solve(Sigma)%*%Sigma.1%*%solve(Sigma))
##             var.control/n.0 + var.case/n.1
##         }))
##     ## vars[['reduced']]
##     Sigma.f <- with(params.f, pi.0*Sigma.0.f + (1-pi.0)*Sigma.1.f)
##     Sigma.r <- Sigma.f[1:p.r,1:p.r]
##     beta.star.f <- with(params.f, solve(Sigma.f)%*%(mu.1.f-mu.0.f))
##     n.0 <- params.f$n.0; n.1 <- params.f$n.1
##     covs <- lapply(list(control=list(Sigma.d.f=Sigma.0.f,pi=pi.0),case=list(Sigma.d.f=Sigma.1.f,pi=pi.1)), function(params) with(params, {
##         A.d <- pi^2*Sigma.d.f%*%beta.star.f%*%t(beta.star.f)
##         solve(Sigma.f)%*%( (1+sum(diag(A.d)))*diag(p.f) + A.d)%*%Sigma.d.f[,1:p.r]%*%solve(Sigma.r)        
##     }))
##     cov <- covs$control/n.0+covs$case/n.1
##     rbind(cbind(vars$reduced, t(cov)), cbind(cov,vars$full))
## }


#' Variance-covariance matrix for LDA influence functions and Hoeffding gradient
#'
#' @param beta.star the probability limit of the LDA estimates
#' @param params parameters of the normal model, which are estimated
#'     if NULL value is passed
#' @return the variance-covariance matrix of the reduced coefficients, full coefficients, and Hoeffding gradient
#' @export
var.normal.lda <- function(beta.star,params) {
    vcov.beta.gamma <- var.beta.gamma.lda(params)# / (n.0+n.1)
    ## var.beta <- vcov.beta.gamma[1:p.f,1:p.f]
    ## var.gamma <- vcov.beta.gamma[(p.f+1):(2*p.f),(p.f+1):(2*p.f)]
    cov.beta.gamma.hajek <- cov.beta.gamma.hajek.normal.lda(params) #/ (n.0+n.1)
    var.hajek <- var.hajek.normal(beta.star,params) #/ (n.0+n.1)
    ## var.try <- t(a)%*%var.beta%*%a
    ## var.try <- t(a)%*%var.gamma%*%a
    ## var.try <- t(b)%*%vcov.beta.gamma%*%b
    ## var.try <- t(a.f)%*%var.hajek%*%a.f
    rbind(cbind(vcov.beta.gamma,cov.beta.gamma.hajek),cbind(t(cov.beta.gamma.hajek),var.hajek))
}

## ## deprecate in favor of qaucdiff
## qdelta.normal.lda <- function(p,beta.star.f,p.r,params.f,reps=1e3) {
##     vcov.beta.gamma.hajek  <-  var.normal.lda(beta.star.f,p.r,params.f) 
##     auc.hessian <- auc.hessian.normal(beta.star.f,params.f)
##     beta.gamma.hajek <- rmvnorm(reps,sigma=vcov.beta.gamma.hajek)
##     beta <- beta.gamma.hajek[,1:p.f]
##     gamma <- beta.gamma.hajek[,(p.f+1):(2*p.f)]
##     hajek <- beta.gamma.hajek[,(2*p.f+1):(3*p.f)]
##     reference.hajek <- rowSums((beta-gamma)*hajek)
##     reference.taylor <-   (rowSums((beta%*%auc.hessian)*beta)/2 - rowSums((gamma%*%auc.hessian)*gamma)/2)
##     reference <- reference.hajek + reference.taylor
##     quantile(reference,p)
## }






## ## specific to mrc betahat. in particular assuming homoskedastic data.

## ## hessian of mrc auc at c.beta.LDA, a simpler formula than
## ## auc.hessian.normal. the "bread" of the sandwich variance of betahat.
## bread.mrc <- function(params) { 
##     mu <- with(params,mu.1-mu.0)
##     Sigma <- params$Sigma
##     q <- function(beta)as.numeric(t(beta)%*%Sigma%*%beta)
##     beta.star <- params$beta.star
##     beta.lda <- solve(Sigma)%*%mu
##     c <- 2/beta.lda[1]
##     V <- dnorm(sqrt(2*q(beta.star))/c)*(Sigma%*%beta.star%*%t(beta.star)%*%Sigma/q(beta.star) - Sigma)*sqrt(2/q(beta.star))/c
##     solve(V[-1,-1])
## }

## ## return value with effective # of params ie p-1
## infl.normal.mrc <- function(x,y,params) {
##     ## browser()
##     mu <- with(params, mu.1-mu.0)
##     Sigma <- params$Sigma
##     q <- function(beta)as.numeric(t(beta)%*%Sigma%*%beta)
##     beta.star <- params$beta.star
##     ## beta.lda <- solve(Sigma)%*%mu
##     ## c <- 1/beta.lda[1]
##     ## V <- dnorm(sqrt(2*q(beta.star))/c)*(Sigma%*%beta.star%*%t(beta.star)%*%Sigma/q(beta.star) - Sigma)*sqrt(2/q(beta.star))/c
##     ## -solve(V[-1,-1]) %*% (hajek.grad.normal(x,y,beta.star,params)[-1])
##     ## bread <- bread.mrc(params)
##     params <- c(params,list(Sigma.0=Sigma,Sigma.1=Sigma))
##     bread <- solve(auc.hessian.normal(beta.star,params)[-1,-1])
##     -bread %*% (hajek.grad.normal(x,y,beta.star,params)[-1])
## }


## vcov.beta.gamma.mrc <- function(beta,params) {
##     ## browser()
##     ## bread.mrc <- function(params) { # hessian of mrc auc
##     ##     mu <- with(params,mu.1-mu.0)
##     ##     Sigma <- params$Sigma
##     ##     q <- function(beta)as.numeric(t(beta)%*%Sigma%*%beta)
##     ##     beta.star <- params$beta.star
##     ##     beta.lda <- solve(Sigma)%*%mu
##     ##     c <- 1/beta.lda[1]
##     ##     V <- dnorm(sqrt(2*q(beta.star))/c)*(Sigma%*%beta.star%*%t(beta.star)%*%Sigma/q(beta.star) - Sigma)*sqrt(2/q(beta.star))/c
##     ##     solve(V[-1,-1])
##     ## }
##     p.r <- params$p.r
##     m <- params$m; n <- params$n
##     Sigma.f <- params$Sigma
##     mu.f  <- params$mu <- with(params,mu.1-mu.0); 
##     beta.f <- params$beta.star; beta.r <- beta.f[1:p.r]
##     params.f <- params
##     params.r <- with(params.f, list(mu.0=mu.0[1:p.r],mu.1=mu.1[1:p.r],Sigma=Sigma[1:p.r,1:p.r],beta.star=beta.star[1:p.r],m=m,n=n))
##     q <- as.numeric(t(beta.f)%*%Sigma.f%*%beta.f)
##     bread.r <- bread.mrc(params.r); bread.f <- bread.mrc(params.f)
##     cov.beta.gamma <- matrix(0,p.f-1,p.f-1)
##     cov.beta.gamma[seq(p.r-1),seq(p.f-1)] <- bread.r%*%J.4(beta.star.f,mu.f,Sigma.f,p.r)[-1,-1]%*%bread.f*(1/m+1/n)
##     meat.f <- var.hajek.normal(beta.star.f,params.f)[-1,-1]
##     var.gamma <- bread.f %*% meat.f %*% bread.f
##     meat.r <- var.hajek.normal(beta.star.r,params.r)[-1,-1]
##     var.beta <- matrix(0,p.f-1,p.f-1)
##     var.beta[seq(p.r-1),seq(p.r-1)] <- bread.r %*% meat.r %*% bread.r
##     rbind(cbind(var.beta,cov.beta.gamma),cbind(t(cov.beta.gamma),var.gamma))
## }



## cov.beta.gamma.hajek.normal.mrc <- function(beta,params) {
##     ## browser()
##     p.r <- params$p.r; p.f <- params$p.f
##     vcov.beta.gamma <- vcov.beta.gamma.mrc(beta,params)
##     bread.f <- bread.mrc(params.f)
##     cov.beta.gamma <- vcov.beta.gamma[seq(p.f-1),p.f:(2*p.f-2)]
##     cov.beta.hajek <- -cov.beta.gamma%*%solve(bread.f)
##     var.gamma <- vcov.beta.gamma[p.f:(2*p.f-2),p.f:(2*p.f-2)]
##     cov.gamma.hajek <- -var.gamma%*%solve(bread.f)
##     rbind(cov.beta.hajek,cov.gamma.hajek)
##     }




## vcov.beta.gamma.mrc <- function(beta,params) {
##     p.r <- params$p.r
##     m <- params$m; n <- params$n
##     Sigma.f <- params$Sigma
##     mu.f  <- params$mu <- with(params,mu.1-mu.0); 
##     beta.f <- params$beta.star; beta.r <- beta.f[1:p.r]
##     params.f <- c(params,list(Sigma.0=Sigma.f,Sigma.1=Sigma.f))
##     params.r <- with(params.f, list(mu.0=mu.0[1:p.r],mu.1=mu.1[1:p.r],Sigma.0=Sigma.f[1:p.r,1:p.r],Sigma.1=Sigma.f[1:p.r,1:p.r],beta.star=beta.star[1:p.r],m=m,n=n))
##     q <- as.numeric(t(beta.f)%*%Sigma.f%*%beta.f)
##     ## bread.r <- bread.mrc(params.r); bread.f <- bread.mrc(params.f)
##     bread.f <- solve(auc.hessian.normal(beta.star.f,params.f)[-1,-1])
##     bread.r <- solve(auc.hessian.normal(beta.star.r,params.r)[-1,-1])
##     cov.beta.gamma <- matrix(0,p.f-1,p.f-1)
##     ## cov.beta.gamma[seq(p.r-1),seq(p.f-1)] <- bread.r%*%J.4(beta.star.f,mu.f,Sigma.f,p.r)[-1,-1]%*%bread.f*(1/m+1/n)
##     cov.beta.gamma[seq(p.r-1),seq(p.f-1)] <- bread.r%*%var.hajek.normal(beta.star.f,params.f)[-1,-1][seq(p.r-1),]%*%bread.f
##     meat.f <- var.hajek.normal(beta.star.f,params.f)[-1,-1]
##     var.gamma <- bread.f %*% meat.f %*% bread.f
##     meat.r <- var.hajek.normal(beta.star.r,params.r)[-1,-1]
##     var.beta <- matrix(0,p.f-1,p.f-1)
##     var.beta[seq(p.r-1),seq(p.r-1)] <- bread.r %*% meat.r %*% bread.r
##     rbind(cbind(var.beta,cov.beta.gamma),cbind(t(cov.beta.gamma),var.gamma))
## }



## cov.beta.gamma.hajek.normal.mrc <- function(beta,params) {
##     ## browser()
##     p.r <- params$p.r; p.f <- params$p.f
##     vcov.beta.gamma <- vcov.beta.gamma.mrc(beta,params)
##     params.f <- c(params,list(Sigma.0=Sigma.f,Sigma.1=Sigma.f))
##     bread.f <- solve(auc.hessian.normal(beta.star.f,params.f)[-1,-1])
##     cov.beta.gamma <- vcov.beta.gamma[seq(p.f-1),p.f:(2*p.f-2)]
##     cov.beta.hajek <- -cov.beta.gamma%*%solve(bread.f)
##     var.gamma <- vcov.beta.gamma[p.f:(2*p.f-2),p.f:(2*p.f-2)]
##     cov.gamma.hajek <- -var.gamma%*%solve(bread.f)
##     rbind(cov.beta.hajek,cov.gamma.hajek)
##     }

## TODO  deprecate dgm.params argument, instead should be passed to vcov.fn and hessian.fn 

#' Fit a model for the difference of AUCs
#'
#' @param w covariates
#' @param d class labels, 0 or 1
#' @param p.reduced the number of initial components of w used in the reduced model
#' @param N the total sample size, taken to be nrow(w) if w isn't NULL
#' @param infl.fn the influence function to use; see details for specification
#' @param vcov.fn a function returning the variance-covariance matrix of the coefficient estimates and hoeffding gradient; see details for specification
#' @param hessian.fn a function returning the hessian of the AUC; see details for specification
#' @param eps.finite.diff the step size for numeric differentiation of the Hessian (if hessian.fn is NULL) and Hoeffding gradient
#' @param dgm.params parameters of the data-generating mechanism to be passed to vov.fn and hessian.fn
#' @return  an object of class "aucdiff"
#' @export
aucdiff <- function(w=NULL,d=NULL,beta,p.reduced,N=NULL,infl.fn,vcov.fn=NULL,hessian.fn=NULL,eps.finite.diff=.05,dgm.params) {
    ## browser()
    p.r <- p.reduced
    p.f <- length(beta)
    if(is.null(N)) N <- nrow(w) ## switch to n.0/n.1/n notation
    ## if(is.null(infl.fn)||is.null(vcov.fn)) {
    ## if(is.null(vcov.fn)) {
    x.0 <- w[d==0,]; x.1 <- w[d==1,]
    m <- sum(d); n <- sum(1-d)
    ## }
    ## vcov.beta.gamma.hajek  <-  var.normal.lda(beta.star.f,p.r,dgm.params) 
    ## auc.hessian <- auc.hessian.normal(beta.star.f,dgm.params)
    if(!is.null(hessian.fn)) {
        auc.hessian <- hessian.fn(beta,dgm.params)
    } else {
        auc.hessian <- nlme::fdHess(beta,function(beta)auc.hat(x.0%*%beta,x.1%*%beta),.relStep=eps.finite.diff,minAbsPar=1)$Hessian
    }
    if(!is.null(vcov.fn)) {
        vcov.beta.gamma.hajek  <-  vcov.fn(beta,dgm.params) 
    } else {
        infl.f <- infl.fn(w,d)#,dgm.params)
        infl.r <- infl.fn(w[,1:p.r],d)#,dgm.params)
        infl.r <- rbind(infl.r,matrix(0,nrow=p.f-p.r,ncol=ncol(infl.r)))
        hajek.grad.case.hat <- apply(x.1, 1, function(y) nlme::fdHess(beta,function(beta)mean(x.0%*%beta <  as.numeric(y%*%beta)),.relStep=eps.finite.diff, minAbsPar=1)$gradient )
        hajek.grad.control.hat <- apply(x.0, 1, function(x) nlme::fdHess(beta,function(beta)mean(as.numeric(x%*%beta) <  x.1%*%beta),.relStep=eps.finite.diff, minAbsPar=1)$gradient )
        ## hajek.grad.hat <- cbind( hajek.grad.control.hat*(m+n)/m,hajek.grad.case.hat*(m+n)/n)
        hajek.grad.hat <- matrix(nrow=p.f,ncol=m+n)
        hajek.grad.hat[,d==0] <- hajek.grad.control.hat*(m+n)/m
        hajek.grad.hat[,d==1] <- hajek.grad.case.hat*(m+n)/n
        var.hajek.hat <- var(t(hajek.grad.hat))#/(m+n)
        cov.gamma.hajek.hat <- cov(t(infl.f),t(hajek.grad.hat))#/(m+n)
        cov.beta.hajek.hat <- cov(t(infl.r),t(hajek.grad.hat))#/(m+n)
        cov.beta.gamma.hajek.hat <- rbind(cov.beta.hajek.hat,cov.gamma.hajek.hat)
        ## cov.beta.hajek.hat <- rbind(cov.beta.hajek.hat,matrix(0,nrow=p.f-p.r,ncol=p.f))
        vcov.beta.gamma.hat <- cov(cbind(t(infl.r),t(infl.f))) #/ (m+n)
        vcov.beta.gamma.hajek <- rbind(cbind(vcov.beta.gamma.hat,cov.beta.gamma.hajek.hat),cbind(t(cov.beta.gamma.hajek.hat),var.hajek.hat))        
    }
    structure(list(vcov.beta.gamma.hajek=vcov.beta.gamma.hajek, auc.hessian=auc.hessian,p.f=p.f,p.r=p.r,N=N), class='aucdiff')
}


#' Generate a random sample of AUC differences
#'
#' @param n the number of samples
#' @param fit an object of class "aucdiff"
#' @return  numeric vector consisting of samples
#' @export
raucdiff <- function(n, fit) {
    p.f <- fit$p.f; p.r <- fit$p.r
    N <- fit$N
    beta.gamma.hajek <- mvtnorm::rmvnorm(n,sigma=fit$vcov.beta.gamma.hajek/N)
    beta <- beta.gamma.hajek[,1:p.f]
    gamma <- beta.gamma.hajek[,(p.f+1):(2*p.f)]
    hajek <- beta.gamma.hajek[,(2*p.f+1):(3*p.f)]
    reference.hajek <- rowSums((beta-gamma)*hajek)
    reference.taylor <-   (rowSums((beta%*%fit$auc.hessian)*beta)/2 - rowSums((gamma%*%fit$auc.hessian)*gamma)/2)
    ## plot(ecdf(reference.hajek+reference.taylor),add=TRUE,lty=2,col=2)
    reference <- reference.hajek+reference.taylor
    ## ecdf(reference.hajek+reference.taylor)(q)
}

#' Approximate CDF of AUC differences
#'
#' @param q  quantile
#' @param fit an object of class "aucdiff"
#' @param mc.reps the number of samples used to estimate the CDF
#' @return  the approximate probability of an AUC difference being less than or equal to q
#' @export
paucdiff <- function(q, fit, mc.reps=5e3) {
    p.f <- fit$p.f; p.r <- fit$p.r
    N <- fit$N
    beta.gamma.hajek <- mvtnorm::rmvnorm(mc.reps,sigma=fit$vcov.beta.gamma.hajek/N)
    beta <- beta.gamma.hajek[,1:p.f]
    gamma <- beta.gamma.hajek[,(p.f+1):(2*p.f)]
    hajek <- beta.gamma.hajek[,(2*p.f+1):(3*p.f)]
    reference.hajek <- rowSums((beta-gamma)*hajek)
    reference.taylor <-   (rowSums((beta%*%fit$auc.hessian)*beta)/2 - rowSums((gamma%*%fit$auc.hessian)*gamma)/2)
    ## plot(ecdf(reference.hajek+reference.taylor),add=TRUE,lty=2,col=2)
    reference <- reference.hajek+reference.taylor
    ecdf(reference)(q)
    ## ecdf(reference.hajek+reference.taylor)(q)
}



#' Approximate quantile function of AUC differences
#'
#' @param p probability
#' @param fit an object of class "aucdiff"
#' @param mc.reps the number of samples used to estimate the quantile
#' @return  an approximation to the number q such that the probability of an AUC difference being less than or equal to q is p
#' @export
qaucdiff <- function(p, fit, mc.reps=5e3) {
    p.f <- fit$p.f; p.r <- fit$p.r
    N <- fit$N
    beta.gamma.hajek <- mvtnorm::rmvnorm(mc.reps,sigma=fit$vcov.beta.gamma.hajek/N)
    beta <- beta.gamma.hajek[,1:p.f]
    gamma <- beta.gamma.hajek[,(p.f+1):(2*p.f)]
    hajek <- beta.gamma.hajek[,(2*p.f+1):(3*p.f)]
    reference.hajek <- rowSums((beta-gamma)*hajek)
    reference.taylor <-   (rowSums((beta%*%fit$auc.hessian)*beta)/2 - rowSums((gamma%*%fit$auc.hessian)*gamma)/2)
    ## plot(ecdf(reference.hajek+reference.taylor),add=TRUE,lty=2,col=2)
    reference <- reference.hajek+reference.taylor
    quantile(reference,p)
}

## ecdf.aucdiff <- function(fit,mc.reps=5e3) {
##     p.f <- fit$p.f; p.r <- fit$p.r
##     beta.gamma.hajek <- mvtnorm::rmvnorm(mc.reps,sigma=fit$vcov.beta.gamma.hajek)
##     beta <- beta.gamma.hajek[,1:p.f]
##     gamma <- beta.gamma.hajek[,(p.f+1):(2*p.f)]
##     hajek <- beta.gamma.hajek[,(2*p.f+1):(3*p.f)]
##     reference.hajek <- rowSums((beta-gamma)*hajek)
##     reference.taylor <-   (rowSums((beta%*%fit$auc.hessian)*beta)/2 - rowSums((gamma%*%fit$auc.hessian)*gamma)/2)
##     ## plot(ecdf(reference.hajek+reference.taylor),add=TRUE,lty=2,col=2)
##     ecdf(reference.hajek+reference.taylor)
##     }


#' Plot the approximate CDF of the difference of AUCs
#'
#' @param fit an object of class "aucdiff"
#' @param mc.reps the number of samples used to estimate the CDF
#' @param ... additional parameters passed down to the graphics
#'     routine
#' @return side effect
#' @export
plot.aucdiff <- function(fit, mc.reps=5e3,...) {
    ecdf.aucdiff <- function(x)paucdiff(x,fit,mc.reps)
    curve(ecdf.aucdiff,...)
    }


#' An approximate CI for the difference of AUCs
#'
#' @param fit an object of class "aucdiff"
#' @param level level of the CI
#' @param mc.reps the number of samples used to estimate the quantile
#' @return a pair of numbers giving the lower and upper endpoints of
#'     the CI
#' @export
confint.aucdiff <- function(fit,level=.95,mc.reps=5e3) {
    alpha <- 1-level
    ## qaucdiff(c(alpha/2,1-alpha/2),fit,mc.reps)
    c( qaucdiff(alpha,fit,mc.reps), Inf)
    ## c( -Inf, qaucdiff(1-alpha,fit,mc.reps))
    }

#' Influence function for GLM estimates
#'
#' @param w covariates
#' @param d class labels, 0 or 1
#' @param params parameters of the GLM, including the link function
#'     and its first and second derivatives
#' @param terms.only boolean, whether to return the terms of the
#'     influence function directly or a class-weighted average
#' @return depending on "terms.only", either the terms of the
#'     influence function directly or a class-weighted average
#' @export
infl.glm <- function(x,d,params,terms.only=TRUE) {
    ## p <- params$p.f
    ## browser()
    if(is.null(params$beta))
        params$beta <- coef(glm(d~x-1,family=binomial(link=params$link.name)))
    score <- function(x,d,params) { # x in model matrix format
        h <- params$link; h.1 <- params$link.deriv
        eta <- as.numeric(x%*%params$beta)
        t( (d-h(eta))/h(eta)/(1-h(eta))*h.1(eta) * x )
    }
    fi <- function(x,d,params) { # x in model matrix format
        h <- params$link; h.1 <- params$link.deriv; h.2 <- params$link.deriv2
        n <- nrow(x)
        eta <- as.numeric(x%*%params$beta)
        denom <- h(eta)*(1-h(eta))
        -t((h.2(eta)*(d-h(eta))/denom - h.1(eta)/denom^2 * (d*(h.1(eta)-2*h(eta)*h.1(eta))+h.1(eta)*h(eta)^2))*x) %*% x / n
    }
    ## terms <- solve(fi(x,g,params)[1:p,1:p])%*%score(x,g,params)[1:p,]
    terms <- solve(fi(x,d,params))%*%score(x,d,params)
    if (terms.only) return(terms) else return(rowMeans(terms))
}
