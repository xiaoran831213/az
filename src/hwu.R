library(CompQuadForm)

hwu.core<-function(y,K,geno,X=NULL,center.geno=F,gsim=c("add","eq","dist"),appx=c("davis","normal"))
{
    if(center.geno) geno<-apply(geno,2,map.std)
    n<-length(as.vector(y));
    y<-rank(y); y<-(y-mean(y))/sd(y);
    
    if(is.null(X))
    {
        X<-cbind(rep(1,n));
    }
    else
    {
        X<-cbind(1,X);
    }
    inX<-solve(t(X)%*%X);
    P<-diag(1,n)-X%*%inX%*%t(X);
    e<-as.vector(P%*%y);
    y<-e/sqrt(sum(e^2)/(length(e)-ncol(X)));
    
    #y<-rank(y); y<-(y-mean(y))/sd(y);
    
    mt.trt<-y%*%t(y);
    
    if(gsim=="add")
    {
        mt.geno=geno%*%t(geno)
    }
    else if(gsim=="eq")
    {
        mt.geno=outer(geno[,1],geno[,1],"==")*1
    }
    else if(gsim=="dist")
    {
        mt.geno=hwu.weight.gaussian(geno);
    }
    
    wt<- ( mt.geno ) * K
    #wt<- K;
    
    diag(wt)<-0;
    #diag(mt.trt)<-0;
    
    U<-sum(mt.trt*wt);
    
    wt<-wt - X%*%inX%*% (t(X)%*%wt);
    wt<-wt - (wt%*%X) %*%inX%*%t(X);
    
    appx<-match.arg(appx);
    
    if(appx=="davis")
    {
        sg<-eigen(wt,symmetric=TRUE, only.values = TRUE)
        coef<-sg$values;
        
        p.value<-davies(U,coef,acc=0.000001)$Qq;
        
    }else if(appx=="normal")
    {
        miu<-sum(diag(wt));
        var<-2*sum(wt^2);
        Z<-(U-miu)/sqrt(var)
        
        p.value<-pnorm(Z,lower.tail=F)
    }
    return (list(U=U,p=p.value,coef=coef))
}

hwu.run<-function(Tn,K,geno,X=NULL,center.geno=F,gsim=c("add","eq","dist"),appx=c("davis","normal"))
{
    appx<-match.arg(appx);
    gsim<-match.arg(gsim);
    hwu.core(Tn,K,geno,X,center.geno,gsim=gsim,appx=appx);
}

## x ---- covariate
## y ---- response variable
## f ---- U kernel
## r ---- residual matrix
## w, ... ---- weight terms.
hwu.dg2 <- function(y, w, x=NULL)
{
    ## response and covariate
    M <- length(y);

    ## standardize y to m=0, s=1
    y <- rank(y);
    y <- (y-mean(y))/sd(y);
    
    ## regression residual matrix, R = I - X(X'X)^X'
    if(is.null(x))
    {
        x = matrix(1, M, 1L)
    }
    else
    {
        x = cbind(1, x)
    }
    hat <- diag(1, M, M) - tcrossprod(x %*% solve(crossprod(x)), x);
    
    ## exclude liner covariant effect on y, leave residual of Y
    y <- hat %*% y;
    y <- y/sqrt(sum(y^2)/(M-ncol(x)));
    
    ## the U kernel is the pair wise similarity between phenotypes
    f <- tcrossprod(y);
    
    diag(w) <- 0; # ??
    
    ## compute U score
    u <- sum(w * f);

    ## exclude coveriant and intercept effect on both
    ## dimensions of w
    w <- tcrossprod(hat %*% w, hat);

    ## calculate p-value of u.
    pval <- tryCatch(
    {
        coef <- eigen(w, symmetric=T, only.values=T)$values;
        p = davies(u, coef, acc=1e-6)$Qq
        p
    }, warning = function(wa)
    {
        print(wa)
        p
    }, error = function(e)
    {
        print(e)
        NA
    })
    pval
}

## weight centeralizer
.wct <- function(w)
{
    w - outer(rowMeans(w), colMeans(w), '+') + mean(w)
}

.map.std.norm <- function(y) qnorm((rank(y)-0.5)/length(y))

.map.std <- function(y)
{
    y<-y-mean(y)
    sdy<-sd(y)
    if(sdy!=0) y<-y/sdy;
    y
}

## the gussian distance kernel
.hwu.GUS <- function(x, w = rep(1, ncol(x)), std = FALSE)
{
    ## normalize features
    x <- apply(x, 2L, .map.std.norm);
    
    ## exp(- weighted gaussian distance) = weight gaussian similiarity
    ## centralize the similarity.
    m <- exp(-dist(scale(x, F, sqrt(2 * sum(w) / w)), method='euclidean')^2)
    m <- as.matrix(m)
    diag(m) <- 1
    if(std)
    {
        a <- min(m)
        b <- max(m)
        m <- (m - a) / (b - a)
    }
    m
}

.hwu.IBS <- function(x, w = rep(1, ncol(x)), lv = 2, std = FALSE)
{
    ## centred weight IBS, scale also coerce dist to matrix
    ## x <- apply(x, 2L, .map.std.norm)
    m <- 1 - dist(scale(x, F, lv * sum(w) / w), method='manhattan')
    m <- as.matrix(m)
    diag(m) <- 1L
    if(std)
    {
        a <- min(m)
        b <- max(m)
        m <- (m - a) / (b - a)
    }
    m
}

## x: variant row major genetic matrix
## p: alternative allele frequency
.hwu.GRM <- function(x, p=NULL, lv = 2, std = c('sigmoid', 'scale01', 'no'))
{
    if(!is.matrix(x))
        stop('x is not a matrix')
    if(is.null(p))
        p <- colMeans(x, na.rm=T)/lv

    x <- scale(x, 2 * p, sqrt(2 * p * (1 - p)))
    m <- tcrossprod(x) / ncol(x)
    std <- match.arg(std)
    if(std == 'sigmoid')
        1 / (1 + exp(-m))
    else if(std == 'scale01')
    {
        a <- min(m)
        b <- max(m)
        (m - a) / (b - a)
    }
    else
        m
}

hwu.weight.cov<-function(x, w = NULL)
{
    if(!is.matrix(x))
        stop('x is not a matrix')

    ## mormalize feature weights
    if(is.null(w))
        w <- rep(1L, ncol(x))
    x<-apply(x, 2L, map.std.norm);
    if(is.null(w)) w<-rep(1,ncol(x));
    w<-w/sum(w);
    
    m<-matrix(0,nrow=nrow(x),ncol=nrow(x));
    for(i in 1:ncol(x))
    {
        m<-m+w[i]*(outer(x[,i],x[,i],"*"))
    }
    m
}

hwu.weight.burden <- function(x, w = NULL)
{
    if(!is.matrix(x))
        stop('x is not a matrix')

    ## normalize all features (column wise)
    x <- apply(x, 2L, .map.std.norm)

    ## normalize feature weights
    if(is.null(w))
        w <- rep(1L, ncol(x))
    w <- w/sum(w)

    ## collapse features into one
    x %*% w
}

hwu.collapse.burden <- function(x)
{
    x<-as.matrix(x);
    w <- apply(x, 2, mean, na.rm = T)/2;
    w <- 1/sqrt(w*(1-w));               # may create NaN
    w[!is.finite(w)]<-0;
    w <- w/sum(w);
    g <- x %*% w;
    g
}

.wt.MAF <- function(x)
{
    if(!is.matrix(x))
        stop('x is not a matrix')

    ## get MAF
    m <- colMeans(x, na.rm=T) / 2
    w <- 1/sqrt(m * (1-m))
    #w[!is.finite(w)] <- 0
    w <- w/sum(w)
    w
}

.wt.MAFlg <- function(x)
{
    if(!is.matrix(x))
        stop('x is not a matrix')

    ## get MAF
    m <- colMeans(x, na.rm=T)
    w <- -log10(m)
    w[!is.finite(w)] <- 0
    w <- w/sum(w)
    w
}
