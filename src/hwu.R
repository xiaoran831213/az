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
hwu.dg2 <- function(y, x=NULL, w, ...)
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
    r <- diag(1, M, M) - tcrossprod(x %*% solve(crossprod(x)), x);
    
    ## exclude liner covariant effect on y, leave residual of Y
    y <- r %*% y;
    y <- y/sqrt(sum(y^2)/(M-ncol(x)));
 
    ## the U kernel is the pair wise similarity between phenotypes
    f <- tcrossprod(y);

    ## get product of all weight terms.
    w <- Reduce(f='*', x=list(w, ...));
    diag(w) <- 0; # ??
    
    ## compute U score
    u <- sum(w * f);

    ## exclude coveriant and intercept effect on both
    ## dimensions of w
    w <- tcrossprod(r %*% w, r);

    ## calculate p-value of u.
    coef <- eigen(w, symmetric=T, only.values=T)$values;
    p <- davies(u, coef, acc=0.000001)$Qq;
    p
}

.map.std.norm<-function(y)
{
    y<-rank(y);
    y<-(y-0.5)/length(y);
    qnorm(y)
}

.map.std<-function(y)
{
    y<-y-mean(y)
    sdy<-sd(y)
    if(sdy!=0) y<-y/sdy;
    y
}

hwu.weight.gaussian <- function(x, w = NULL)
{
    stopifnot(is.matrix(x))

    ## normalize features
    x <- apply(x, 2, .map.std.norm);
    if(is.null(w))
        w<-rep(1,ncol(x));
    w<-w / sum(w) / 2;

    ## subject pairwise measure
    m <- matrix(0,nrow=nrow(x),ncol=nrow(x))

    ## go through all feature to measure gaussian distance
    for(i in 1:ncol(x))
    {
        m <- m + w[i] * (outer(x[,i],x[,i],"-")) ^ 2
    }
    
    ## exp(- gaussian distance) = gaussian similiarity
    exp(-m)
}

hwu.weight.IBS <- function(x, w = NULL, lv = 2L)
{
    stopifnot(is.matrix(x))

    ## normalize feature weights
    if(is.null(w))
        w <- rep(1L, ncol(x))
    w <- w / sum(w) / lv
    
    ## subject pairwise similarity weights
    m <- matrix(0, nrow = nrow(x), ncol = nrow(x))

    ## go through all featurem to measure similarity
    for(i in 1L:ncol(x))
    {
        m <- m + w[i] * (lv - abs(outer(x[,i], x[,i], '-')))
    }
    m
}

hwu.weight.cov<-function(x, w = NULL)
{
    stopifnot(is.matrix(x))

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
    stopifnot(is.matrix(x))

    ## normalize all features (column wise)
    x <- apply(x, 2L, .map.std.norm)

    ## normalize feature weights
    if(is.null(w))
        w <- rep(1L, ncol(x))
    w <- w/sum(w)

    ## collapse features into one
    crossprod(x, w)
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

hwu.w.MAFsd <- function(x)
{
    stopifnot(is.matrix(x))

    ## get MAF
    m <- colMeans(x, na.rm=T) / 2
    w <- 1/sqrt(m * (1-m))
    #w[!is.finite(w)] <- 0
    w <- w/sum(w)
    w
}

hwu.w.MAFlg <- function(x)
{
    stopifnot(is.matrix(x))

    ## get MAF
    m <- colMeans(x, na.rm=T)
    w <- -log10(m)
    w[!is.finite(w)] <- 0
    w <- w/sum(w)
    w
}
