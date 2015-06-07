library(CompQuadForm)
HWU<-new.env();

HWU$core<-function(y,K,geno,X=NULL,center.geno=F,gsim=c("add","eq","dist"),appx=c("davis","normal"))
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
        mt.geno=HWU$weight.gaussian(geno);
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

HWU$run<-function(Tn,K,geno,X=NULL,center.geno=F,gsim=c("add","eq","dist"),appx=c("davis","normal"))
{
    appx<-match.arg(appx);
    gsim<-match.arg(gsim);
    HWU$core(Tn,K,geno,X,center.geno,gsim=gsim,appx=appx);
}

## x ---- covariate
## y ---- response variable
## f ---- U kernel
## r ---- residual matrix
## w, ... ---- weight terms.
HWU$dg2 <- function(y, x=NULL, w, ...)
{
    ## response and covariate
    M <- length(y);

    ## standardize y to m=0, s=1
    y <- rank(y);
    y <- (y-mean(y))/sd(y);
    
    ## regression residual matrix, R = I - X(X'X)^X'
    if(is.null(x))
    {
        x = rep(1, M)
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

HWU$map.std.norm<-function(y)
{
    y<-rank(y);
    y<-(y-0.5)/length(y);
    qnorm(y)
}

HWU$map.std<-function(y)
{
    y<-y-mean(y)
    sdy<-sd(y)
    if(sdy!=0) y<-y/sdy;
    y
}

HWU$weight.gaussian<-function(x,weight=NULL)
{
    x<-apply(x,2, HWU$map.std.norm);
    if(is.null(weight)) weight<-rep(1,ncol(x));
    weight<-weight/sum(weight)/2;
    
    wt<-matrix(0,nrow=nrow(x),ncol=nrow(x));
    for(i in 1:ncol(x))
    {
        wt<-wt-weight[i]*(outer(x[,i],x[,i],"-"))^2
    }
    #without exp is also okay
    wt<-exp(wt);
    wt
}

HWU$weight.cov<-function(x,weight=NULL)
{
    x<-apply(x,2,map.std.norm);
    if(is.null(weight)) weight<-rep(1,ncol(x));
    weight<-weight/sum(weight);
    
    wt<-matrix(0,nrow=nrow(x),ncol=nrow(x));
    for(i in 1:ncol(x))
    {
        wt<-wt+weight[i]*(outer(x[,i],x[,i],"*"))
    }
    wt
}

HWU$collapse.burden<-function(x)
{
    x<-as.matrix(x);
    w <- apply(x, 2, mean, na.rm = T)/2;
    w <- 1/sqrt(w*(1-w));               # may create NaN
    w[!is.finite(w)]<-0;
    w <- w/sum(w);
    g <- x %*% w;
    g
}
