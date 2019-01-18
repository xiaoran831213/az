library(Matrix)
library(CompQuadForm)
library(mgcv)

## trace function
tr<-function(A)
{
    return(sum(diag(A)))
}

VarScoreTest<-function(trait,geno,covariate=NULL,weights=1)
{
    y.tilde <- trait;
    n <- length(y.tilde)
    if(is.null(covariate))
    {
        X <- as.matrix(rep(1,n));
    }else
    {
        X <- covariate;
    }

    S<-getIBS(geno,weights)
    if(weights==1)
    {
        S <- getIBS(geno,weights)/(2*ncol(geno))
    }
    ## S<-getDPS(geno,weights, od=2)
    ## if(weights==1)
    ## {
    ##     S <- getDPS(geno,weights, od=2)
    ## }

    diag(S)<-0;
    diag<-rowSums(S);
    D<-as.matrix(Diagonal(n,diag));
    gamma<-mean(cor(geno));

    Va<-D-gamma*S;
    
    VaL<-chol(Va);
    inv.Va<-chol2inv(VaL);

                                        #transform y using QR decomposition
    qrX <- qr(X)
    y <- qr.qty(qrX,y.tilde)[(rankMatrix(X)+1):n];
    K <-t(qr.Q(qrX,complete=TRUE)[,(rankMatrix(X)+1):n]);

                                        #observed score statistic
    s1 <- ((n-rankMatrix(X))/2) * t(y) %*% K %*% inv.Va %*% t(K) %*% y / (t(y) %*% y) - (1/2) * tr(inv.Va)

    const0 <- (1 / (n - rankMatrix(X))) * (2 * s1 + tr(inv.Va));
    const0 <- as.numeric(const0);
    B0 <- K %*% inv.Va %*% t(K) - const0 * diag(rep(1,n - rankMatrix(X)));	
    
    p.value <- davies(0,eigen(B0,symmetric=TRUE,only.value=TRUE)$values)$Qq;;
    p.value
}
