#### source('./GGRF.R')
library(Matrix)
#### library(genetics)
#### library(fda)
#### library(fda.usc)
library(CompQuadForm)
library(mgcv)

##trace function
tr<-function(A)
{
    return(sum(diag(A)))
}


##LD coefficient
getLD<-function(geno)
{
    geno.type.tmp<-apply(geno,2,as.genotype.allele.count,alleles=c("a","A"));
    geno.type<-as.data.frame(geno.type.tmp);
    geno.type<-makeGenotypes(geno.type);
    ld.tab<-LD(geno.type);
    ld.coef<-mean(ld.tab[[4]],na.rm=TRUE);
    return(ld.coef)
}

##Variance Component Score Test
VarScoreTest<-function(trait,geno,covariate=NULL,weights,type=c("normal","binary"))##,pos,maf)
{
    y<-trait;
    n<-length(y);
    if(is.null(covariate))
    {
        X<-rep(1,n);
    }else
    {
        X<-covariate;
    }

    ##knots<-pos;
    ##norder<-4;
    ##nbasis<-length(knots)+norder-2;
    ##genobasis<-create.bspline.basis(range(knots),nbasis,norder,knots);
    ##Lfdobj<-2; 
    ##genomaf<-t(apply(geno,1,function(x){x/sqrt(maf)}));
    ##genofd<-GetGenoFdn(nN,genomaf,pos,genobasis,Lfdobj,nbasis);
    ##rownames(genofd$data)<-c(1:nN);

    S<-getIBS(geno,weights)/(2*sum(weights));
    if(weights==1){S<-getIBS(geno,weights)/(2*ncol(geno))}
    ##S<-exp(-metric.lp(genofd,lp=1,w=1));
    diag(S)<-0;
    diag<-rowSums(S);
    D<-as.matrix(Diagonal(n,diag));
    ##gamma<-getLD(geno); ##tuning parameter gamma
    gamma<-mean(cor(geno));

    Va<-D-gamma*S;
    eig.Va<-eigen(Va,only.values=TRUE)$values;
    VaL<-chol(Va);
    inv.Va<-solve(VaL) %*% t(solve(VaL));

    ##Normal Trait
    if(type=="normal")
    {
        ##Estimates under H0
        beta0<-solve(t(X)%*%X)%*%t(X)%*%y;
        p<-length(beta0);
        sigma0<-(1/n)*t(y-X%*%beta0)%*%(y-X%*%beta0);
        
        ##Score Test Statistic
        P.hat<-X %*% solve(t(X) %*% X) %*% t(X);
        H<-diag(1-diag(P.hat));
        U_score<-(1/(2*sigma0^2)) * (t(y-X%*%beta0)%*%inv.Va%*%(y-X%*%beta0) - tr(H %*% inv.Va));

        ##variance of score test statistic
        ##I.Fisher<-(1/sigma0^2)*tr(Va%*%Va)-(1/(2*sigma0^2))*tr(inv.Va%*%inv.Va)-(1/(n*sigma0^2))*tr(Va)*(tr(Va)-(1/2)*tr(inv.Va));
        I.Fisher<-(1/(2*sigma0^2))*(sum(eig.Va^(-2))-(1/n)*(sum(eig.Va^(-1)))^2);

        ##Calculate p-Value
        ##p.value<-pchisq((U_score*I.Fisher^{-1/2})^2,df=1,lower=FALSE);
        p.value<-pnorm(U_score*I.Fisher^{-1/2},lower=FALSE);
        
    }

    if(type=="binary")
    {
        ##Necessary matrices and vectors estimated under H0
        mu0.hat<-glm(y~X,family=quasibinomial)$fitted.values;
        V.hat<-as.matrix(Diagonal(n,mu0.hat*(1-mu0.hat)));
        phi<-1/(n-1)*t(y-mu0.hat)%*%V.hat%*%(y-mu0.hat)
        k3<-phi^2*(1-2*mu0.hat)*mu0.hat*(1-mu0.hat);
        k4<-phi^3*(1-6*mu0.hat+6*mu0.hat^2)*mu0.hat*(1-mu0.hat);
        A<-inv.Va;
        r<-sqrt(2)*mu0.hat*(1-mu0.hat)/phi;
        R<-r%o%r;
        diag(R)<-phi^(-4)*k4+2*phi^(-2)*mu0.hat^2*(1-mu0.hat)^2;
        C<-as.matrix(Diagonal(n,phi^(-3)*k3));
        W<-as.matrix(Diagonal(n,phi^(-1)*mu0.hat*(1-mu0.hat)));
        ones<-rep(1,n);

        ##Efficient score
        U_score<-(1/(2*phi^2))*t(y-mu0.hat)%*%A%*%(y-mu0.hat)-(1/2)*tr(W%*%A);

        ##Fisher information
        I11<-(1/4)*t(ones)%*%(A%*%R%*%A)%*%ones;
        I21<-(1/2)*t(X)%*%C%*%diag(A);
        I22<-t(X)%*%W%*%X;
        I.Fisher<-I11-t(I21)%*%solve(I22)%*%I21;

        ##Calculate p-Value
        p.value<-pnorm(U_score*I.Fisher^{-1/2},lower=FALSE);
    }

    p.value
}
