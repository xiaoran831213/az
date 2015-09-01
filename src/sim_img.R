source('src/utl.R')
source('src/hwu.R')
source('src/hlp.R')
source('src/img.R')
library(Matrix)

ini.img <- function(img)
{
    ## append dimension names to the surface data
    names(dimnames(img$sfs)) <- c('ftr', 'vtx', 'sbj')

    ## basic decoration
    img <- within(
        img,
    {
        names(dimnames(cmx)) <- list('a', 'b')
        sbj <- dimnames(sfs)$sbj
        vtx <- dimnames(sfs)$vtx
        enc <- lapply(enc, function(u)
        {
            dimnames(u) <- list(sbj=sbj, vtx=sprintf('C%04X', 1L:ncol(u)))
            u
        })
    })
    img
}

pck.sbj.img <- function(img, sbj)
{
    I <- match(sbj, img$sbj)
    img <- within(
        img,
    {
        sbj <- sbj[I]
        sfs <- sfs[, , I]
        enc <- lapply(enc, function(u)
        {
            u[I, ]
        })
    })
}

img.sim <- function(img, n.s = 50L, ft = 'tck', seed = NULL)
{
    set.seed(seed)
    
    ## pick subjects
    img <- pck.sbj.img(img, sample(img$sbj, n.s))
    
    n.s <- length(img$sbj)
    n.v <- length(img$vtx)
    enc <- img$enc

    ## surface feature encoding
    vc <- subset(enc, grepl(ft, names(enc)))
    names(vc) <- paste('e', 1L:length(vc) - 1L, sep = '')
    vt <- t(vc[[1]])                 # vertex at column major

    ## assign effect to vertices
    ve.mu <- 0.0
    ve.sd <- 1.0
    ve.fr <- .05
    ve = rnorm(n.v, ve.mu, ve.sd) * rbinom(n.v, 1L, ve.fr)
    
    ## vertex contributed phenotype
    z1 <- apply(ve * vt, 'sbj', sum)   # vertex effect * vertex value
    
    ## noise effect
    ns.rt <- 3.0                        # noise ratio
    ne <- rnorm(n = n.s, mean = 0, sd = ns.rt * sd(z1))

    ## another phenotype is not affected by vertices
    z0 <- rnorm(n = n.s, mean = 0, sd = 1)

    ## Derive U statistics, get P values of all encoding levels
    pv <- lapply(vc, function(e)
    {
        w <- .wct(.hwu.GUS(e))
        list(
            p0=hwu.dg2(y=z0, w=w),
            p1=hwu.dg2(y=z1+ne, w=w))
    })

    ## resume R random stream
    set.seed(NULL)
    
    c(.record(), unlist(pv))
}

img.main <- function(n.itr = 10L, d.dat = .az.sm2, ...)
{
    fns <- pck.img(d.dat, size = n.itr, replace = T, ret='file')
    sim.rpt <- lapply(fns, function(fn)
    {
        img <- readRDS(fn)
        img <- ini.img(img)
        cat(fn, '\n')
        img.sim(img, ...)
    })
    HLP$mktab(sim.rpt)
}
.az.ec2 <- Sys.getenv('AZ_EC2')
.az.ec3 <- Sys.getenv('AZ_EC3')
.az.ec4 <- Sys.getenv('AZ_EC4')
.az.ec5 <- Sys.getenv('AZ_EC5')

img.test <- function(n.s = 200)
{
    n.i <- 2000

    t2 <- img.main(n.i, n.s=n.s, d.dat=.az.ec2, seed = NULL)
    t2
    ## t4 <- img.main(n.i, n.s=n.s, d.dat=.az.ec4, seed = NULL)
    ## t5 <- img.main(n.i, n.s=n.s, d.dat=.az.ec5, seed = NULL)
    ## t2$ec <- '-1/2'
    ## t4$ec <- '-3/4'
    ## t5$ec <- '-2/3'
    ## rt <- rbind(t2, t4, t5)
    ## rt
}

img.pwr1 <- function(rpt, t = 0.05, ret = 2)
{
    if(ret == 0)
        pvl.rgx <- 'p0$'
    else if(ret == 1)
        pvl.rgx <- 'p1$'
    else
        pvl.rgx <- 'p[01]$'

    pvl <- rpt[, grep(pvl.rgx, names(rpt))]
    cfg <- rpt[, grep('p[01]', names(rpt), invert = T)]
    
    rej <- function(p)
    {
        sum(p < t, na.rm = T)/sum(!is.na(p))
    }

    n.i <- aggregate(cfg$ft, by=cfg, length, simplify = T)
    n.i <- n.i[, ncol(n.i)]
    p.v <- aggregate(pvl, by=cfg, rej)
    
    cbind(n.i, p.v)
}
