source('src/utl.R')
source('src/hwu.R')
source('src/hlp.R')
library(Matrix)

img.sim <- function(img, n.s = 50L, ft = 'tck', seed = NULL)
{
    set.seed(seed)
    
    ## pick subjects
    img <- img.sbj.pck(img, sample(img$sbj, n.s))
    
    n.s <- length(img$sbj)
    n.v <- length(img$vtx)
    enc <- img$enc

    ## surface feature encoding
    vc <- subset(enc, grepl(ft, names(enc)))
    names(vc) <- paste('e', 1L:length(vc) - 1L, sep = '')
    vt <- t(vc[[1]])                 # vertex at column major

    ## assign effect to vertices
    ## ve.mu <- 0.0
    ## ve.sd <- 1.0
    ## ve.fr <- .05
    ## ve = rnorm(n.v, ve.mu, ve.sd) * rbinom(n.v, 1L, ve.fr)
    
    ## ## vertex contributed phenotype
    ## z1 <- apply(ve * vt, 'sbj', mean)   # vertex effect * vertex value
    
    ## ## noise effect
    ## ns.rt <- 4.0                        # noise
    ## ne <- rnorm(n = n.s, mean = 0, sd = ns.rt * sd(z1))

    ## another phenotype is not affected by vertices
    z0 <- rnorm(n = n.s, mean = 0, sd = 2)

    ## Derive U statistics, get P values of all encoding levels
    pv <- lapply(vc, function(e)
    {
        w <- .wct(.hwu.GUS(e))
        list(
            p0=hwu.dg2(y=z0, w=w))
            #p1=hwu.dg2(y=z1+ne, w=w))
    })

    ## resume R random stream
    set.seed(NULL)
    
    c(.record(), unlist(pv))
}

img.main <- function(n.itr = 10L, d.dat = .az.img, ...)
{
    fns <- img.pck(d.dat, size = n.itr, ret='file', seed = 150L)
    sim.rpt <- lapply(fns, function(fn)
    {
        img <- .img.read(fn, vbs = T)
        img.sim(img, ...)
    })
    HLP$mktab(sim.rpt)
}

img.test <- function(n.s = 200)
{
    n.i <- 1000

    t2 <- img.main(n.i, n.s=n.s, d.dat=.az.ec2)
    t4 <- img.main(n.i, n.s=n.s, d.dat=.az.ec4)
    t5 <- img.main(n.i, n.s=n.s, d.dat=.az.ec5)
    t2$ec <- '-1/2'
    t4$ec <- '-3/4'
    t5$ec <- '-2/3'
    rt <- rbind(t2, t4, t5)
    rt
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
