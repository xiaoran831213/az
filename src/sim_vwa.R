source('src/img.R')
source('src/utl.R')
source('src/hwu.R')
source('src/hlp.R')
library(igraph)

sim.vwa <- function(img, n.s = 50L, ft = 'tck', seed = NULL)
{
    library(abind)
    set.seed(seed)
    
    ## pick subjects
    img <- img.sbj.pck(img, sample(img$sbj, n.s))
    n.s <- length(img$sbj)
    n.v <- length(img$vtx)

    ## surface feature of various level of gaussian blur
    gb <- abind(sdN=img$sfs[ft,,], img$gsb[ft,,,], along = 2L)
    names(dimnames(gb)) <- names(dimnames(img$gsb)[2:4])
    
    ## scale to [0, 1]
    gb1 <- apply(gb, 'sdv', function(g)
    {
        (g - min(g))/(max(g) - min(g))
    })
    
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
        w <- .hwu.GUS(e)
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

