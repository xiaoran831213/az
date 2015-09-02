source('src/img.R')
source('src/utl.R')
source('src/hwu.R')
source('src/hlp.R')
source('src/sim_img.R')
library(igraph)

sim.vwa <- function(img, n.s = 100L, ft = 'tck', gb.sd = 1:7, seed = NULL)
{
    set.seed(seed)
    ## pick subjects
    sbj <- sample(img$sbj, n.s)
    img <- pck.sbj.img(img, sbj)
    n.s <- length(img$sbj)

    ## surface vertex from various level of gaussian bluring
    gb <- aperm(img$gsb[ft, , ,], perm=c('sbj', 'vtx', 'sdv'))
    
    ## assign effect to vertices. Lv.0  0 blur is the original
    ## use column major surface to make sure ve * ve will work.
    vt <- t(gb[, , 1L])
    ve.mu <- 0.0
    ve.sd <- 1.0
    ve.fr <- .05
    n.v <- length(img$vtx)
    ve = rnorm(n.v, ve.mu, ve.sd) * rbinom(n.v, 1L, ve.fr)
    
    ## ## vertex contributed phenotype
    z1 <- apply(ve * vt, 'sbj', sum)   # vertex effect * vertex value
    
    ## ## noise effect
    ns.rt <- 3.0                        # noise
    ne <- rnorm(n.s, 0, ns.rt * sd(z1))

    ## another phenotype is not affected by vertices
    z0 <- rnorm(n.s, 0, 1)
    set.seed(NULL)
    
    ## Derive U statistics, get P values of all encoding levels
    gb <- gb[, , gb.sd, drop = F]
    pv <- apply(gb, c('vtx', 'sdv'), function(v)
    {
        w <- .wct(.hwu.GUS(as.matrix(v)))
        c(
            p0=hwu.dg2(y=z0, w=w),
            p1=hwu.dg2(y=z1+ne, w=w))
    })
    names(dimnames(pv))[1] <- 'pvl'
    pv <- apply(pv, c('pvl', 'sdv'), min)
    nm <- expand.grid(dimnames(pv))
    pv <- unlist(pv)
    names(pv) <- do.call(paste, list(nm$sdv, nm$pvl, sep='.'))
    c(.record(), pv)
}

main.vwa <- function(n.itr = 10L, d.dat = .az.sm2, ...)
{
    fns <- pck.img(d.dat, size = n.itr, ret='file', seed = 150L)
    sim.rpt <- lapply(fns, function(fn)
    {
        img <- readRDS(fn)
        cat(fn, '\n')
        sim.vwa(img, ...)
    })
    HLP$mktab(sim.rpt)
}

