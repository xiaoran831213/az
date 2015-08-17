source('src/gno.R')
source('src/utl.R')
source('src/hwu.R')
source('src/hlp.R')
source('src/sim_img.R')

## randomly pick encoded image data from a folder
mix.sim <- function(
    img, gno, N.S = .Machine$integer.max,
    ve.mu=0, ve.sd=.4, ve.fr=.05, vt.nm='tck',
    ge.mu=0, ge.sd=.4, ge.fr=.05,
    ep=.5, wp=.5, ne.rt=1.0)
{
    ## number of vertices and g-variants
    n.iv <- length(img$vtx)
    
    ## pick out subjects
    sbj <- intersect(img$sbj, gno$sbj)
    N.S <- min(N.S, length(sbj))
    sbj <- sample(sbj, N.S)
    img <- img.sbj.pck(img, sbj)
    gno <- gno.sbj.pck(gno, sbj)
    
    ## * -------- [vertex effect] -------- *
    ## for now we only use 1 vertex feature
    vt.ec <- subset(img$enc, grepl(vt.nm, names(img$enc)))

    ## encoding level 0 is in fact the unencoded vertices, make sure
    ## subjects are of column major so (ve * vt) could work!
    vt <- t(vt.ec[[1]])
    
    ## assign effect to some vertices
    ve = rnorm(n.iv, ve.mu, ve.sd) * rbinom(n.iv, 1L, ve.fr)
    
    ## vertex contributed phenotype
    y1.ve <- apply(ve * vt, 'sbj', mean)
    ## * -------- [vertex effect] -------- *


    ## * -------- [genome effect] -------- *
    gt <- gno$gmx                       # genomic matrix
    gt = GNO$clr.dgr(gt)                # clean degeneration
    if(length(gt) == 0)                 # give up the trial
        return(NA)
    gt = GNO$imp(gt)                    # imput missing g-variant
    n.gv <- nrow(gt)                    # survived g-variant

    ## assign effect to some genome
    ge <- rnorm(n.gv, ge.mu, ge.sd) * rbinom(n.gv, 1L, ge.fr)
    
    ## genome contributed response
    y1.ge <- apply(ge * gt, 'sbj', mean)
    ## * -------- [genome effect] -------- *

    ## sanity check: image and genome must share subjects
    stopifnot(colnames(gt) == colnames(vt))

    ## * -------- [joint effect] -------- *
    y1 <- y1.ve * ep + y1.ge * (1 - ep)
    y1.mu <- mean(y1)
    y1.sd <- sd(y1)
    
    ## noise effect
    ne <- rnorm(N.S, 0, ne.rt * y1.sd)

    ## a null response not affected by any variables
    y0 <- rnorm(N.S, y1.mu, y1.sd)

    ## * -------- U sta and P val --------*
    ## HWU requires subjes be of row major
    gt <- t(gt)
    
    ## apply weight proportions of genome and vertex
    vt <- lapply(vt.ec, function(v) v * wp)
    gt <- gt * (1 - wp)
    
    ## find genome weight, and vertex weight of all encoding depth
    wv <- lapply(vt, hwu.weight.gaussian)
    wg <- list(g.0=hwu.weight.IBS(gt))
    
    pv <- try(mapply(wv, wg, FUN = function(v, g)
    {
        list(
            #p0=hwu.dg2(y = y0 + ne, w=v),
            p1=hwu.dg2(y = y1 + ne, w=list(v, g)))
    }))

    if(inherits(pv, 'try-error'))
    {
        cat(gno$str, pv)
        return(NA)
    }
    
    c(.record(), pv)
}

## check is an object is a scalar
.scalar <- function(obj)
{
    if(is.list(obj))
        return(FALSE)
    if(!is.vector(obj))
        return(FALSE)
    if(!is.null(dim(obj)) || length(obj) > 1L)
        return(FALSE)
    TRUE
}

## collect object in a function environment, by default only
## visible scalars are collected
.record <- function(pass=.scalar)
{
    ret <- list()
    env <- parent.frame()
    for(nm in ls(env))
    {
        obj <- env[[nm]]
        if(!pass(obj))
            next
        ret[[nm]] <- obj
    }
    ret
}

mix.main <- function(n.itr = 5, n.sbj = 200)
{
    ## pick genotypes and images
    
    gno.dir <- seg.pck(vcf.dir = Sys.getenv('AZ_WGS'), size=n.itr, drop=F)
    img.dir <- img.pck(src = Sys.getenv('AZ_EC2'), size=n.itr, replace = T)
    sim.rpt <- mapply(img.dir, gno.dir, FUN = function(img.url, gno.seg)
    {
        img <- img.get(img.url)
        gno <- seg.get(gno.seg)
        mix.sim(img=img, gno=gno, N.S=200)
    }, SIMPLIFY = F)
    names(sim.rpt) <- sprintf('s%03X', 1L:length(sim.rpt))
    HLP$mktab(sim.rpt)
}
