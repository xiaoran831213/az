source('src/utl.R')
source('src/hwu.R')
source('src/hlp.R')

## randomly pick encoded image data from a folder
img.pck <- function(src, size = 1, replace = FALSE, drop = TRUE, vbs = FALSE)
{
    ## pick image set (a white matter surface region)
    fns <- sample(dir(src, '*.rds', full.names = T), size, replace)
    ids <- sub('[.]rds', '', basename(fns))

    ret <- mapply(fns, ids, FUN = function(fn, id)
    {
        img <- readRDS(fn)
        if(vbs)
            cat(fn, '\n')

        ## append dimension names to the surface data
        names(dimnames(img$sfs)) <- c('ftr', 'vtx', 'sbj');

        within(
            img,
        {
            ## append dimension names to the encodings
            src <- src
            sbj <- dimnames(sfs)$sbj
            vtx <- dimnames(sfs)$vtx
            enc <- lapply(enc, function(u)
            {
                dimnames(u) <- list(sbj=sbj, cdx=sprintf('C%04X', 1L:ncol(u)))
                u
            })
            ssn <- id
        })
    }, SIMPLIFY = FALSE)
    names(ret) <- ids

    if(drop & length(ret) < 2L)
        ret <- ret[[1]]
    ret
}

img.sbj.pck <- function(img, sbj)
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

    ## in case sb. think {sbj} means sample size instead of
    ## the IDs of wanted subject
    if(length(img$sbj) == 0L)
        warning('no subject ID matches image data source.')
    img
}

img.sim <- function(img, n.s = 200L, f1.nm = 'tck')
{
    ## pick subjects
    img <- img.sbj.pck(img, sample(img$sbj, n.s))
    
    N <- length(img$sbj)
    M <- length(img$vtx)
    enc <- img$enc

    ## for now we only use 1 feature, also rescaled it to [0, 1]
    f1.nm <- 'tck'
    f1.ec <- subset(enc, grepl(f1.nm, names(enc)))
    f1 <- t(f1.ec[[1]])                 # vertex at column major

    ## assign effect to each vertex
    ve.mu <- 0.0
    ve.sd <- 0.4
    ve = rnorm(n = M, mean = ve.mu, sd = ve.sd)
    
    ## mask functional vertices
    ve.fr <- 0.05
    ve.mk <- rbinom(n = M, size = 1, prob = ve.fr)
    ve <- ve * ve.mk                    # vertex effect * vertex mask
    
    ## vertex contributed phenotype
    z1 <- apply(ve * f1, 'sbj', mean)   # vertex effect * vertex value
    
    z1.mu <- mean(z1)
    z1.sd <- sd(z1)
    
    ## noise effect
    ne.sr <- 3.0                        # noise to vertex sd ratio
    ne <- rnorm(n = N, mean = 0, sd = ne.sr * z1.sd)

    ## another phenotype is not affected by vertices
    z0 <- rnorm(n = N, mean = z1.mu, sd = z1.sd)

    ## Derive U statistics, get P values of all encoding levels
    pv.ec <- lapply(f1.ec, function(e)
    {
        w <- .hwu.GUS(e)
        list(
            p.0=hwu.dg2(y=z0+ne, w=w),
            p.1=hwu.dg2(y=z1+ne, w=w))
    })
    c(.record(), unlist(pv.ec))
}

.az.img <- Sys.getenv('AZ_EC2')
img.main <- function(n.itr = 10L, n.sbj = 200L, v.dat = NULL)
{
    if(is.null(v.dat))
    {
        sim.rpt <- replicate(n.itr,
        {
            v <- img.pck(.az.img, vbs = T)
            img.sim(v, n.s=n.sbj)
        }, simplify = FALSE)
    }
    else
    {
        sim.rpt <- sapply(v.dat, function(v)
        {
            cat(v$ssn, '\n')
            img.sim(v, n.s=n.sbj)
        }, simplify = FALSE)
        rm(v.dat)
    }
    
    HLP$mktab(sim.rpt)
}

img.pwr <- function(rpt, t = 0.05)
{
    n.itr <- nrow(rpt)
    p.hdr <- grepl('p[0-9]*[.][01]$', colnames(rpt))
    p.val <- subset(rpt, select=p.hdr)
    lapply(p.val, function(p) sum(p < t) / n.itr)
}

img.pw0 <- function(rpt, t = 0.05)
{
    n.itr <- nrow(rpt)
    p.hdr <- grepl('p[0-9]*[.]0$', colnames(rpt))
    p.val <- subset(rpt, select=p.hdr)
    lapply(p.val, function(p) sum(p < t) / n.itr)
}

img.pw1 <- function(rpt, t = 0.05)
{
    n.itr <- nrow(rpt)
    p.hdr <- grepl('p[0-9]*[.]1$', colnames(rpt))
    p.val <- subset(rpt, select=p.hdr)
    lapply(p.val, function(p) sum(p < t) / n.itr)
}
