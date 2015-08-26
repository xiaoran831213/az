source('src/utl.R')
source('src/hwu.R')
source('src/hlp.R')

.img.read <- function(fn, verbose = FALSE)
{
    img <- readRDS(fn)
    if(verbose)
        cat(fn, '\n')
    
    ## append dimension names to the surface data
    names(dimnames(img$sfs)) <- c('ftr', 'vtx', 'sbj');

    within(
        img,
    {
        ## append dimension names to the encodings
        names(dimnames(cmx)) <- list('a', 'b')
        src <- dirname(fn)              # source folder
        sbj <- dimnames(sfs)$sbj
        vtx <- dimnames(sfs)$vtx
        enc <- lapply(enc, function(u)
        {
            dimnames(u) <- list(sbj=sbj, vtx=sprintf('C%04X', 1L:ncol(u)))
            u
        })
        ssn <- sub('[.]rds', '', basename(fn)) # center vertex
    })
}

## randomly pick encoded image data from a folder
img.pck <- function(
    src, size = 1, replace = FALSE, drop = TRUE, vbs = FALSE, ret = c('data', 'file'))
{
    ## pick out images by file name
    fns <- file.path(src, dir(src, '*.rds'))
    if(replace | size < length(fns))
        fns <- sample(fns, size, replace)
    
    ## only return file nemas
    if(ret[1] == 'file')
        return(fns)
    
    ret <- sapply(fns, .img.read, verbos = vbs, simplify = F, USE.NAMES = F)
    names(ret) <- sub('[.]rds', '', basename(fns))

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

.az.img <- Sys.getenv('AZ_EC2')         # 1/2 encoding
.az.ec3 <- Sys.getenv('AZ_EC3')         # super fitted encoding
.az.ec4 <- Sys.getenv('AZ_EC4')         # 3/4 encoding
.az.ec5 <- Sys.getenv('AZ_EC5')         # 2/3 encoding
img.main <- function(n.itr = 10L, n.sbj = 200L, v.dat = NULL, d.dat = .az.img)
{
    if(is.null(v.dat))
    {
        fns <- img.pck(d.dat, size = n.itr, ret='file')
        sim.rpt <- sapply(fns, function(fn)
        {
            img <- .img.read(fn, verbose = T)
            img.sim(img, n.s=n.sbj)
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

img.pwr <- function(rpt, t = 0.05, ret = 2)
{
    n.itr <- nrow(rpt)
    if(ret == 0)
        rgx <- 'p[0-9]*[.]0$'
    else if(ret == 1)
        rgx <- 'p[0-9]*[.]1$'
    else
        rgx <- 'p[0-9]*[.][01]$'
        
    p.hdr <- grepl(rgx, colnames(rpt))
    p.val <- subset(rpt, select=p.hdr)
    lapply(p.val, function(p) sum(p < t) / n.itr)
}
