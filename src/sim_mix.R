source('src/gno.R')
source('src/utl.R')
source('src/hwu.R')
source('src/hlp.R')

## randomly pick encoded image data from a folder
img.pck <- function(src, N = NULL)
{
    if(!file.exists(src))
        stop(paste(src, "not exists"))

    ## pick image set (a white matter surface region)
    wms <- sample(x = dir(src, '*.rds'), size = 1)
    img <- readRDS(paste(src, wms, sep='/'))
    img$wms <- wms
    img$src <- src
    
    ## pick N subjects
    if (is.null(N))
        N <- .Machine$integer.max
    N <- min(N, dim(img$sfs)[3L])
    idx <- sort(sample.int(n = N))

    img <- within(
        img,
        {
            sfs <- sfs[, , idx, drop = F];
            enc <- lapply(enc, function(x) x[idx, , drop = F]);
            N <- N;
            M <- dim(sfs)[2L];
        })

    dmn <- dimnames(img$sfs);
    names(dmn) <- c('ftr', 'vtx', 'sbj')
    dimnames(img$sfs) <- dmn
    return(img)
}

## rescale every feature to [0,1]
img.sc1 <- function(x, MARGIN = NULL)
{
    ## basic rescalar
    sc1 <- function(u) (u-min(u))/(max(u)-min(u))
    
    ## deal with trivial scenario, e.g. vectors, null margin.
    if(is.null(dim(x)) || is.null(MARGIN))
        return(sc1(x))

    d <- dim(x)
    m <- sort(MARGIN)
    
    ## permute the margins to higher dimensions
    j <- length(d) - length(m) + 1
    perm <- 1L:length(d)
    for(i in m)
    {
        perm[i] <- j; perm[j] <- i
        j <- j + 1
    }
    x <- aperm(x, perm)
    m <- (length(d)-length(m) + 1) : length(d)
    
    ## apply operation to lower dimensions, here the lower
    ## dimensions will collapse
    y <- apply(x, m, sc1)

    ## restore lower dimensions
    dim(y) <- dim(x)
    dimnames(y) <- dimnames(x)

    ## permute the dimensions back
    y <- aperm(y, perm)
    y
}

img.sim <- function()
{
    ## randomly pick WM surface sample
    img <- img.pck(Sys.getenv('AZ_EC1'), 200L)
    N <- img$N  # actual subject count
    M <- img$M  # vertex count
    enc <- img$enc

    ## for now we only use thickness feature, rescaled to [0, 1]
    f1.nm <- 'tck'
    f1 <- img.sc1(img$sfs[f1.nm, , ])
    f1.ec <- c(
        list(t(f1)),
        subset(enc, grepl(f1.nm, names(enc))))
    names(f1.ec)[1L] <- paste(f1.nm, 0, sep='.')
    
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

    ## a correlated covariant
    c1.mr <- 1.0                        # c1 to vertex mu ratio
    c1.sr <- 1.0                        # c1 to vertex sd ratio
    c1 <- rnorm(n = N, mean = c1.mr * z1.mu, sd = c1.sr * sd(z1))

    ## Derive U statistics, get P values of all encoding levels
    pv.ec <- unlist(lapply(f1.ec, function(e)
    {
        wi <- hwu.weight.gaussian(e)
        list(
            p0=hwu.dg2(y=z0+ne, x=c1, w=wi),
            p1=hwu.dg2(y=z1+ne, x=c1, w=wi))
    }))
    
    ## p value of per-vertex test. to conserve time, only functional
    ## vertices are picked
    ## permute dimensions to: (subject, feature, vertex)
    ## pv <- apply(aperm(x[ve.mk, ,], c(3, 2, 1)), 3L, function(v)
    ## {
    ##     m <- glm.fit(x = cbind(v, c1), y = z1)
    ##     s <- summary(m)
    ## })
    ## simulation record
    c(.record(), pv.ec)
}

## check is an object is a scalar
.scalar <- function(obj)
{
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

mix.sim <- function(N = NULL)
{
    ## pick encoded image data
    img = img.pck('hpc/encoded_wms/09_0078')
    imx <- img$imx
    
    ## pick genomic segment, and only pick those with image available
    gno <- GNO$pck(vcf = 'dat/wgs/gno', seg='dat/wgs/gen', sbj=img$sbj, wnd=5000, n = 1)[[1]]
    gmx <- gno$gmx
    
    ## randomly select common subject
    sbj = intersect(gno$sbj, img$sbj)
    if(!is.null(N))
    {
        N <- min(N, length(sbj))
    }
    else
    {
        N <- length(sbj)
    }
    sbj = sort(sbj[sample.int(length(sbj), size = N)])
    imx <- imx[, match(sbj, img$sbj)]
    gmx <- gmx[, match(sbj, gno$sbj)]
    
    ## clean up genotype
    gmx = GNO$clr.dgr(gmx)
    ## sometimes cleaup will uncheck all variants, we have
    ## to skip this type I trial
    if(length(gmx) == 0)
        return(NA)

    ## guess missing values
    gmx = GNO$imp(gmx)
    
    ## tranpose data matrices to design matrices
    gmx = t(gmx)
    imx = t(imx)
    stopifnot(nrow(gmx) == nrow(imx))
    
    ## make one binary phenotype and normal covariate
    y <- rbinom(n = N, size = 1, prob = 0.5)
    x <- rnorm(n = N, mean = 0, sd = 1.5)
    o <- list()
    
    ## g <- hwu.collapse.burden(gmx)
    ## wg <- hwu.weight.gaussian(gmx)
    wg <- hwu.weight.IBS(gmx)
    
    ## i <- hwu.collapse.burden(imx)
    ## wi <- hwu.weight.gaussian(imx)
    wi <- hwu.weight.gaussian(imx)

    o$g1IBS.i1GSN <- hwu.dg2(y, x, wg, wi)
    
    wg <- hwu.weight.IBS(gmx, w = hwu.w.MAFsd(gmx))
    o$gWIBS.i1GSN <- hwu.dg2(y, x, wg, wi)

    o
}

main <- function(n.itr = 1000)
{
    p <- replicate(n.itr, img.sim(), simplify = F)
    HLP$mktab(p)
}