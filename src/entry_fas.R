source('src/gno.R')
source('src/utl.R')
source('src/hwu.R')
source('src/hlp.R')

## randomly pick encoded image data from a folder
img.pck <- function(src, N = NULL)
{
    if(!file.exists(src))
        stop(paste(src, "not exists"))

    img <- sample(x = dir(src, '*.rds'), size = 1)
    img <- readRDS(paste(src, img, sep='/'))

    ## tranfrom possible list structure to pure array
    x <- unlist(img$x)
    y <- unlist(img$y)
    dim(x) <- dim(img$x)
    dim(y) <- dim(img$y)
    dimnames(x) <- dimnames(img$x)
    dimnames(y) <- dimnames(img$y)
    
    ## pick subjects
    if (is.null(N))
        N <- dim(x)[3]
    if (N < dim(x)[3])
    {
        idx <- sample.int(n = N)
        x <- x[, , idx, drop = F]
        y <- y[idx, , drop = F]
    }

    img$x <- x
    img$y <- y
    img$N <- N
    img
}

## rescale every feature to [0,1]
img.sc1 <- function(x, MARGIN)
{
    ## permute the margins to higher dimensions
    m <- sort(MARGIN)
    d <- dim(x)
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
    y <- apply(x, m, function(u)
    {
        (u-min(u))/(max(u)-min(u))
    })
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
    N <- 200L
    img <- img.pck('hpc/trained_sda/09_0078', N = N)
    x <- img$x    # image data
    y <- img$y    # image code
    N <- img$N  # actual subject count

    ## rescale to [0,1] over each feature
    x <- img.sc1(x, 2L)

    ## assign effect to each vertex feature(s)
    ve.mu <- 0.0
    ve.sd <- 0.4
    ve.fr <- 0.05
    ve = mapply(
        rnorm, n = dim(x)[1L],
        mean = rep(ve.mu, times = dim(x)[2L]),    # effect mean(s)
        sd = rep(ve.sd, times = dim(x)[2L]))      # effect sd(s)
    ## mask functional surface vertices
    ve.mk <- rbinom(n = dim(x)[1L], size = 1, prob = ve.fr)

    ## calculate subject phenotype affected by vertices
    ve <- apply(ve, 2L, function(fe)
    {
        fe * ve.mk    # feture effect * vertex mask
    })
    z1 <- apply(x, 3L, function(sv)
    {
        mean(sv * ve)                   # subject vertex * vertex effect
    })
    z1.mu <- mean(z1)
    z1.sd <- sd(z1)

    ## noise effect
    ne.sr <- 1.0                      # noise to vertex sdev ratio
    ne <- rnorm(n = N, mean = 0, sd = ne.sr * z1.sd)

    ## another phenotype is not affected by vertices
    z0 <- rnorm(n = N, mean = z1.mu, sd = z1.sd)

    ## a correlated covariant
    c1.mr <- 0.3                      # c1 to vertex mean ratio
    c1.sr <- 1.5                      # c1 to vertex sdev ratio
    c1 <- rnorm(n = N, mean = c1.mr * z1.mu, sd = c1.sr * sd(z1))

    ## Derive U statistics

    ## p value of image encoding test
    wi <- hwu.weight.gaussian(y)
    iy.p0 <- hwu.dg2(y=z0+ne, x=c1, w=wi)
    iy.p1 <- hwu.dg2(y=z1+ne, x=c1, w=wi)

    ## p value of image raw feature test (golden standard)
    x3 <- aperm(x, c(3, 1, 2))
    dim(x3) <- c(N, length(x3) / N)
    wi <- hwu.weight.gaussian(x3)
    ix.p0 <- hwu.dg2(y =z0+ne, x=c1, w=wi)
    ix.p1 <- hwu.dg2(y =z1+ne, x=c1, w=wi)

    ## p value of per-vertex test. to conserve time, only functional
    ## vertices are picked
    ## permute dimensions to: (subject, feature, vertex)
    ## pv <- apply(aperm(x[ve.mk, ,], c(3, 2, 1)), 3L, function(v)
    ## {
    ##     m <- glm.fit(x = cbind(v, c1), y = z1)
    ##     s <- summary(m)
    ## })
    ## simulation record
    .record()
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
