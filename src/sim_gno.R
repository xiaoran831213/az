source('src/gno.R')
source('src/utl.R')
source('src/hwu.R')
source('src/hlp.R')
source('src/sim_img.R')

## randomly pick encoded image data from a folder
mix.sim <- function(
    img, gno, N.S = .Machine$integer.max,
    ve.mu=0, ve.sd=.4, ve.fr=.05,
    ge.mu=0, ge.sd=.4, ge.fr=.05,
    ne.rt=3)
{
    ## number of vertices and g-variants
    N.V <- length(img$vtx)
    N.G <- nrow(gno$map)
    
    ## pick out subjects
    N.S <- min(N.S, length(img$sbj), length(gno$sbj))
    img <- img.sbj.pck(img, sample(img$sbj, N.S))
    N.S <- length(img$sbj)
    
    ## for now we only use 1 feature
    ft.ec <- subset(img$enc, grepl(ft.nm, names(img$enc)))
    ft <- ft.ec[[1]]
    
    ## assign effect to some vertices
    ve = rnorm(N.V, ve.mu, ve.sd) * rbinom(N.V, 1L, ve.fr)
    
    ## vertex contributed phenotype
    y1.ve <- apply(ve * ft, 'sbj', mean)
    y1.ve.mu <- mean(y1.ve)
    y1.ve.sd <- sd(y1.ve)
    
    ## noise effect
    ne <- rnorm(N, 0, ne.rt * y1.ve.sd)

    ## a phenotype not affected by vertices nor genome
    y0 <- rnorm(N, y1.ve.mu, y1.ve.sd)

    ## Derive U statistics, get P values of all encoding levels
    pv.ec <- unlist(lapply(ft.ec, function(e)
    {
        wi <- hwu.weight.gaussian(e)
        list(
            p0=hwu.dg2(y = y0 + ne, w=wi),
            p1=hwu.dg2(y = y1.ve + ne, w=wi))
    }))
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

mix.sim <- function(gno, N.S = .N, ve.mu=0, ve.sd=.4, ve.fr=.05, ne.rt=3)
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

main <- function(n.itr = 5, n.sbj = 200)
{
    ## pick genotypes and images
    
    gno.dir <- seg.pck(.gen, .hkg, wnd=5000, size=n.itr)
    img.dir <- img.pck(Sys.getenv('AZ_EC2'), size=n.itr)
    sim.rpt <- mapply(img.dir, gno.dir, FUN = function(img.url, gno.seg)
    {
        img <- img.get(img.url)
        gno <- seg.get(gno.seg)
        list(img=img, gno=gno)
    }, SIMPLIFY = F)
    names(sim.rpt) <- sprintf('sim%03X', 1L:length(sim.rpt))
    ## HLP$mktab(sim.rpt)
    sim.rpt
}

