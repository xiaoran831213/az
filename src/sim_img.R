source('src/utl.R')
source('src/hwu.R')
source('src/hlp.R')
library(Matrix)

.img.read <- function(fn, vbs = FALSE)
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
    src, size = 1, replace = FALSE, seed = NULL,
    drop = TRUE, vbs = FALSE, ret = c('data', 'file'))
{
    set.seed(seed)
    
    ## pick out images by file name
    fns <- file.path(src, dir(src, '*.rds'))
    if(replace | size < length(fns))
        fns <- sample(fns, size, replace)
    
    ## only return file nemas
    if(ret[1] == 'file')
        return(fns)
    
    ret <- sapply(fns, .img.read, vbs = vbs, simplify = F, USE.NAMES = F)
    names(ret) <- sub('[.]rds', '', basename(fns))

    if(drop & length(ret) < 2L)
        ret <- ret[[1]]

    set.seed(NULL)
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
    ve.mu <- 0.0
    ve.sd <- 1.0
    ve.fr <- .05
    ve = rnorm(n.v, ve.mu, ve.sd) * rbinom(n.v, 1L, ve.fr)
    
    ## vertex contributed phenotype
    z1 <- apply(ve * vt, 'sbj', mean)   # vertex effect * vertex value
    
    ## noise effect
    ns.rt <- 3.0                        # noise
    ne <- rnorm(n = n.s, mean = 0, sd = ns.rt * sd(z1))

    ## another phenotype is not affected by vertices
    z0 <- rnorm(n = n.s, mean = 0, sd = 1)

    ## Derive U statistics, get P values of all encoding levels
    pv <- lapply(vc, function(e)
    {
        w <- .hwu.GUS(e)
        list(
            p0=hwu.dg2(y=z0+ne, w=w),
            p1=hwu.dg2(y=z1+ne, w=w))
    })

    ## resume R random stream
    set.seed(NULL)
    
    c(.record(), unlist(pv))
}

.az.ec2 <- Sys.getenv('AZ_EC2')         # 1/2 encoding
.az.ec3 <- Sys.getenv('AZ_EC3')         # super fitted encoding
.az.ec4 <- Sys.getenv('AZ_EC4')         # 3/4 encoding
.az.ec5 <- Sys.getenv('AZ_EC5')         # 2/3 encoding
img.main <- function(n.itr = 10L, n.sbj = 200, d.dat = .az.img)
{
    fns <- img.pck(d.dat, size = n.itr, ret='file', seed = 150L)
    sim.rpt <- lapply(fns, function(fn)
    {
        img <- .img.read(fn, vbs = T)
        img.sim(img, n.s = n.sbj, seed=120L)
    })
    HLP$mktab(sim.rpt)
}

img.test <- function()
{
    n.i <- 1000

    n.s <- 50
    t2 <- img.main(n.i, n.s, d.dat=.az.ec2)
    t3 <- img.main(n.i, n.s, d.dat=.az.ec3)
    t4 <- img.main(n.i, n.s, d.dat=.az.ec4)
    t5 <- img.main(n.i, n.s, d.dat=.az.ec5)
    t2$ec <- '-1/2'
    t3$ec <- '+1/2'
    t4$ec <- '-3/4'
    t5$ec <- '-2/3'
    rt1 <- rbind(t2, t3, t4, t5)

    n.s <- 100
    t2 <- img.main(n.i, n.s, d.dat=.az.ec2)
    t3 <- img.main(n.i, n.s, d.dat=.az.ec3)
    t4 <- img.main(n.i, n.s, d.dat=.az.ec4)
    t5 <- img.main(n.i, n.s, d.dat=.az.ec5)
    t2$ec <- '-1/2'
    t3$ec <- '+1/2'
    t4$ec <- '-3/4'
    t5$ec <- '-2/3'
    rt2 <- rbind(t2, t3, t4, t5)

    n.s <- 200
    t2 <- img.main(n.i, n.s, d.dat=.az.ec2)
    t3 <- img.main(n.i, n.s, d.dat=.az.ec3)
    t4 <- img.main(n.i, n.s, d.dat=.az.ec4)
    t5 <- img.main(n.i, n.s, d.dat=.az.ec5)
    t2$ec <- '-1/2'
    t3$ec <- '+1/2'
    t4$ec <- '-3/4'
    t5$ec <- '-2/3'
    rt3 <- rbind(t2, t3, t4, t5)

    list(rt1=rt1, rt2=rt2, rt3=rt3)
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
