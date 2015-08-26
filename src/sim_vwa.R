source('src/gno.R')
source('src/utl.R')
source('src/hwu.R')
source('src/hlp.R')
source('src/sim_img.R')
source('src/sim_gno.R')
library(igraph)

.vwa.proc <- function(img)
{
    ## vertex coordinate for each subject
    xyz <- aperm(img$sfs[c('x', 'y', 'z'), ,], c(2, 1, 3))
    cmx <- img$cmx
    sbj <- img$sbj
    
    ## vertex euclidean distance for each subject
    vds <- apply(xyz, 3L, function(p)
    {
        as.matrix(dist(p))
    })
    dim(vds) <- c(dim(cmx), length(sbj))
    dimnames(vds) <- c(dimnames(cmx), list(sbj=sbj))

    ## vertex graph, weight by geodesic distance
    vgp <- apply(vds, 3L, function(d)
    {
        adj <- d * cmx
        graph_from_adjacency_matrix(adj, mode='upper', weighted=T)
    })
    
    img$vds <- vds
    img$vgp <- vgp
    img
}

## vertex wise analyais
.vwa.wt <- function(wt.gt, vt, ys, q = 32)
{
    ps <- lapply(ys, function(y) apply(vt, 2L, function(v)
    {
        wt.vt <- .hwu.GUS(as.matrix(v))
        hwu.dg2(y, wt.vt * wt.gt)
    }))

    lapply(ps, min)
}

vwa.img.sim <- function(img, n.s = 200L, vt.ft = 'tck')
{
    ## pick subjects
    img <- img.sbj.pck(img, sample(img$sbj, n.s))
    
    N <- length(img$sbj)
    M <- length(img$vtx)
    enc <- img$enc

    ## for now we only use 1 feature, also rescaled it to [0, 1]
    vt <- subset(enc, grepl(vt.ft, names(enc)))[[1L]]
    
    ## assign effect to each vertex
    ve.mu <- 0.0
    ve.sd <- 0.5
    ve.fr <- 0.05
    ve = rnorm(n = M, ve.mu, ve.sd) * rbinom(n = M, 1L, ve.fr)
    
    ## vertex contributed phenotype
    z1 <- apply(ve * vt, 'sbj', mean)   # vertex effect * vertex value
    z1.mu <- mean(z1)
    z1.sd <- sd(z1)
    
    ## noise effect
    ne.sr <- 3.0                        # noise to vertex sd ratio
    ne <- rnorm(n = N, mean = 0, sd = ne.sr * z1.sd)

    ## another phenotype is not affected by vertices
    z0 <- rnorm(n = N, mean = z1.mu, sd = z1.sd)

    ## Derive U statistics, get P values of all encoding levels
    ys <- list(y1=z1, y0=z0)
    
    pval <- lapply(ys, function(y)
    {
        ps <- apply(vt, 2L, function(v)
        {
            w <- .hwu.GUS(as.matrix(v))
            hwu.dg2(y, w)
        })
        min(ps)
    })
    
    c(.record(), unlist(pval))
}

## randomly pick encoded image data from a folder
vwa.mix.sim <- function(
    img, gno, n.s = .Machine$integer.max,
    ve.sd=1, ve.fr=.05, vt.nm='tck', vt.ec=4,
    ge.sd=1, ge.fr=.05, ne.rt=.5)
{
    ## number of vertices and g-variants
    n.v <- length(img$vtx)

    ## guess  genomic NA
    gt <- GNO$imp(gno$gmx)              # genomic matrix
    
    ## pick out subjects
    sbj <- intersect(img$sbj, gno$sbj)
    n.s <- min(n.s, length(sbj))
    sbj <- sample(sbj, n.s)
    img <- img.sbj.pck(img, sbj)
    gt <- gt[,sbj, drop = F]

    ## check genotype degeneration
    gt = GNO$clr.dgr(gt)
    if(length(gt) == 0)
    {
        cat('null genome\n')
        return(NA)
    }
    n.g <- nrow(gt)
    
    ## * -------- [vertex effect] -------- *
    ## all vertex encoding
    vc <- subset(img$enc, grepl(vt.nm, names(img$enc)))
    
    ## encoding level 0 is in fact the unencoded vertices, transpose
    ## the code to subject column major so (ve * vt) could work!
    vt <- t(vc[[1]])
    
    ## vertex contributed phenotype
    ve <- rnorm(n.v, 0, ve.sd) * rbinom(n.v, 1L, ve.fr)
    y1.ve <- apply(ve * vt, 'sbj', mean)
    
    ## * -------- [genome effect] -------- *
    ge <- rnorm(n.g, 0, ge.sd) * rbinom(n.g, 1L, ge.fr)
    y1.ge <- apply(ge * gt, 'sbj', mean)
    ## if(any(is.na(y1.ge)))
    ##     return(NA)

    ## sanity check: image and genome must share subjects
    stopifnot(colnames(gt) == colnames(vt))
    y1.ge <- y1.ge + rnorm(n.s, 0, 3.1 * sd(y1.ge))
    y1.ve <- y1.ve + rnorm(n.s, 0, 3.1 * sd(y1.ve))

    ## * -------- [joint effect(s)] -------- *
    y1 <- list(V_=y1.ve, G_=y1.ge)

    sd.G_ <- sd(y1$G_)
    sd.V_ <- sd(y1$V_)
    #VG.sd <- sd(y1$VG)
    y1.mu <- lapply(y1, mean)
    y1.sd <- lapply(y1, sd)
    
    ## null response(s) not affected by any variables
    ## y0 <- rnorm(n.s, mean(unlist(y1.mu)), mean(unlist(y1.sd)))
    y0 <- rnorm(n.s, 0, 1)
    
    ## concatinate all responses
    y <- c(y1, list(NL=y0))
    y.sd <- c(y1.sd, list(NL=sd(y0)))
    
    ## add noise effect(s)
    ## y <- mapply(function(b, s)
    ## {
    ##     b + rnorm(n.s, 0, s * ne.rt)
    ## }, y, y.sd, SIMPLIFY = F)

    ## * -------- U sta and P val --------*
    gt <- t(gt)                         # HWU use row major subjet
    vt <- vc[[vt.ec + 1L]]              # pick out vertex codes
    
    ## get genomic and vertex weight of various proportion
    wt.gt <- .hwu.IBS(gt)
    wt.vt <- list(
        V_=wt.vt,
        G_=wt.gt,
        VG=wt.vt*wt.gt)

    pv <- lapply(y, function(y) lapply(w, function(w) hwu.dg2(y, w)))
    
    if(inherits(pv, 'try-error'))
    {
        cat(gno$str, pv)
        return(NA)
    }
    
    c(.record(), unlist(pv))
}

vwa.img.main <- function(n.itr = 5, n.sbj = 200)
{
    ## pick genotypes and images
    sim.rpt <- replicate(
        n.itr,
    {
        img <- img.pck(.az.img, replace = F)
        cat(img$ssn, '\n')
        vwa.img.sim(img=img, n.s=n.sbj)
    }, simplify = FALSE)
    
    ## report
    names(sim.rpt) <- sprintf('s%03X', 1L:length(sim.rpt))
    HLP$mktab(sim.rpt)
}
