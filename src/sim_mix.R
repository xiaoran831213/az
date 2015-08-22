source('src/gno.R')
source('src/utl.R')
source('src/hwu.R')
source('src/hlp.R')
source('src/sim_img.R')
source('src/sim_gno.R')

## randomly pick encoded image data from a folder
mix.sim <- function(
    img, gno, n.s = .Machine$integer.max,
    ve.sd=.4, ve.fr=.05, vt.nm='tck', vt.ec=4,
    ge.sd=.4, ge.fr=.05,
    pe=c(G_=.0, IG=.5, I_=1), pw=c(G_=.0, IG=.5, I_=1),
    ne.rt=1.0)
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
    gt <- gt[,sbj]

    ## check genotype degeneration
    gt = GNO$clr.dgr(gt)
    if(length(gt) == 0)
    {
        cat('null genotype.')
        return(NA)
    }
    n.g <- nrow(gt)
    
    ## * -------- [vertex effect] -------- *
    ## all vertex encoding
    vc <- subset(img$enc, grepl(vt.nm, names(img$enc)))
    
    ## encoding level 0 is in fact the unencoded vertices, transpose
    ## the code to subject column major so (ve * vt) could work!
    vt <- t(vc[[1]])
    
    ## assign effect to some vertices
    ve = rnorm(n.v, 0, ve.sd) * rbinom(n.v, 1L, ve.fr)
    
    ## vertex contributed phenotype
    y1.ve <- apply(ve * vt, 'sbj', mean)
    ## * -------- [vertex effect] -------- *


    ## * -------- [genome effect] -------- *
    ## assign effect to some genome
    ge <- rnorm(n.g, 0, ge.sd) * rbinom(n.g, 1L, ge.fr)
    
    ## genome contributed response
    y1.ge <- apply(ge * gt, 'sbj', mean)
    ## * -------- [genome effect] -------- *

    ## sanity check: image and genome must share subjects
    stopifnot(colnames(gt) == colnames(vt))

    ## * -------- [joint effect(s)] -------- *
    y1 <- lapply(pe, function(p) y1.ve * p + y1.ge * (1 - p))
    y1.mu <- lapply(y1, mean)
    y1.sd <- lapply(y1, sd)
    
    ## noise effect(s)
    ne <- lapply(y1.sd, function(s) rnorm(n.s, 0, ne.rt * s))
    
    ## null response(s) not affected by any variables
    y0.mu <- mean(unlist(y1.mu))
    y0.sd <- mean(unlist(y1.sd))
    y0 <- list(NL=rnorm(n.s, y0.mu, y0.sd))

    ## concatinate all responses
    y <- c(y1, y0)
    
    ## * -------- U sta and P val --------*
    gt <- t(gt)                         # HWU use row major subjet
    vt <- vc[[vt.ec + 1L]]              # pick out vertex codes
    
    ## get genomic and vertex weight of various proportion
    wt.vt <- .hwu.GUS(vt)
    wt.gt <- .hwu.IBS(gt)
    wt <- list(wG_=wt.gt, wIG=wt.vt*wt.gt, wI_=wt.vt)
    
    p0 <- try(mapply(y0, ne, wt, FUN = function(y, n, w)
    {
        hwu.dg2(y = y + n, w=w)
    }, SIMPLIFY = F))

    if(inherits(p0, 'try-error'))
    {
        cat(gno$str, p0)
        return(NA)
    }
    
    c(.record(), unlist(p0))
}

mix.main <- function(n.itr = 5, n.sbj = 200)
{
    ## pick genotypes and images
    sim.rpt <- replicate(
        n.itr,
    {
        gno <- gno.pck(.az.wgs.bin, replace = F)
        img <- img.pck(.az.img, replace = F)
        cat(gno$ssn, img$ssn, '\n')
        mix.sim(img=img, gno=gno, n.s=n.sbj)
    })

    ## report
    names(sim.rpt) <- sprintf('s%03X', 1L:length(sim.rpt))
    HLP$mktab(sim.rpt)
}
