source('src/gno.R')
source('src/utl.R')
source('src/hwu.R')
source('src/hlp.R')
source('src/img.R')
source('src/sim_gno.R')

## randomly pick encoded image data from a folder
mix.sim <- function(
    img, gno, n.s = .Machine$integer.max,
    ve.sd=1, ve.fr=.05, vt.nm='tck', vt.ec=c(1, 5), vt.gb=c(1, 3),
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
    img <- pck.sbj.img(img, sbj)
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
    names(vc) <- paste('E', 0:(length(vc)-1), sep='')
    
    ## vertex blur (gaussian)
    vb <- aperm(img$gsb[vt.nm, , ,], perm=c('sbj', 'vtx', 'sdv'))
    dimnames(vb)$sdv <- paste('B', 0:(dim(vb)[3]-1), sep='')
    
    ## blur level 0 is original vertices, transpose the matrix to
    ## subject column major so (ve * vt) could work!
    vt <- t(vb[,,'B0'])
    
    ## vertex contributed phenotype
    ve <- rnorm(n.v, 0, ve.sd) * rbinom(n.v, 1L, ve.fr)
    y1.ve <- apply(ve * vt, 'sbj', sum)
    
    ## * -------- [genome effect] -------- *
    ge <- rnorm(n.g, 0, ge.sd) * rbinom(n.g, 1L, ge.fr)
    y1.ge <- apply(ge * gt, 'sbj', sum)
    
    ## sanity check: image and genome must share subjects
    stopifnot(colnames(gt) == colnames(vt))
    #y1.vg <- 0.5 * y1.ge + 0.5 * y1.ve
    y1.vg <- y1.ve + y1.ge
    
    ## * -------- [joint effect(s)] -------- *
    y <- list(
        V__=y1.ve + rnorm(n.s, 0, 3.0 * sd(y1.ve)), # vertex
        G__=y1.ge + rnorm(n.s, 0, 3.0 * sd(y1.ve)), # genome
        V_G=y1.vg + rnorm(n.s, 0, 2.0 * sd(y1.vg)), # mix
        NUL=rnorm(n.s, 0, 1))                       # null effect

    ## avoid collecting scalars degenarated from vectors
    rm(y1.ve, y1.ge, y1.vg)

    ## * -------- U sta and P val --------*
    ## the shared genomic weights should be computed only onece
    wt.gt <- .hwu.IBS(t(gt))            # HWU use row major subject
    wt.gt.ct <- .wct(wt.gt)

    ## regional tests
    p.rgn <- lapply(vc[vt.ec], function(vt)
    {
        wt.vt <- .hwu.GUS(vt)               # coded vertex
        c(
            V_=sapply(y, hwu.dg2, w=.wct(wt.vt)),
            G_=sapply(y, hwu.dg2, w=wt.gt.ct),
            VG=sapply(y, hwu.dg2, w=.wct(wt.vt * wt.gt)))
    })
    p.rgn <- unlist(p.rgn)
    
    ## vertex wise tests, pick out some levels of g-blur
    p.vwa <- apply(vb[,, vt.gb], c('vtx', 'sdv'), function(v)
    {
        wv <- .hwu.GUS(as.matrix(v))
        c(
            V_=sapply(y, hwu.dg2, w=.wct(wv)),
            VG=sapply(y, hwu.dg2, w=.wct(wv * wt.gt)))
    })

    names(dimnames(p.vwa))[1] <- 'mdl'
    p.vwa <- apply(p.vwa, c('mdl', 'sdv'), function(p)
    {
        c(
            NN=min(p),
            FD=min(p.adjust(p, 'fdr')),
            BF=min(p.adjust(p, 'bon')))
    })
    names(dimnames(p.vwa))[1] <- 'adj'
    p.vwa <- aperm(p.vwa, c('sdv', 'mdl', 'adj'))
    nm <- expand.grid(dimnames(p.vwa))
    p.vwa <- as.vector(p.vwa)
    names(p.vwa) <- do.call(paste, c(nm, sep='.'))
    
    ## vertex wise test (vwa)
    ## if(inherits(p.rgn, 'try-error'))
    ## {
    ##     cat(gno$str, rgn=p.rgn)
    ##     return(NA)
    ## }
    
    c(.record(), p.rgn, p.vwa)
}

mix.main <- function(gno = .az.wgs.bin, img = .az.img.sm2, n.i = 5, ...)
{
    gno <- gno.pck(gno, size = n.i, replace = F)
    img <- pck.img(img, size = n.i, replace = F)
    rpt <- mapply(mix.sim, gno, img, MareArgs=..., SIMPLIFY = F)
    ##     n.itr,
    ## {
    ##     g <- gno.pck(gno, replace = F)
    ##     i <- pck.img(gno, replace = F)
    ##     cat(gno$ssn, img$vtx[1], '\n')
    ##     mix.sim(img=img, gno=gno, ...)
    ## }, simplify = FALSE)
    
    ## report
    names(sim.rpt) <- sprintf('s%03X', 1L:length(sim.rpt))
    HLP$mktab(sim.rpt)
}

mix.pwr <- function(rpt, t = 0.05, ret=3)
{
    n.itr <- nrow(rpt)
    if(ret == 0)
        rgx <- '[N0_]*[.][CGV_]*$'       # type 1 error
    else if(ret == 1)                   
        rgx <- '[GV_]*[.][CGV_]*$'       # power
    else
        rgx <- '[GV_NL0]*[.][CGV_]*$'     # both

    p.hdr <- grepl(rgx, colnames(rpt))
    p.val <- subset(rpt, select=p.hdr)
    
    lapply(p.val, function(p) sum(p < t, na.rm = T) / sum(!is.na(p)))
}

.mix.settings <- function()
{
    
}
