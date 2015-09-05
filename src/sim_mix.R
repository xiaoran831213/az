#!/usr/bin/env Rscript
source('src/gno.R')
source('src/utl.R')
source('src/hwu.R')
source('src/hlp.R')
source('src/img.R')

## randomly pick encoded image data from a folder
sim.mix <- function(
    img, gno, n.s = 200,
    ve.sd=1, ve.fr=.05, vft='tck', vt.ec=c(1, 5), vt.gb=c(1, 3),
    ge.sd=1, ge.fr=.05, ne.rt=3)
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
    vc <- subset(img$enc, grepl(vft, names(img$enc)))
    names(vc) <- paste('E', 0:(length(vc)-1), sep='')
    
    ## vertex blur (gaussian)
    vb <- aperm(img$gsb[vft, , ,], perm=c('sbj', 'vtx', 'sdv'))
    dimnames(vb)$sdv <- paste('B', 0:(dim(vb)[3]-1), sep='')
    
    ## blur level 0 is raw vertices, transpose the matrix to
    ## one subject per column so (ve * vt) could work!
    vt <- t(vb[,,'B0'])
    
    ## vertex contributed phenotype
    ve <- rnorm(n.v, 0, ve.sd) * rbinom(n.v, 1L, ve.fr)
    y1.ve <- apply(ve * vt, 'sbj', sum)
    
    ## * -------- [genome effect] -------- *
    ge <- rnorm(n.g, 0, ge.sd) * rbinom(n.g, 1L, ge.fr)
    y1.ge <- apply(ge * gt, 'sbj', sum)
    
    ## sanity check: image and genome must share subjects
    stopifnot(colnames(gt) == colnames(vt))
    y1.vg <- y1.ve + y1.ge
    y1.mx <- y1.vg + y1.ve * y1.ge
    
    ## * -------- [joint effect(s)] -------- *
    y <- list(
        V_=y1.ve + rnorm(n.s, 0, ne.rt * sd(y1.ve)), # vertex
        G_=y1.ge + rnorm(n.s, 0, ne.rt * sd(y1.ve)), # genome
        VG=y1.vg + rnorm(n.s, 0, 2 * sd(y1.vg)),     # additive
        MX=y1.mx + rnorm(n.s, 0, 2 * sd(y1.mx)),     # mix
        NL=rnorm(n.s, 0, 1))                         # null effect

    ## avoid collecting scalars degenarated from vectors
    rm(ve, ge, y1.ve, y1.ge, y1.vg)

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
    
    c(.record(), p.rgn)#, p.vwa)
}

.az.wgs <- Sys.getenv('AZ_WGS')
.az.gno <- paste(.az.wgs, 'bin', sep='.')
main.mix <- function(img=.az.img, gno=.az.gno, n.i=5, ...)
{
    ## get extra arguments
    ## TODO: whty empty ... casuse error?
    arg <- expand.grid(...) 

    ## get sample lists, truncate to the shortest one
    gno <- pck.gno(gno, size = n.i, replace=F, ret='f')
    img <- pck.img(img, size = n.i, replace=F, ret='f')
    n.i <- min(length(gno), length(img))

    ## repeatative simulation
    rpt <- mapply(function(f.i, f.g)
    {
        d <- list(img=readRDS(f.i), gno=readRDS(f.g))
        cat(f.i, f.g, '\n')
        rpt <- apply(arg, 1L, function(a)
        {
            do.call(sim.mix, c(d, a))
        })
        HLP$mktab(rpt)
    }, img[1:n.i], gno[1:n.i], SIMPLIFY = F, USE.NAMES = F)
    do.call(rbind, rpt)
}

pwr.mix <- function(rpt, t = 0.05, ret=2)
{
    rpt$n.g <- NULL                 # don't group over g-variant count
    hdr <- names(rpt)
    if(ret == 0)
        pvl <- rpt[, grepl('[.][VG_]+[.]NL$', hdr)] # type 1 error
    else if(ret == 1)
        pvl <- rpt[, grepl('[.][VG_]+[.][GV_MX]+$', hdr)] # power
    else
        pvl <- rpt[, grepl('[.][VG_]+[.]', hdr)] # both

    cfg <- rpt[, !grepl('[.][VG_]+[.]', hdr)]
    
    rej <- function(p) sum(p<t, na.rm=T)/sum(!is.na(p))
    n.i <- as.vector(tapply(rownames(rpt), cfg, length))
    p.v <- aggregate(pvl, by=cfg, rej)

    cbind(n.i=n.i, p.v)
}

library(argparser)
p <- arg_parser('AZ image genetics sumulation.')
p <- add_argument(
    p, 'img', help = 'image source.')
p <- add_argument(
    p, 'gno', help = 'gnome source.')
p <- add_argument(
    p, '--dst', help = 'target to store reports.',
    default = 'sim_mix.rds')
p <- add_argument(
    p, '--itr', help = 'iteration to run.',
    default=5)
p <- add_argument(
    p, '--n.s', help = 'comma seperated sample sizes.',
    default='100,200,400')
p <- add_argument(
    p, '--vft', help = 'comma seperated vertex features',
    default='tck,slc')

argv <- commandArgs(trailingOnly = TRUE)
if(length(argv) > 0)
{
    opt <- parse_args(p, argv)
    attach(opt)
    n.s <- as.integer(unlist(strsplit(n.s, ',')))
    vft <- unlist(strsplit(vft, ','))
    print(n.s)
    print(vft)
    rpt <- main.mix(img, gno, n.i=itr, n.s=n.s, vft=vft)
    saveRDS(rpt, dst)
    detach(opt)
}
