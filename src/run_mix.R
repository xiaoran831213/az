#!/usr/bin/env Rscript
source('src/gno.R')
source('src/utl.R')
source('src/hwu.R')
source('src/hlp.R')
source('src/img.R')

phe.mix <- function(pf = 'dat/sbj.csv')
{
    d <- read.csv(pf, as.is = T)

    # use baseline case and control
    d <- subset(
        d,
        subset=(VISCODE == 'bl' & PTETHCAT == "Not Hisp/Latino"
        & DX_bl %in% c('AD', 'CN') & PTRACCAT == "White"),
        select = c(
            PTID,
            VISCODE,                    # visit code
            DX_bl,                      # baseline Dx
            AGE,                        # baseline age
            PTGENDER,                   # sex
            PTEDUCAT,                   # edudation
            PTETHCAT,                   # ethnicity
            PTRACCAT,                   # race
            PTMARRY,                    # mariage
            APOE4))
    
    d <- with(
        d,
    {
        data.frame(
            sbj=PTID,
            dgn=ifelse(DX_bl=='CN', 0, 1),
            age=AGE,
            sex=ifelse(PTGENDER=='Male', 0, 1),
            edu=PTEDUCAT,
            ap4=APOE4, stringsAsFactors = F)
    })
    rownames(d) <- d$sbj
    d <- na.omit(d)
    d
}

## randomly pick encoded image data from a folder
run.mix <- function
(
    im, gn, pf,                         # image, genome, profile
    vft='tck',                          # vertex feature
    rsv='dgn',                          # response variable
    cfv=c('age', 'sex', 'edu', 'ap4'),  # confound variable
    vt.ec=c(4), lwr=c(G=0, X=.5, V=1))
{
    gn <- if(is.character(gn)) readRDS(gn) else gn
    im <- if(is.character(im)) readRDS(im) else im
    pf <- if(is.character(pf)) phe.mix(pf) else pf
    ## number of vertices
    n.v <- length(im$vtx)
    
    ## guess  genomic NA
    ## gt <- GNO$imp(gn$gmx)               # genomic matrix
    gt <- gn$gmx
    vt <- im$enc

    ## pick out subjects
    sbj <- Reduce(intersect, list(im$sbj, gn$sbj, pf$sbj))
    im <- sbj.img(im, sbj)
    gt <- gt[, sbj, drop = F]
    pf <- pf[sbj, ]
    n.s <- length(sbj)
    
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
    vc <- sprintf('%s[.][%s]', vft, paste(vt.ec, sep='', collapse = ''))
    vc <- subset(im$enc, grepl(vc, names(im$enc)))
    names(vc) <- sub(paste(vft, '[.]', sep = ''), 'E', names(vc))
    
    ## * -------- U sta and P val --------*
    ## the shared genomic weights should be computed only onece
    wg <- .hwu.IBS(t(gt), std = T)      # HWU use row major subject
    y <- as.matrix(pf[, rsv, drop = F])
    x <- as.matrix(pf[, cfv])
    
    ## regional tests
    p.rgn <- unlist(lapply(vc, function(vt)
    {
        wv <- .hwu.GUS(vt, std = T)     # coded vertex
        pv <- lapply(lwr, function(r)
        {
            if(r==0)
                wt <- wg                # bypass NaN from log(wv)
            else if(r==1)
                wt <- wv                # bypass NaN from log(wg)
            else
                wt <- .sc1(wv * wg)
            hwu.dg2(y, w = .wct(wt), x = x)
        })
        pv
    }))
    c(.record(), p.rgn,
      gsn=gn$ssn, gnm=.gen[gn$ssn, 'GEN'],
      wsn=im$wsn, wnm=im$wms)
}

.az.wgs <- Sys.getenv('AZ_WGS')
.az.gno <- paste(.az.wgs, 'bin', sep='.')
.az.img <- Sys.getenv('AZ_AENC')
.az.out <- Sys.getenv('AZ_REAL')
.gen <- .ls.seg()
main.mix <- function(img=.az.img, gno=.az.gno, n.i = NULL, ...)
{
    ## get sample lists, truncate to the shortest one
    gno <- pck(gno, size = n.i, replace=F, ret='f')
    img <- pck(img, size = n.i, replace=F, ret='f')
    phe <- 'dat/sbj.csv'

    ## get extra arguments
    dot <- list(...)
    
    if(length(dot) < 1L)
        args <- data.frame(im=img, gn=gno, pf=phe, stringsAsFactors = F)
    else
    {
        args <- expand.grid(
            ..., idx=1L:n.i, KEEP.OUT.ATTRS = F, stringsAsFactors = F)
        args <- args[do.call(order, args[, names(dot)]), ]

        args <- within(
            args,
        {
            im <- img[idx]
            gn <- gno[idx]
            rm(idx)
        })
    }
    
    ## reset row numbers so split(.) won't bring back unsorted table
    args <- sapply(1L:nrow(args), function(i) as.list(args[i,]), simplify = F)
    
    ## repeatative simulation
    simple_fn<- function(fn) sub('[.]rds$', '', basename(fn))
    rpt <- lapply(args, function(a)
    {
        print(simple_fn(unlist(a)))
        do.call(run.mix, a)
    })
    HLP$mktab(rpt)
}

main.cmb <- function(chr, wsn, img=.az.img, gno=.az.gno, n.i=NULL)
{
    gs <- rownames(subset(.gen, CHR %in% chr))
    gs <- sprintf('%s/%s.rds', gno, gs)
    gs <- gs[file.exists(gs)]
    if(length(n.i) >= 2L)
        gs <- gs[n.i]
    if(length(n.i) == 1L && n.i > 0L)
        gs <- sample(gs, n.i)

    im <- readRDS(sprintf('%s/%s.rds', img, wsn))
    im$wsn <- wsn
    pf <- phe.mix()

    rpt <- lapply(gs, function(gn)
    {
        gn <- readRDS(gn)
        cat(im$wsn, gn$ssn, '\n')
        run.mix(im, gn, pf, lwr=c(G=0, X=.5, V=1))
    })
    rpt <- HLP$mktab(rpt)
    rpt
}

cml.mix <- function()
{
    library(argparser)
    p <- arg_parser('AZ image genetics sumulation.')
    p <- add_argument(
        p, 'img', help = 'image source.')
    p <- add_argument(
        p, 'gno', help = 'gnome source.')
    p <- add_argument(
        p, '--dst', help = 'target to store reports.')
    p <- add_argument(
        p, '--chr', help = 'chromosome.',
        default='ch03')
    p <- add_argument(
        p, '--wms', help = 'white matter surface.',
        default='lh04')
    p <- add_argument(
        p, '--n.i', help = 'number of iterations.',
        default = 0L)
    p <- add_argument(
        p, '--vft', help = 'comma seperated vertex features',
        default='tck')
    p <- add_argument(
        p, '--rlw', help = 'ratio of logged weight kernels',
        default='0,0.5,1')

    argv <- commandArgs(trailingOnly = TRUE)
    if(length(argv) > 0)
    {
        opt <- parse_args(p, argv)
        print(opt)
        with(
            opt,
        {
            chr <- as.integer(sub('[chr]*', '', chr))
            rpt <- main.cmb(chr, wms, img, gno, n.i=n.i)
            saveRDS(rpt, dst)
        })
        cat('xt: success\n')
    }
}

cml.mix()
