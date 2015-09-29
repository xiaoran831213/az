#!/usr/bin/env Rscript
source('src/gno.R')
source('src/utl.R')
source('src/hwu.R')
source('src/hlp.R')
source('src/img.R')
source('src/sim_rpt.R')

phe.mix <- function()
{
    d <- read.csv('dat/sbj.csv', as.is = T)

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
            ap4=APOE4)
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
    cfv=c('age', 'sex', 'ap4'),         # confound variable
    vt.ec=c(4), lwr=c(G=0, X=.5, V=1))
{
    gn <- if(is.character(gn)) readRDS(gn) else gn
    im <- if(is.character(im)) readRDS(im) else im
    
    ## number of vertices
    n.v <- length(im$vtx)
    
    ## guess  genomic NA
    gt <- GNO$imp(gn$gmx)               # genomic matrix

    ## pick out subjects
    sbj <- intersect(im$sbj, gn$sbj, pf$sbj)
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
    vc <- subset(im$enc, grepl(vft, names(im$enc)))
    names(vc) <- paste('E', 0:(length(vc)-1), sep='')
    
    ## * -------- U sta and P val --------*
    ## the shared genomic weights should be computed only onece
    wg <- .hwu.IBS(t(gt), std = T)      # HWU use row major subject
    y <- pf[, rsv]
    x <- pf[, cfv]
    ## regional tests
    ##:ess-bp-start::browser@nil:##
browser(expr=is.null(.ESSBP.[["@3@"]]))##:ess-bp-end:##
    
    p.rgn <- unlist(lapply(vc[vt.ec], function(vt)
    {
        wv <- .hwu.GUS(vt, std = T)     # coded vertex
        pv <- lapply(lwr, function(r)
        {
            if(r==0)
                wt <- wg          # bypass NaN from log(wv)
            else if(r==1)
                wt <- wv          # bypass NaN from log(wg)
            else
                wt <- .sc1(wv * wg)
            hwu.dg2(y, w = .wct(wt), x = x)
        })
        pv
    }))
    
    c(.record(), p.rgn)
}

.az.wgs <- Sys.getenv('AZ_WGS')
.az.gno <- paste(.az.wgs, 'bin', sep='.')
.az.img <- Sys.getenv('AZ_AENC')
main.mix <- function(img=.az.img, gno=.az.gno, ...)
{
    ## get sample lists, truncate to the shortest one
    gno <- pck(gno, size = n.i, replace=F, ret='f')
    img <- pck(img, size = n.i, replace=F, ret='f')
    n.i <- min(length(gno), length(img))

    ## get extra arguments
    dot <- list(...)
    
    if(length(dot) < 1L)
        args <- data.frame(im=img, gn=gno, stringsAsFactors = F)
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

cml.mix <- function()
{
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
        default=1)
    p <- add_argument(
        p, '--n.s', help = 'comma seperated sample sizes.',
        default='50')
    p <- add_argument(
        p, '--vft', help = 'comma seperated vertex features',
        default='tck,slc')
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
            n.s <- as.integer(unlist(strsplit(n.s, ',')))
            rlw <- as.double(unlist(strsplit(rlw, ',')))
            vft <- unlist(strsplit(vft, ','))
            dst <- paste(dst, 'rds', sep='.')

            ## ve.sd <- c(.10, .15, .20, .25, .30)
            ## ns.sd <- seq(.10, 1.0, by=.200)
            ve.fr <- c(0.01, 0.02, 0.03, 0.04, 0.05)
            rpt <- main.mix(img, gno, n.i=itr, n.s=n.s, vft=vft)
            saveRDS(rpt, dst)
        })
        cat('xt: success\n')
    }
}

cml.mix()
