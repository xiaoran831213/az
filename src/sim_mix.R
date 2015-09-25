#!/usr/bin/env Rscript
source('src/gno.R')
source('src/utl.R')
source('src/hwu.R')
source('src/hlp.R')
source('src/img.R')
source('src/sim_rpt.R')

## randomly pick encoded image data from a folder
sim.mix <- function(
    im, gn, n.s = 200,
    ve.sd=.06, ve.fr=.05, vft='tck', vt.ec=c(1,5), vt.gb=c(1, 3),
    ge.sd=.06, ge.fr=.05, lwr=c(G=0, X=.5, V=1),
    ns.sd = 4)
{
    gn <- if(is.character(gn)) readRDS(gn) else gn
    im <- if(is.character(im)) readRDS(im) else im
    
    ## number of vertices
    n.v <- length(im$vtx)
    
    ## guess  genomic NA
    gt <- GNO$imp(gn$gmx)               # genomic matrix

    ## pick out subjects
    sbj <- intersect(im$sbj, gn$sbj)
    n.s <- min(n.s, length(sbj))
    sbj <- sample(sbj, n.s)
    im <- sbj.img(im, sbj)
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
    vc <- subset(im$enc, grepl(vft, names(im$enc)))
    names(vc) <- paste('E', 0:(length(vc)-1), sep='')
    
    ## vertex blur (gaussian)
    vb <- aperm(im$gsb[vft, , ,], perm=c('sbj', 'vtx', 'sdv'))
    dimnames(vb)$sdv <- paste('B', 0:(dim(vb)[3]-1), sep='')
    
    ## blur level 0 is raw vertices, transpose the matrix to
    ## one subject per column so (ve * vt) could work!
    vt <- t(vb[,,'B0'])

    ## vertex and genome effect
    ve <- rnorm(n.v, 0, ve.sd) * rbinom(n.v, 1L, ve.fr)
    ge <- rnorm(n.g, 0, ge.sd) * rbinom(n.g, 1L, ge.fr)
    
    ## phenotype composistion, including noise
    ns <- function(v, sd) v + rnorm(length(v), sd = sd)
    y <- within(
        list(),
    {
        ## generate dependent variable
        V <- apply(ve * vt, 'sbj', sum)  # vertex
        G <- apply(ge * gt, 'sbj', sum)  # genome
        G <- .sc1(G) #G <- if(is.na(sd(G))) 0 else scale(G)
        V <- .sc1(V) #V <- scale(V)
        A <- V + G
        X <- V + G + V * G
        
        ## assign noise
        V <- ns(V, ns.sd)
        G <- ns(G, ns.sd)
        A <- ns(A, ns.sd)
        X <- ns(X, ns.sd)

        N <- rnorm(n.s, 0, 1)
    })
    #y <- lapply(y, function(e) e + rnorm(n.s, 0, ne.rt * sd(e)))

    ## avoid recording scalar dropped from vector
    rm(ve, ge)
    
    ## * -------- U sta and P val --------*
    ## the shared genomic weights should be computed only onece
    wg <- .hwu.IBS(t(gt), std = T)      # HWU use row major subject
    
    ## regional tests
    p.rgn <- unlist(lapply(vc[vt.ec], function(vt)
    {
        wv <- .hwu.GUS(vt, std = T)     # coded vertex
        pv <- lapply(lwr, function(r)
        {
            if(r==0)
                wt <- wg                # bypass NaN from log(wv)
            else if(r==1)
                wt <- wv                # bypass NaN from log(wg)
            else
            {
                wt <- wv * wg
                wt <- (wt-min(wt))/(max(wt)-min(wt))
            }
            lapply(y, hwu.dg2, w = .wct(wt))
        })
        pv
    }))
    
    ## vertex wise tests, pick out some levels of g-blur
    ## p.vwa <- apply(vb[,, vt.gb], c('vtx', 'sdv'), function(v)
    ## {
    ##     lwv <- log(.hwu.GUS(as.matrix(v)))
    ##     pv <- unlist(lapply(c(X=lwr, V=1), function(r)
    ##     {
    ##         if(r==1)
    ##             lwt <- lwv       # bypass NaN from log(wg)
    ##         else
    ##             lwt <- r * lwv + (1-r) * lwg
    ##         wt <- exp(lwt)
    ##         lapply(y, hwu.dg2, w = .wct(wt))
    ##     }))
    ## })
    ## names(dimnames(p.vwa))[1] <- 'mdl'
    
    ## p.vwa <- apply(p.vwa, c('mdl', 'sdv'), function(p)
    ## {
    ##     c(
    ##         NN=min(p),
    ##         FD=min(p.adjust(p, 'fdr')),
    ##         BF=min(p.adjust(p, 'bon')),
    ##         AT=min(1, p*64))
    ## })
    ## names(dimnames(p.vwa))[1] <- 'adj'
    
    ## p.vwa <- aperm(p.vwa, c('adj', 'sdv', 'mdl'))
    ## nm <- expand.grid(dimnames(p.vwa))
    ## p.vwa <- as.vector(p.vwa)
    ## names(p.vwa) <- do.call(paste, c(nm, sep='.'))

    c(.record(), p.rgn)#, p.vwa)
}

.az.wgs <- Sys.getenv('AZ_WGS')
.az.gno <- paste(.az.wgs, 'bin', sep='.')
main.mix <- function(img=.az.img, gno=.az.gno, n.i=5L, ...)
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
        do.call(sim.mix, a)
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
