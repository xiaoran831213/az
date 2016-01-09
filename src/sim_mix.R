#!/usr/bin/env Rscript
library(boot)
library(parallel)
source('src/gno.R')
source('src/utl.R')
source('src/hwu.R')
source('src/hlp.R')
source('src/img.R')
source('src/sim_rpt.R')
source('src/dsg.R')

## vertex wise tests with some levels of gaussian-blur
## y  : list of responses
## vb : vertex data with a number of blur levels
## wg : genomic weight matrix
## lwr: weight ratio of the joint kernel 0=G, 0.5=X, 1=V
sim.vwa <- function(y, vb, wg, lwr, timer = F)
{
    ## the iteration over all vertices
    time0 <- proc.time()

    cl <- makeCluster(8)
    p.vwa <- parApply(cl,  vb, c('vtx', 'sdv'), function(v)
    {
        source('src/hwu.R')
        wv <- .hwu.GUS(as.matrix(v))
        pv <- lapply(lwr, function(r)
        {
            if(r==0)
                wt <- wg                # bypass NaN from log(wv)
            else if(r==1)
                wt <- wv                # bypass NaN from log(wg)
            else
            {
                wt <- wv * wg
                wt <- (wt - min(wt))/(max(wt) - min(wt))
            }
            lapply(y, hwu.dg2, w = .wct(wt))
        })
        unlist(pv)
    })
    stopCluster(cl)
    
    names(dimnames(p.vwa))[1] <- 'mdl'
    p.vwa <- apply(p.vwa, c('mdl', 'sdv'), function(p)
    {
        c(P=min(p),                   # raw p-value
          F=min(p.adjust(p, 'fdr')),  # false discovery rate
          B=min(p.adjust(p, 'bon')))  # bonferroni correction
    })
    names(dimnames(p.vwa))[1] <- 'adj'
    
    p.vwa <- aperm(p.vwa, c('adj', 'sdv', 'mdl'))
    nm <- expand.grid(dimnames(p.vwa))
    p.vwa <- as.vector(p.vwa)
    names(p.vwa) <- do.call(paste, c(nm, sep='.'))

    if(timer)
        print(proc.time() - time0)
    p.vwa
}
## randomly pick encoded image data from a folder
sim.mix <- function(im, gn, n.s = 300,
                    ve.sd=1, ve.fr=.15, vft='tck', vt.ec=c(1, 5), 
                    ge.sd=1, ge.fr=.15, lwr=c(X=.5, V=1),
                    vt.gb=NULL, ns.sd = 5)
{
    gn <- dosage(gn)
    im <- vimage(im)
    
    ## guess genomic NA
    gn <- impute(gn)

    ## pick out common subjects, then take the sample
    cs <- sample(intersect(sbj(im), sbj(gn)), min(nsb(im), nsb(gn), n.s))
    im <- sbj(im, cs)
    gn <- sbj(gn, cs)
    
    ## check genotype degeneration
    gn = rmDgr(gn)
    if(length(gn) == 0)
        return('NG')

    ## number of vertices & g-variants
    n.v <- nvt(im); n.g <- ngv(gn)
    
    ## all vertex code of one feature
    vc <- subset(im$enc, grepl(vft, names(im$enc)))
    names(vc) <- paste('E', 0:(length(vc)-1), sep='')
    
    ## * -------- [vertex & genomic effect] -------- *
    ## get raw surface vertices, scale the value to [0, 1] 
    vx <- .sc1(im$sfs[vft, , ])         # vertex value
    gt <- gn$gmx                        # genome value
    
    ve <- rnorm(n.v, 0, ve.sd) * rbinom(n.v, 1L, ve.fr)
    ge <- rnorm(n.g, 0, ge.sd) * rbinom(n.g, 1L, ge.fr)
    
    ## composisit effect, and assign noise
    .e <- within(list(),
    {
        ## generate dependent variable
        V <- apply(ve * vx, 'sbj', sum)  # vertex
        G <- apply(ge * gt, 'sbj', sum)  # genome
        G <- .std(G)
        V <- .std(V)
        A <- V + G
        X <- V + G + V * G

        ## irrelevant effect for type I error check
        ## N <- rnorm(length(cs))
    })
    .e <- lapply(.e, function(x) x + rnorm(length(x), 0, ns.sd))
    
    ## avoid recording scalar dropped from vector
    rm(ve, ge, cs)

    ## generate linear response: 
    y.lnr <- .e
    names(y.lnr) <- paste(names(.e), "L", sep = '')
    
    ## generate binary response: 
    y.bin <- lapply(.e, function(et)
    {
        p <- inv.logit(et)
        b <- rbinom(length(et), 1, p)
        b
    })
    names(y.bin) <- paste(names(.e), "B", sep = '')

    ## y.bin.rsd <- lapply(y.bin, function(b)
    ## {
    ##     reg <- glm(b ~ 1, family = binomial)
    ##     residuals(reg)
    ## })
    ## names(y.bin.rsd) <- paste(names(.e), "R", sep = '')

    ## combin all responses
    y <- c(y.lnr, y.bin) #, y.bin.rsd)
    
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
                wt <- .sc1(wv * wg)
            lapply(y, hwu.dg2, w = .wct(wt))
        })
        pv
    }))
    names(p.rgn) <- paste('', names(p.rgn), sep = '.')
    
    ## vertex wise tests, pick out some levels of g-blur
    p.vwa <- if(!is.null(vt.gb))
    {
        vb <- im$gsb[vft, , , ,drop = T] # only one feature at a time
        vb <- vb[, vt.gb, , drop = F]    # may be one or more g-blur
        vb <- aperm(vb, perm=c('sbj', 'vtx', 'sdv'))
        dimnames(vb)$sdv <- paste('B', vt.gb - 1, sep='')
        p.vwa <- sim.vwa(y, vb, wg, lwr)
    }
    else NULL
    
    c(.record(), p.rgn, p.vwa)
}

main.mix <- function(img='dat/img', gno='dat/gno', n.i=5L, ...)
{
    ## get sample lists, truncate to the shortest one
    gno <- pck(gno, size = n.i, replace=F, ret='f')
    img <- pck(img, size = n.i, replace=F, ret='f')
    n.i <- min(length(gno), length(img))

    ## get extra arguments, also remove NULL arguments so they won't
    ## cause havoc in expand.grid
    dot <- list(...)
    dot <- dot[!sapply(dot, is.null)]
    
    if(length(dot) < 1L)
        args <- data.frame(im=img, gn=gno, stringsAsFactors = F)
    else
    {
        ## idx = 1L:n.i expands the repetitions.
        args <- list(idx = 1L:n.i, KEEP.OUT.ATTRS = F, stringsAsFactors = F)
        args <- do.call(expand.grid, c(dot, args))
        args <- args[do.call(order, args[, names(dot), drop = F]), ]
        args <- within(
            args,
        {
            im <- img[idx]
            gn <- gno[idx]
            rm(idx)
        })
    }
    
    ## turn data.frame to a list of lists
    args <- tab2lol(args)
    
    ## repeatative simulation
    simple_fn<- function(fn) sub('[.]rds$', '', basename(fn))
    rpt <- lapply(args, function(a)
    {
        print(simple_fn(unlist(a)))
        do.call(sim.mix, a)
    })
    lol2tab(rpt)
}

## command line parser
.cmd <- function(argv = NULL)
{
    if(is.null(argv))
       argv <- commandArgs(trailingOnly = TRUE)

    ## code develplent mode
    if(length(argv) < 1L)
        return(invisible(NULL))

    ## actual run
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
        p, '--n.s', help = '"," seperated sample sizes.',
        default='50')
    p <- add_argument(
        p, '--vft', help = '"," seperated vertex features',
        default='tck,slc')
    p <- add_argument(
        p, '--vwa', help = 'vertex-wise analysis (VWA)',
        default=FALSE)
    p <- add_argument(
        p, '--vt.gb', help = '"," delimeted vertex gaussian blur level. (Def=1, no blur)',
        default=1)
    p <- add_argument(
        p, '--lwr', help = 'logged weight ratio of U kernels',
        default='0,0.5,1')

    opt <- parse_args(p, argv)
    opt <- within(opt, {rm(help); rm(opts)})
    
    print(opt)
        
    with(
        opt,
    {
        n.s <- as.integer(unlist(strsplit(n.s, ',')))
        lwr <- as.double(unlist(strsplit(lwr, ',')))
        vft <- unlist(strsplit(vft, ','))
        
        if(!vwa)
            vt.gb <- NULL
        if(!grepl('.rds$', dst, T))
            dst <- paste(dst, 'rds', sep='.')

        rpt <- main.mix(img, gno, n.i=itr, n.s=n.s, vft=vft, vt.gb=vt.gb)
        saveRDS(rpt, dst)
    })
    cat('xt: success\n')
}

.cmd()
