source('src/gno.R')
source('src/utl.R')
source('src/hwu.R')
source('src/hlp.R')

## randomly pick encoded image data from a folder
gno.sim <- function(gno, n.s=200L, ge.sd=.5, ge.fr=.25, ne.rt=3.0)
{
    ## pick subjects
    gno <- gno.sbj.pck(gno, sample(gno$sbj, n.s))

    ## * -------- [genome effect] -------- *
    gt <- gno$gmx                       # genomic matrix
    gt = GNO$clr.dgr(gt)                # clean degeneration
    
    if(length(gt) == 0)                 # give up the trial
    {
        cat('Empty: ', gno$str)
        return(NA)
    }
    gt = GNO$imp(gt)                    # imput missing g-variant
    n.g <- nrow(gt)

    ## assign effect to some genome
    ge <- rnorm(n.g, 0, ge.sd) * rbinom(n.g, 1L, ge.fr)
    
    ## genome contributed response
    y1 <- apply(ge * gt, 'sbj', mean)
    y1.mu <- mean(y1)
    y1.sd <- sd(y1)
    
    ## noise effect
    ne <- rnorm(n.s, 0, ne.rt * y1.sd)

    ## a null response not affected by any variables
    y0 <- rnorm(n.s, y1.mu, y1.sd)
    
    ## * -------- U sta and P val --------*
    ## HWU requires subjes be of row major
    gt <- t(gt)
    ## wg1 <- .hwu.IBS1(gt)
    ## wg2 <- .hwu.IBS2(gt)
    wg3 <- .hwu.IBS3(gt)
    
    ## p1.0 <- try(hwu.dg2(y = y0 + ne, w = wg1))
    ## if(inherits(p1.0, 'try-error'))
    ##      p1.0 <- NA

    ## p1.1 <- try(hwu.dg2(y = y1 + ne, w = wg1))
    ## if(inherits(p1.1, 'try-error'))
    ##      p1.1 <- NA

    ## p2.0 <- try(hwu.dg2(y = y0 + ne, w = wg2))
    ## if(inherits(p2.0, 'try-error'))
    ##     p2.0 <- NA

    ## p2.1 <- try(hwu.dg2(y = y1 + ne, w = wg2))
    ## if(inherits(p2.1, 'try-error'))
    ##     p2.1 <- NA

    p3.0 <- try(hwu.dg2(y = y0 + ne, w = wg3))
    if(inherits(p3.0, 'try-error'))
        p3.0 <- NA

    p3.1 <- try(hwu.dg2(y = y1 + ne, w = wg3))
    if(inherits(p3.1, 'try-error'))
        p3.1 <- NA

    #c(.record(), p1.0=p1.0, p1.1=p1.1, p2.0=p2.0, p2.1=p2.1)
    c(.record(), p3.0=p3.0, p3.1=p3.1)
}

.az.wgs <- Sys.getenv('AZ_WGS')
.az.wgs.bin <- paste(.az.wgs, 'bin', sep='.')
gno.main <- function(n.itr = 5, n.sbj = 200, g.dat = NULL)
{
    if(is.null(g.dat))
    {
        cat('load', n.itr, 'genome data.\n')
        g.dat <- gno.pck(.az.wgs.bin, size = n.itr, vbs = T)
    }
    
    sim.rpt <- lapply(g.dat, function(g)
    {
        cat(gno.str(g), '\n')
        gno.sim(g, n.s=n.sbj)
    })

    ## report
    names(sim.rpt) <- sprintf('s%03X', 1L:length(sim.rpt))
    HLP$mktab(sim.rpt)
}

gno.pwr <- function(rpt, t = 0.05)
{
    n.itr <- nrow(rpt)
    p.hdr <- grepl('p[0-9].[01]$', colnames(rpt))
    p.val <- subset(rpt, select=p.hdr)
    lapply(p.val, function(p) sum(p < t, na.rm = T) / n.itr)
}
