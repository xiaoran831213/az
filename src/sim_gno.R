source('src/gno.R')
source('src/utl.R')
source('src/hwu.R')
source('src/hlp.R')

## randomly pick encoded image data from a folder
gno.sim <- function(gno, n.s=200L, ge.sd=.5, ge.fr=.25, ns.rt=2.1)
{
    ## * -------- [genome effect] -------- *
    gt <- gno$gmx                       # genomic matrix

    gt = GNO$imp(gt)                    # imput missing g-variant
    n.g <- nrow(gt)
    
    ## pick subjects
    n.s <- min(n.s, ncol(gt))
    gt <- gt[, sample.int(ncol(gt), n.s), drop = F]

    gt = GNO$clr.dgr(gt)                # clean degenerated g-variant
    if(length(gt) == 0)                 
    {
        cat('Empty segment')            # give up the trial if all
        return(NA)                      # g-variants are invalid
    }

    ## assign effect to some genome
    ge <- rnorm(n.g, 0, ge.sd) * rbinom(n.g, 1L, ge.fr)
    
    ## genome contributed response
    y1 <- apply(ge * gt, 'sbj', mean)
    y1.mu <- mean(y1)
    y1.sd <- sd(y1)
    rm(ge)                     # prevent recoding of univariable GE
    
    ## noise effect
    ne <- rnorm(n.s, 0, ns.rt * y1.sd)

    ## a null response not affected by any variables
    y0 <- rnorm(n.s, y1.mu, y1.sd)
    
    ## * -------- U sta and P val --------*
    ## HWU requires subjes be of row major
    gt <- t(gt)
    wg <- .hwu.IBS(gt, runif(n = ncol(gt)))
    
    p.0 <- hwu.dg2(y = y0 + ne, w = wg)
    p.1 <- hwu.dg2(y = y1 + ne, w = wg)
    
    c(.record())
}

.az.wgs <- Sys.getenv('AZ_WGS')
.az.gno <- paste(.az.wgs, 'bin', sep='.')
main.gno <- function(n.itr = 5, n.sbj = 200, g.dat = NULL)
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

pwr.gno <- function(rpt, t = 0.05)
{
    n.itr <- nrow(rpt)
    p.hdr <- grepl('p[0-9]*.*[01]$', colnames(rpt))
    p.val <- subset(rpt, select=p.hdr)
    lapply(p.val, function(p) sum(p < t, na.rm = T) / n.itr)
}
