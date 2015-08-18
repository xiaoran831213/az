source('src/gno.R')
source('src/utl.R')
source('src/hwu.R')
source('src/hlp.R')

## randomly pick encoded image data from a folder
gno.sim <- function(gno, n.s=200L, ge.sd=.5, ge.fr=.05, ne.rt=3.0)
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
    wg <- hwu.weight.IBS(gt)
    
    p0=try(hwu.dg2(y = y0 + ne, w = wg))
    if(inherits(p0, 'try-error'))
    {
        saveRDS(gno, 'gno.sim.err.rds')
        stop(p0)
    }
    
    p1=try(hwu.dg2(y = y1 + ne, w = wg))
    if(inherits(p1, 'try-error'))
    {
        saveRDS(gno, 'gno.sim.err.rds')
        stop(p1)
    }

    c(.record(), p0=p0, p1=p1)
}

gno.main <- function(n.itr = 5, n.sbj = 200)
{
    ## pick genotypes and images
    gno.dir <- seg.pck(vcf.dir = Sys.getenv('AZ_WGS'), size=n.itr, drop=F)

    ## run through simulations
    sim.rpt <- lapply(gno.dir, FUN = function(seg)
    {
        gno <- seg.get(seg)
        cat(gno$str, '\n')
        gno.sim(gno=gno, n.s=n.sbj)
    })

    ## report
    names(sim.rpt) <- sprintf('s%03X', 1L:length(sim.rpt))
    HLP$mktab(sim.rpt)
}

.power <- function(rpt, t = 0.05)
{
    n.itr <- nrow(rpt)
    p.hdr <- grepl('p[01]$', colnames(rpt))
    p.val <- subset(rpt, select=p.hdr)
    lapply(p.val, function(p) sum(p < t, na.rm = T) / n.itr)
}
