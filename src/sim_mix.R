source('src/gno.R')
source('src/utl.R')
source('src/hwu.R')
source('src/hlp.R')
source('src/sim_img.R')
source('src/sim_gno.R')

## randomly pick encoded image data from a folder
mix.sim <- function(
    img, gno, n.s = .Machine$integer.max,
    ve.sd=1, ve.fr=.05, vt.nm='tck', vt.ec=4,
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
    img <- img.sbj.pck(img, sbj)
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
    
    ## encoding level 0 is in fact the unencoded vertices, transpose
    ## the code to subject column major so (ve * vt) could work!
    vt <- t(vc[[1]])
    
    ## vertex contributed phenotype
    ve <- rnorm(n.v, 0, ve.sd) * rbinom(n.v, 1L, ve.fr)
    y1.ve <- apply(ve * vt, 'sbj', sum)
#    y1.ve <- .map.std(y1.ve)
    
    ## * -------- [genome effect] -------- *
    ge <- rnorm(n.g, 0, ge.sd) * rbinom(n.g, 1L, ge.fr)
    y1.ge <- apply(ge * gt, 'sbj', sum)
#    y1.ge <- .map.std(y1.ge)
    
    ## sanity check: image and genome must share subjects
    stopifnot(colnames(gt) == colnames(vt))
    #y1.vg <- 0.5 * y1.ge + 0.5 * y1.ve
    y1.vg <- y1.ve + y1.ge
    
    ## * -------- [joint effect(s)] -------- *
    y <- list(
        V_=y1.ve + rnorm(n.s, 0, 3.0 * sd(y1.ve)), # vertex
        G_=y1.ge + rnorm(n.s, 0, 3.0 * sd(y1.ve)), # genome
        VG=y1.vg + rnorm(n.s, 0, 2.0 * sd(y1.vg)), # mix
        NL=rnorm(n.s, 0, 1))                       # null effect

    ## * -------- U sta and P val --------*
    gt <- t(gt)                         # HWU use row major subjet
    vt <- vc[[vt.ec + 1L]]              # pick out vertex codes
    
    ## get genomic and vertex weight of various proportion
    wt.vt <- .hwu.GUS(vt)
    wt.gt <- .hwu.IBS(gt)
    wt.vg <- wt.vt * wt.gt

    w <- list(
        CV__=.wct(wt.vt), CG__=.wct(wt.gt), CVG_=.wct(wt.vg))
    w <- within(w, CVCG <- CV__* CG__)
    
    pv <- lapply(y, function(y) lapply(w, function(w) hwu.dg2(y, w)))
    
    if(inherits(pv, 'try-error'))
    {
        cat(gno$str, pv)
        return(NA)
    }
    
    c(.record(), unlist(pv))
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
    }, simplify = FALSE)
    
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

