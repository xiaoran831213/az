source('src/hlp.R')
cat.rpt <- function(src, ...)
{
    ## pick out images by file name
    dirs <- c(src, ...)
    
    fns <- unlist(lapply(dirs, function(d) file.path(d, dir(d, '*.rds'))))
    dat <- lapply(fns, readRDS)
    do.call(rbind, dat)
}

pwr <- function(rpt, t = 0.05, ret=2)
{
    rpt$n.g <- NULL
    hdr <- names(rpt)
    ptw <- '[VGAXW0-9]+'              # pattern of weight symbol
    pte <- '[VGAXN]+$'                # pattern of effect symbol
    pts <- paste('[.]', ptw, '[.]', pte, sep='') 
    cfg <- rpt[, !grepl(pts, hdr)]
    
    if(ret == 0)                        # patten of effect symbol
        pte <- '[N]+$'
    if(ret == 1)                   
        pte <- '[VGAX]+$'

    ## pattern of simulation symbol
    pts <- paste('[.]', ptw, '[.]', pte, sep='')
    pvl <- rpt[, grepl(pts, hdr)]

    rej <- function(p)
    {
        round(sum(p<t, na.rm=T) / sum(!is.na(p)), digits=3L)
    }
    
    n.i <- as.vector(tapply(rownames(rpt), cfg, length))
    p.v <- aggregate(pvl, by=cfg, rej)

    cbind(n.i=n.i, p.v)
}

.cb <- function(...) expand.grid(..., stringsAsFactors = F)
.mt <- function(pt, tx) regmatches(tx, regexec(pt, tx))
.pk <- function(l, f) sapply(l, `[[`, f)
pic <- function(pwr, xts='n.s', et='01234567', wt="VGX", yt="VGX")
{
    nm <- names(pwr)
    cfg <- .mt(sprintf('^E([%s]).([%s]).([%s])$', et, wt, yt),  nm)
    names(cfg) <- nm
    cfg <- lapply(Filter(length, cfg), function(m)
    {
        m <- c(as.list(m), list(pwr[,m[1]]))
        names(m) <- list('nm', 'et', 'wt', 'yt', 'yv')
        m
    })
    cfg <- Filter(function(f) f$et == '0' || f$wt != 'G', cfg)

    yts <- unique(.pk(cfg, 'yt'))
    yr <- c(0,1)

    graphics.off()
    with(.cb(xts, yts), mapply(function(xt, yt)
    {
        png(sprintf('%s_Y.%s.png', xt, yt), width=1050, height=1050, res=144)
        gt <- sprintf(
            'Test against\n %s effect',
            switch(EXPR=yt, A='Additve', X='Interaction', G='Genetic', V='Vertex'))
        xv <- pwr[,xt]
        xr <- range(xv)

        xl <- 'sample size' ## should not hard code this!
        yl <- 'power'
        par(cex.lab = 1.3, cex.axis = 1.3, mar = c(4,4,4,1), mgp = c(2.5, 1, 0), lwd = 3)
        plot(1, type = 'n', main = gt,  xlim = xr, ylim = yr, xlab = xl, ylab = yl)
        abline(h = pretty(yr), v = pretty(xr), col = "lightgray", lwd = 1)

        ## functions to pick color, line type and line width for a
        ## weight kernel and encoding.
        pty <- function(w, e)
        {
            col <- switch(EXPR = w, X=2, G=3, V=4)
            if(e == '0' || w == 'G')
                lty = 1
            else
                lty = 3
            lwd = 0.5 * as.integer(e) + 1
            list(col = col, lty = lty) 
        }

        plt <- lapply(Filter(function(f) f$yt == yt, cfg), function(f)
        {
            types <- pty(f$wt, f$et)
            lines(x=xv, y=f$yv, col = types$col, lwd = types$lwd, lty = types$lty)
            list(w=f$wt, e=f$et)
        })

        ## the lengend
        lgd <- lapply(plt, function(p)
        {
            
            txt <- sprintf(
                '%s kernel%s',
                switch(EXPR=p$w, X='joint', G='genetc', V='vertex'),
                switch(EXPR=p$e, '0'='', ', vertex code'))
            typ <- pty(p$w, p$e)
            c(list(txt=txt), typ)
        })
        
        legend(
            'topleft', ncol = 1,
            legend=.pk(lgd, 'txt'), col=.pk(lgd, 'col'), lty=.pk(lgd, 'lty'), pt.cex = 0.5)
        dev.off()
    }, xts, yts))

    invisible()
}

tab <- function()
{
    rr1 <- readRDS('dat/rr1.rds')

    rv <- with(rr1, data.frame(ssn=wsn, region=wnm, size=n.v, pvl=E4.V))
    rv <- na.omit(unique(rv))
    rv$bon <- p.adjust(rv$pvl, 'bon')
    rv$fdr <- p.adjust(rv$pvl, 'fdr')
    rv <- with(rv, rv[order(pvl),])
    write.csv(rv, '~/Dropbox/rv.csv', row.names = F)

    rg <- with(rr1, data.frame(ssn=gsn, region=gnm, size=n.g, pvl=E4.G))
    rg <- unique(rg)
    rg <- na.omit(rg)
    rg$bon <- p.adjust(rg$pvl, 'bon')
    rg$fdr <- p.adjust(rg$pvl, 'fdr')
    rg <- with(rg, rg[order(pvl),])
    write.csv(rg, '~/Dropbox/rg.csv', row.names = F)

}
