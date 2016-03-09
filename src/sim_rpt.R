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
    ## remove the column of number of genomic variants
    rpt$n.g <- NULL

    ## replace missing symbols with '#'
    names(rpt) <- sub('^[.]', '#.', names(rpt))
    names(rpt) <- sub('[.]$', '#.', names(rpt))
    hdr <- names(rpt)

    ## symbolic patterns
    pp <- '[NF]'                        # p-value correction [NPFB]
    pd <- '[BE][0-9]'                   # data type
    pw <- '[VGAX]'                      # u-kernel
    pe <- '[NVGAX][LB]'                 # effect

    ## pick out configurations
    cfg <- rpt[, !grepl(pe, hdr)]

    ## patten of effect symbol
    if(ret == 0)                        
        pe <- '[NN][NLB]'
    if(ret == 1)                   
        pe <- '[NVGAX][NLB]'

    ## pattern of simulation symbol
    pts <- paste(pp, pd, pw, pe, sep='[.]')
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
    cfg <- Filter(function(f) f$wt == 'G' || f$et > 0, cfg)
    yts <- unique(.pk(cfg, 'yt'))
    yr <- c(0,1)

    graphics.off()
    with(.cb(xts, yts), mapply(function(xt, yt)
    {
        png(sprintf('%s_%s.png', xt, yt), width=1050, height=950, res=144)
        xv <- pwr[,xt]
        xr <- range(xv)

        xl <- 'sample size' ## should not hard code this!
        yl <- 'power'
        par(cex.lab = 1.3, cex.axis = 1.3, mar = c(4,4,1,1), mgp = c(2.5, 1, 0), lwd = 3)
        plot(1, type = 'n', main = NA,  xlim = xr, ylim = yr, xlab = xl, ylab = yl)
        abline(h = pretty(yr), v = pretty(xr), col = "lightgray", lwd = 1)

        ## functions to pick color, line type and line width for a
        ## weight kernel and encoding.
        pty <- function(w, e)
        {
            col <- switch(EXPR = w, X=2, G=3, V=4)
            if(e == '0' && w != 'G')
                lty = 3
            else
                lty = 1
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
            txt <- switch(
                EXPR=p$w,
                X=expression('U'[J]),
                G=expression('U'[G]),
                V=expression('U'[V]))
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

tab <- function(rr1)
{
    rr1 <- na.omit(rr1)

    rv <- with(rr1, data.frame(sn=wsn, nm=wnm, sz=n.v, pvl=E4.V))
    rv <- aggregate(formula = pvl ~ sn + nm + sz, FUN = mean, data = rv)
    rownames(rv) <- rv$sn
    rv$sn <- NULL
    rv <- with(rv, rv[order(pvl),])
    ##write.csv(rv, '~/Dropbox/rv.csv', row.names = F)
    
    rg <- with(rr1, data.frame(sn=gsn, nm=gnm, sz=n.g, pvl=E4.G))
    rg <- aggregate(formula = pvl ~ sn + nm + sz, FUN = mean, data = rg)
    rownames(rg) <- rg$sn
    rg$sn <- NULL
    rg <- with(rg, rg[order(pvl),])
    ##write.csv(rg, '~/Dropbox/rg.csv', row.names = F)

    rx <- with(
        rr1,
        data.frame(
            row.names = paste(wsn, gsn, sep='.'),
            wnm, gnm, nv=n.v, ng=n.g, pg=E4.G, pv=E4.V, px=E4.X))
    rx <- with(rx, rx[order(px),])

    list(rv=rv, rg=rg, rx=rx)
}

pix <- function(rx, pch=0, np = 1000)
{
    ## top 20
    library(xtable)
    tp20 <- subset(head(r1$rx, 20), select = c(wnm, gnm, pv, pg, px))
    tp20 <- format(tp20, digits = 3)
    names(tp20) <- c('cortical surface', 'gene', '$P_V$', '$P_G$', '$P_X$')
    tp20 <- xtable(tp20, 'Top 20 combinations', 'tb:tp20')
    print(tp20, file = 'rpt/tp20.tex', include.rownames = F, sanitize.text.function = identity)
    
    r2 <- rx[seq(1, nrow(rx), l=np), ]
    lpg <- -log10(r2$pg)
    lpv <- -log10(r2$pv)
    lpx <- -log10(r2$px)

    ## library(tikzDevice)
    ## tikz('rpt/q1a.fig.tex', width = 5, height = 5)
    png('rpt/lgp_xgv.png', width=1050, height=750, res=144)
    yl <- '-LG(P)'
    yr <- c(0, max(lpg, lpv, lpx))
    par(pch = pch, cex.lab = 1.3, cex.axis = 1.3, mar = c(4,4,1,1), mgp = c(2.5, 1, 0))

    x <- seq_along(lpx)
    plot(x=x, y=lpx, col='red', ylab=yl, xlab=NA)
    abline(h = pretty(yr), col = "lightgray", lwd = 1)
    points(x=x, y=lpg, col='green')
    points(x=x, y=lpv, col='blue')

    lgd <- c(X=expression('U'[J]), G=expression('U'[G]), V=expression('U'[V]))
    legend('topright', 
        legend=lgd, col=c('red', 'green', 'blue'), pch=pch, pt.cex = 1.5)

    dev.off()
}

## r1 <- readRDS('bin/rr1.rds')
## pix(r1$rx)
getReport <- function()
{
    ## load simulation output
    d0 <- readRDS('dat/s4~8.rds')
    d0 <- sapply(names(d0)[-(1:10)], function(x)
    {
        cf <- regexec('^(.)[.](..)[.](.)[.](.)(.)$', x)
        cf <- regmatches(x, cf)[[1]][-1]
        rt <- data.frame(
            ssz = d0[, 4],              # sample size
            vtx = cf[2],
            knl = unname(c(G='G', V='V', X='J')[cf[3]]),
            src = cf[4],
            typ = cf[5],
            mtd = unname(c(E0='RGN', E4='RGN', B2='VWA')[cf[2]]),
            adj = cf[1],
            pvl = d0[, x],
            stringsAsFactors = FALSE)
        rt
    }, simplify = F, USE.NAMES = F)
    d0 <- do.call(rbind, d0)
    d0 <- within(d0,
    {
        vtx[knl == 'G' ] <- 'NA'
    })

    ## power calculation
    pw <- function(x, t=0.05) sum(x<t)/length(x)
    fw <- pvl ~ ssz + vtx + knl + src + typ + mtd + adj
    d1 <- aggregate(formula = fw, FUN = pw, data = d0)
    
    d1 <- subset(d1, adj %in% c('N', 'F'))
    d1 <- within(d1,
    {
        pwr <- pvl
        rm(pvl, adj)
    })
    d1
}

pic <- function(d0)
{
    library(ggplot2)
    graphics.off()
    
    ## basic plot elements
    p0 <- ggplot() + xlab('N') + ylab('P') + xlim(100, 800) + ylim(0, 1)
    gp <- geom_point
    gl <- geom_line

    ## continuouse & binary responses
    ld <- expand.grid(
        src = unlist(strsplit('GVAXN', '')),
        typ = c('L', 'B'),
        stringsAsFactors = F)
    ds <- c(
        G = 'Genetic effect only',
        V = 'Vertex effect only',
        A = 'Additive (G+V)',
        X = 'Interactive (G+V+GV)',
        N = 'Irrelevent effect' )
    dt <- c(
        L = 'continuous response',
        B = 'binary response')
    
    ld <- apply(ld, 1L, function(x)
    {
        s <- x[1]; t <- x[2]
        d <- subset(d0, src==s & typ==t, select = -c(src, typ, mtd))
        p <- p0
        p <- p + gp(data = d, aes(ssz, pwr, shape = knl), size = 2L)
        p <- p + gl(data = d, aes(ssz, pwr, linetype = vtx, group = paste(knl, vtx)))
        tt <- ggtitle(paste(ds[s], dt[t], sep = ', '))
        p <- p + tt + theme(plot.title = element_text(lineheight=.8, face="bold"))
        dev.new(); print(p)
    })
}
