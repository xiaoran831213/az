source('src/hlp.R')
cat.rpt <- function(src, ...)
{
    ## pick out images by file name
    dirs <- c(src, ...)
    
    fns <- unlist(lapply(dirs, function(d) file.path(d, dir(d, '*.rds'))))
    dat <- lapply(fns, readRDS)
    do.call(rbind, dat)
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

## Simulation report
getSIM <- function()
{
    rds <- 'dat/sim_rpt.rds'
    if(file.exists(rds))
        return(invisible(readRDS(rds)))

    ## list simulation output files
    f0 <- expand.grid('sim/t', 1:8, letters[1:8], stringsAsFactors = F)
    f0 <- do.call(paste, c(f0, sep=''))
    f0 <- unlist(sapply(f0, dir, 'rds$', full.name = T, USE.NAMES = F))

    ## load simulation output data
    d0 <- lapply(f0, readRDS)
    d0 <- do.call(rbind, d0) 

    ## rearrange the table
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
        vtx[knl == 'G' ] <- 'NA'        # a gnomic test needs no vertex
    })

    ## save and return
    saveRDS(d0, rds)
    invisible(d0)
}

## power calculation of simulation
powSIM <- function(d0)
{
    ## power calculation
    pw <- function(x, t = 0.05) sum(x < t) / length(x)
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

## picture of power from simulation reports
picSIM <- function(p0)
{
    library(ggplot2)
    graphics.off()
    
    ## basic plot elements
    g <- ggplot() + xlab('N') + ylab('P') + xlim(100, 800) + ylim(0, 1)
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
        d <- subset(p0, src==s & typ==t, select = -c(src, typ, mtd))
        p <- g + gp(data = d, aes(ssz, pwr, shape = knl), size = 2L)
        p <- p + gl(data = d, aes(ssz, pwr, linetype = vtx, group = paste(knl, vtx)))
        tt <- ggtitle(paste(ds[s], dt[t], sep = ', '))
        p <- p + tt + theme(plot.title = element_text(lineheight=.8, face="bold"))
        dev.new(); print(p)
    })
}
