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
getSIM <- function(recache = FALSE)
{
    rds <- 'dat/sim_rpt.rds'
    if(file.exists(rds) && !recache)
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
            src = unname(c(G='G', V='V', A='A', X='I')[cf[4]]),
            typ = unname(c(L='C', B='D')[cf[5]]),
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
picSIM <- function(pwr)
{
    library(ggplot2)
    graphics.off()

    ## reorder power source
    pwr <- within(pwr,
    {
        odr <- c(N=1, G=2, V=3, A=4, I=5)[src]
        src <- as.factor(src)
        src <- reorder(src, odr)
    })
    
    ## basic plot elements
    g <- ggplot() + xlab('N') + ylab('P') + xlim(100, 800) + ylim(0, 1)
    g <- g + theme(legend.position = "bottom", legend.box = "horizontal")
    gp <- geom_point
    gl <- geom_line
    fw <- facet_wrap

    ## Continuous & Dichotomous responses
    dt <- list(
        C = subset(pwr, typ == 'C' & src != 'N'),
        D = subset(pwr, typ == 'D' & src != 'N'))
    rt <- lapply(dt, function(d)
    {
        g <- g + gp(data = d, aes(ssz, pwr, shape = knl), size = 1.7)
        g <- g + gl(data = d, aes(ssz, pwr, linetype = vtx, group = paste(knl, vtx)))

        ## facet label
        lb <- label_bquote(Y[.(as.character(src))]^.(typ))
        g <- g + fw(~ src + typ, 2, 2, labeller = lb)
        g <- g + theme(strip.text.x = element_text(family = 'times'))

        ## lengend
        ## vertex similarity
        lg.V <- guides(
            linetype = guide_legend(
                title = expression(S[..]^V),
                label.position = "bottom",
                label.hjust = 0.5))
        ## U statistics
        lg.U <- guides(
            shape = guide_legend(
                title = 'U',
                label.position = 'bottom',
                label.hjust = 0.5))
            
        g <- g + lg.V + lg.U

        f <- paste('PWR_', c(C='CNT', D='DCT')[d[1, 'typ']], '.png', sep = '')
        ggsave(f, g, 'png', 'rpt/img', dpi = 400)
        g
    })

    invisible(rt)
}

main <- function()
{
    rpt <- getSIM(T)
    pwr <- powSIM(rpt)
    pic <- picSIM(pwr)
    invisible(pic)
}
