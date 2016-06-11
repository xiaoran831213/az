source('src/hlp.R')
cat.rpt <- function(src, ...)
{
    ## pick out images by file name
    dirs <- c(src, ...)
    
    fns <- unlist(lapply(dirs, function(d) file.path(d, dir(d, '*.rds'))))
    dat <- lapply(fns, readRDS)
    do.call(rbind, dat)
}

getRDA <- function(rr1)
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
powSIM <- function(recache = FALSE)
{
    ## try the cached report first
    rds <- 'dat/sim_pwr.rds'
    if(file.exists(rds) && !recache)
        return(invisible(readRDS(rds)))
 
    ## get sumulation report first
    d0 <- getSIM()
    
    ## power calculation
    pw <- function(x, t = 0.05) sum(x < t) / length(x)
    fw <- pvl ~ ssz + vtx + knl + src + typ + adj
    d1 <- aggregate(formula = fw, FUN = pw, data = d0)

    ## number of repetition
    rp <- aggregate(formula = fw, FUN = length, data = d0)$pvl
    d1 <- within(d1, rep <- rp)
    
    ## either unajusted or FDR adjusted p value
    d1 <- subset(d1, adj %in% c('N', 'F'))
    d1 <- within(d1,
    {
        pwr <- pvl
        rm(pvl, adj)
    })

    ## reorder power source
    d1 <- within(d1,
    {
        odr <- c(N=1, G=2, V=3, A=4, I=5)[src]
        src <- as.factor(src)
        src <- reorder(src, odr)

        odr <- c(N=1, G=2, V=3, J=4)[knl]
        knl <- as.factor(knl)
        knl <- reorder(knl, odr)

        vtx[is.na(vtx)] <- 'NL'
        odr <- c(E0=1, E4=2, B2=3, NL=4)[vtx]
        vtx <- as.factor(vtx)
        vtx <- reorder(vtx, odr)
        
        rm(odr)
    })

    ## save to the cache and return
    saveRDS(d1, rds)
    invisible(d1)
}

.pic.pow.basic <- function(d)
{
    ## basic plot elements
    g <- ggplot(d, aes(x = ssz, y = pwr, group = paste(knl, vtx)))
    g <- g + xlab('sample size') + ylab('power') + xlim(100, 800) + ylim(0, 1)
    g <- g + theme(legend.position = "bottom", legend.box = "horizontal")

    ## U kernel composition is represented by point
    g <- g + geom_point(aes(shape = knl), size = 2L)
    g <- g + scale_shape_discrete(
        name = "",
        breaks = c("J", "G", "V"),
        labels = c("Joint", "Genomic", "Vertices"))
    g <- g + guides(
        shape = guide_legend(
            title = 'U kernels',
            label.position = 'bottom',
            label.hjust = 0.5))

    ## theme of facet titles
    g + theme(strip.text.x = element_text(family = 'times'))
    g
}

.pic.pow.facet <- function(type = c('Continuous', 'Dichotomous'))
{
    tp <- match.arg(type, c('Continuous', 'Dichotomous'))
    ## use facets to separate effect types
    lb <- function(labels, multi_line = TRUE)
    {
        .e <- expression
        if(tp == 'Continuous')
            ef = list(
                G = .e(Y[G] == G + epsilon),
                V = .e(Y[V] == V + epsilon),
                A = .e(Y[A] == G + V + epsilon),
                I = .e(Y[I] == G + V + G * symbol("*") * V + epsilon))
        else
            ef = list(
                G = .e(Pr(Y[G] ==1) == logit^-1 * (G + epsilon)),
                V = .e(Pr(Y[V] ==1) == logit^-1 * (V + epsilon)),
                A = .e(Pr(Y[A] ==1) == logit^-1 * (G + V + epsilon)),
                I = .e(Pr(Y[I] ==1) == logit^-1 * (G + V + G * symbol('*') * V + epsilon)))
        
        re1 <- list(ef[unlist(labels)])
        re1
    }
    
    r <- facet_wrap(~ src, NULL, NULL, labeller = lb)
    r
}

## compare three types of U kernel consitution
pic.KNL <- function(pwr, phe.type = 'C')
{
    library(ggplot2)

    ## continuous response, original vertices, regional test
    d <- subset(pwr, typ == phe.type & !vtx %in% c('E4', 'B2'), -c(typ))
    
    ## basic plot elements
    g <- .pic.pow.basic(d)
    g <- g + geom_line()
    
    ## divide into facets by effect composition
    g <- g + .pic.pow.facet(phe.type)
    g
}

## compare region test with vertex-wise analysis
pic.VWA <- function(pwr, phe.type = 'C')
{
    library(ggplot2)
    ## continuous response, original vertices, regional test
    d <- subset(pwr, typ == phe.type & vtx != 'E4' & knl != "G", -c(typ))
    
    ## basic plot elements
    g <- .pic.pow.basic(d)

    ## vertex kernel is represented type by line type
    g <- g + geom_line(aes(linetype = vtx))
    g <- g + scale_linetype_discrete(
        name = "",
        breaks = c("E0", "B2"),
        labels = c("signal\naggregation", "vertex-wise\nanalysis"))
    g <- g + guides(
        linetype = guide_legend(
            title = 'algorithm',
            label.position = 'bottom',
            label.hjust = 0.5))

    ## facets for effect composition
    g <- g + .pic.pow.facet(phe.type)
    
    g
}

pic.SAE <- function(pwr, phe.type = 'C')
{
    library(ggplot2)

    ## continuous response, original/encoded vertices, regional test
    d <- subset(pwr, typ == phe.type & vtx != 'B2' & knl != "G", -c(typ))
    
    ## basic plot elements
    g <- .pic.pow.basic(d)

    ## U kernel is represented by point, vertex type by line
    g <- g + geom_line(aes(linetype = vtx))
    g <- g + scale_linetype_discrete(
        name = "",
        breaks = c("E0", "B2", "E4"),
        labels = c("original", "original", "encoded"))
    g <- g + guides(
        linetype = guide_legend(
            title = 'type of\nvertex',
            label.position = 'bottom',
            label.hjust = 0.5))

    ## use facets to separate effect types
    g <- g + .pic.pow.facet(phe.type)
    g
}

## picture of power from simulation reports
picSIM <- function(pwr)
{
    library(ggplot2)
    lp <- list(
        PWR.CNT.KNL = pic.KNL(pwr, 'C'),
        PWR.CNT.VWA = pic.VWA(pwr, 'C'),
        PWR.CNT.SAE = pic.SAE(pwr, 'C'),
        PWR.BIN.KNL = pic.KNL(pwr, 'D'),
        PWR.BIN.VWA = pic.VWA(pwr, 'D'),
        PWR.BIN.SAE = pic.SAE(pwr, 'D'))

    for(nm in names(lp))
    {
        fn <- paste('rpt/img/', gsub('[.]', '_', nm), '.png', sep = '')
        ggsave(fn, lp[[nm]])
    }

    invisible(lp)
}

main <- function()
{
    rpt <- getSIM(T)
    pwr <- powSIM(rpt)
    pic <- picSIM(pwr)
    invisible(pic)
}

qqplot <- function(dat)
{
    library(ggplot2)
    dat <- by(dat, dat[, c('alg', 'wgt')], within,
    {
        ## negative log 10 of p-values, actual
        lpvl1 <- -log10(pvl)
        ## negative log 10 of p-values, theoratical
        lpvl0 <- -log10(seq(0, 1, length.out = length(pvl) +1 ))[-1]

        label <- sprintf('ALG=%s, KNL=%s, W=%s', alg, krn, wgt)
    })
    dat <- lapply(dat, na.omit)
    dat <- do.call(rbind, dat)

    dat$lpvl1 <- with(dat,
    {
        ix <- is.infinite(lpvl1)
        mx <- max(lpvl1[!ix]) + 1
        lpvl1[ix] <- mx
        lpvl1
    })
    x.rng <- with(dat, c(0, max(lpvl1)))
    
    nc <- sqrt(length(unique(dat$label)))
    qp <- ggplot(dat)
    qp <- qp + geom_point(aes(x=lpvl0, y = lpvl1))
    qp <- qp + xlab(expression(Theoretical~~-log[10](italic(p))))
    qp <- qp + ylab(expression(Observed~~-log[10](italic(p))))
    qp <- qp + geom_abline(slope = 1, intercept = 0)
    qp <- qp + facet_wrap(~ label, ncol = ceiling(nc))
    qp <- qp + xlim(x.rng)
    qp
}
