source('src/hlp.R')
library(ggplot2)

cat.rpt <- function(src, ...)
{
    ## pick out images by file name
    dirs <- c(src, ...)
    
    fns <- unlist(lapply(dirs, function(d) file.path(d, dir(d, '*.rds'))))
    dat <- lapply(fns, readRDS)
    do.call(rbind, dat)
}

## get real data analysis report
getRDA <- function(recache = F)
{
    ## check the cached version
    rs <- 'dat/rda_clr.rds'
    if(!recache && file.exists(rs))
        return(readRDS(rs))

    ## clean up
    d0 <- readRDS('dat/rda_raw.rds')
    d0 <- na.omit(d0)

    d0 <- with(d0,
    {
        wnm <- paste(sub('h.*$', '', wsn), wnm, sep = '.')
        data.frame(
            row.names = paste(wsn, gsn, sep='.'),
            GEN = as.factor(gnm),       # gene
            CTX = as.factor(wnm),       # cortex
            NV = n.v,                   # number of vertices
            NG = n.g,                   # number of genomic variants
            PG = E4.G,                  # p value of U_G
            PV = E4.V,                  # p value of U_V
            PJ = E4.X)                  # p value of U_J
    })
    d0 <- with(d0, d0[order(PJ),])

    ## multiple testing correction
    d0 <- within(d0,
    {
        BJ <- pmin(1, PJ * length(PJ))
        BV <- pmin(1, PV * length(levels(CTX)))
        BG <- pmin(1, PG * length(levels(GEN)))

        FJ <- p.adjust(PJ, 'fdr')
        FV <- p.adjust(PV, 'fdr')
        FG <- p.adjust(PG, 'fdr')
    })

    ## save to cache and return
    saveRDS(d0, rs)
    d0
}

## picture of real data analysis
pic.RDA.PVL <- function(dt, np = 500L, out = NULL)
{
    dt <- dt[seq(1, nrow(dt), l=np), ]
    dt <- within(dt,
    {
        LG <- -log10(PG)
        LV <- -log10(PV)
        LJ <- -log10(PJ)
        rm(PG, PV, PJ)
    })
    
    if(!is.null(out))
        png(out, width=2400, height=1200, res=300, pointsize = 10)

    ## start the plot
    par(cex.lab = 1.3, cex.axis = 1.3, mar = c(4,4,1,1), mgp = c(2.5, 1, 0))
    y.lm <- with(dt, c(0, max(LG, LV, LJ)))
    plot("",
         xlim = c(0, 100),
         xlab = expression(rank(-log[10]*P[J])*'%'),
         ylim = y.lm,
         ylab = expression(-log[10] * P))
    
    abline(h = pretty(y.lm), col = "lightgray", lwd = 1)
    with(dt,
    {
        x <- 1:nrow(dt) / nrow(dt) * 100
        points(x, LG, pch = 43)
        points(x, LV, pch = 23)
        points(x, LJ, pch = 20)
    })
    lgd <- c(
        G=expression(-log[10]*P[G]),
        V=expression(-log[10]*P[V]),
        J=expression(-log[10]*P[J]))
    legend('topright', 
        legend=lgd, pch=list(43, 23, 20), pt.cex = 1.5)

    if(!is.null(out))
        dev.off()
}

pic.RDS.qq <- function(dt, np = 500L)
{
    dt <- dt[seq(1, nrow(dt), l=np), ]
    dt <- within(dt,
    {
        LG <- -log10(PG)
        LV <- -log10(PV)
        LJ <- -log10(PJ)
        rm(PG, PV, PJ)
    })

    library(ggplot2)
    qp <- ggplot(dt)
    qp <- qp + xlab(expression(-log[10](italic(p[J]))))
    qp <- qp + ylab(expression(-log[10](italic(p[V]))))
    qp <- qp + geom_point(aes(x = LJ, y = LV), shape = 1)
    qp <- qp + geom_abline(slope = 1, intercept = 0)
    ##    qp <- qp + facet_wrap(~ label, ncol = ceiling(nc))
    qp <- qp + ylim(0, max(dt$LJ))
    qp
}

## mark statistical significance
mk <- function(pvl, fdr, bon)
{
    p <- format(pvl, digits = 3, scientific = T)
    m <- character(length(pvl))
    f <- fdr < 0.01
    b <- bon < 0.05
    m[f] <- paste(m[f], "_+", sep = "")
    m[b] <- paste(m[b], "^*", sep = "")
    sprintf('$%s%s$', p, m)
}
mk.note <- function()
{
    '\\hline
    \\multicolumn{7}{l}{\\texttt{*: below 0.05 after Bonferroni correction}} \\\\ \n
    \\multicolumn{7}{l}{\\texttt{+: below 0.01 after FDR correction}}        \\\\'
}
mk.head <- function()
{
    'GENE & CORTEX & $|V|$ & $|G|$ & $P_G$ & $P_V$ & $P_J$'
}

## top 20 most significant
tab.RDA.T20 <- function(dt, out = "rpt/tbl/RDA_T20.tex")
{
    dt <- subset(dt, PJ < pmin(PG,PV))[1:20, ]
    dt <- within(dt,
    {
        PG <- mk(PG, FG, BG)
        PV <- mk(PV, FV, BV)
        PJ <- mk(PJ, FJ, BJ)
        rm(BG, BV, BJ, FG, FV, FJ)    
    })
    
    rownames(dt) <- NULL
    colnames(dt) <- c('GENE', 'CORTEX', '$|V|$', '$|G|$',
                      '$\\qquad P_G$', '\\qquad $P_V$', '\\qquad $P_J$')
    
    library(xtable)
    ds <- strsplit('sssddgee', '')[[1]]
    al <- strsplit('lllcclll', '')[[1]]
    tb <- xtable(dt, digits = 3, align = al, display = ds)
    atr <- list(pos = list(nrow(tb)), command = mk.note()) # add to rows

    print(tb, file=out, floating = F, include.rownames=F, 
          sanitize.text.function = identity,
          add.to.row = atr)
}

## top 20 most significant, per cortex region
tab.RDA.JNT <- function(dt, out = "rpt/tbl/RDA_JNT.tex")
{
    dt <- subset(dt, PJ < pmin(PG,PV))
    
    ## split to every cortex region
    dt <- by(dt, dt$CTX, function(ctx)
    {
        ctx[1,]
    })
    dt <- do.call(rbind, dt)
    dt <- with(dt, dt[order(PJ), ])[1:20, ]

    dt <- within(dt,
    {
        PJ <- mk(PJ, FJ, BJ)
        PV <- mk(PV, FV, BV)
        PG <- mk(PG, FG, BG)
        rm(BG, BV, BJ, FG, FV, FJ)
    })
    
    rownames(dt) <- NULL
    colnames(dt) <- c('GENE', 'CORTEX', '$|V|$', '$|G|$',
                      '$\\qquad P_G$', '\\qquad $P_V$', '\\qquad $P_J$')
    
    library(xtable)
    ds <- strsplit('sssddgee', '')[[1]] # display
    al <- strsplit('lllcclll', '')[[1]] # alignment
    tb <- xtable(dt, digits = 3, align = al, display = ds)
    atr <- list(pos = list(nrow(tb)), command = mk.note()) # add to rows

    print(tb, file = out, include.rownames = F, floating = F,
          sanitize.text.function = identity,
          add.to.row = atr)
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
            src = unname(c(G='G', V='V', A='A', X='I', N='N')[cf[4]]),
            typ = unname(c(L='C', B='D')[cf[5]]),
            adj = cf[1],
            pvl = d0[, x],
            stringsAsFactors = FALSE)
        rt
    }, simplify = F, USE.NAMES = F)
    d0 <- do.call(rbind, d0)
    d0 <- within(d0,
    {
        vtx[knl == 'G' ] <- 'NL'        # a gnomic test needs no vertex
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
    
    ## either unajusted or fdr adjusted p value
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
    .e <- expression
    g <- g + geom_point(aes(shape = knl), size = 2.5)
    g <- g + scale_shape_discrete(
        name = "  Statistics:",
        breaks = c("J", "G", "V"),
        labels = c("joint  ", "genomic  ", "imaging  "))

    g <- g + guides(
        shape = guide_legend(label.hjust = 0.5))

    ## theme of facet titles
    g <- g + theme(strip.text.x = element_text(family = 'times'))
    g
}

.pic.pow.facet <- function(type = c('Continuous', 'Dichotomous'))
{
    .t <- match.arg(type, c('Continuous', 'Dichotomous'))
    .e <- expression
    
    ## use facets to separate effect types
    if(.t == 'Continuous')
        .l <- list(
            N = .e(Y[0] == epsilon),
            G = .e(Y[G] == G + epsilon),
            V = .e(Y[V] == V + epsilon),
            A = .e(Y[A] == G + V + epsilon),
            I = .e(Y[I] == G + V + G * symbol("*") * V + epsilon))
    else
        .l <- list(
            N = .e(Pr(Y[0] ==1) == epsilon),
            G = .e(Pr(Y[G] ==1) == logit^-1 * (G + epsilon)),
            V = .e(Pr(Y[V] ==1) == logit^-1 * (V + epsilon)),
            A = .e(Pr(Y[A] ==1) == logit^-1 * (G + V + epsilon)),
            I = .e(Pr(Y[I] ==1) == logit^-1 * (G + V + G * symbol('*') * V + epsilon)))

    lb <- function(labels, multi_line = TRUE)
    {
        list(.l[unlist(as.character(unlist(labels)))])
    }
    
    r <- facet_wrap(~ src, NULL, NULL, labeller = lb)
    r
}

## compare three types of U kernel consitution
pic.KNL <- function(pwr, phe.type = 'C')
{
    ## continuous response, original vertices, regional test
    d <- subset(pwr, typ == phe.type & !vtx %in% c('E4', 'B2'), -c(typ))
    
    ## basic plot elements
    g <- .pic.pow.basic(d)
    g <- g + geom_line()
    
    ## divide into facets by effect composition
    g <- g + .pic.pow.facet(phe.type)

    ## decorate legends
    g <- g + theme(legend.position = 'top', legend.margin = unit(0, 'cm'))
    g <- g + theme(legend.title = element_text(face="bold"))

    g
}

## compare region test with vertex-wise analysis
pic.VWA <- function(pwr, phe.type = 'C')
{
    ## continuous response, original vertices, regional test
    d <- subset(pwr, typ == phe.type & vtx != 'E4' & knl != "G", -c(typ))
    
    ## basic plot elements
    g <- .pic.pow.basic(d)

    ## vertex kernel is represented type by line type
    g <- g + geom_line(aes(linetype = vtx))
    g <- g + scale_linetype_discrete(
        name = '  Algorithm:',
        breaks = c("E0", "B2"),
        labels = c("grouping &\naggregation  ", "vertex-wise\nanalysis  "))

    ## facets for effect composition
    g <- g + .pic.pow.facet(phe.type)
    
    ## decorate legends
    g <- g + theme(legend.position = 'top', legend.margin = unit(0, 'cm'))
    g <- g + theme(legend.title = element_text(face="bold"))

    g
}

pic.SAE <- function(pwr, phe.type = 'C')
{
    ## continuous response, original/encoded vertices, regional test
    d <- subset(pwr, typ == phe.type & vtx != 'B2' & knl != "G", -c(typ))
    
    ## basic plot elements
    g <- .pic.pow.basic(d)

    ## U kernel is represented by point, vertex type by line
    g <- g + geom_line(aes(linetype = vtx))
    g <- g + scale_linetype_discrete(
        name = "  Imaging Profile:",
        breaks = c("E0", "B2", "E4"),
        labels = c("raw  ", "raw  ", "high order ft.  "))

    ## use facets to separate effect types
    g <- g + .pic.pow.facet(phe.type)

    ## decorate legends
    g <- g + theme(legend.position = 'top', legend.margin = unit(0, 'cm'))
    g <- g + theme(legend.title = element_text(face="bold"))

    g
}

## picture of power from simulation reports
picSIM <- function(pwr)
{
    library(ggplot2)
    pwr <- subset(pwr, src != 'N')
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
    sim <- getSIM()
    sim <- powSIM()
    pic <- picSIM(sim)

    rda <- getRDA()
    tab.RDA.T20(rda, out = 'rpt/tbl/RDA_T20.tex')
    tab.RDA.JNT(rda, out = 'rpt/tbl/RDA_JNT.tex')
    pic.RDA.PVL(rda, out = 'rpt/img/RDA_PVL.png')    

    NULL
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

    ## number of columns of the facets
    nc <- sqrt(length(unique(dat$label)))

    ## plot
    qp <- ggplot(dat)
    qp <- qp + geom_point(aes(x=lpvl0, y = lpvl1))
    qp <- qp + xlab(expression(Theoretical~~-log[10](italic(p))))
    qp <- qp + ylab(expression(Observed~~-log[10](italic(p))))
    qp <- qp + geom_abline(slope = 1, intercept = 0)
    qp <- qp + facet_wrap(~ label, ncol = ceiling(nc))
    qp <- qp + xlim(x.rng)
    qp
}
