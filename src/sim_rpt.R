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

pic <- function(pwr, xt='n[.]s', et='012345678', wt="VGX", yt="VGAX")
{
    et <- paste('[', et, ']', sep='')
    wt <- paste('[', wt, ']', sep='')
    yt <- paste('[', yt, ']', sep='')

    nm <- names(pwr)
    yt <- paste(et, wt, yt, sep='[.]')

    xns <- nm[grepl(xt, nm)]
    yns <- nm[grepl(yt, nm)]

    
    xn <- xns[1]
    xv <- pwr[,xn]
          
    par(mfcol=c(1, length(yns)))
    for(j in 1L:length(yns))
    {
        yn <- yns[j]
        yv <- pwr[,yn]
        plot(x=xv, y=yv, main = 'Power v.s. N', type = 'l', ylim = c(0,1))
        
    }
    par(mfcol=c(1,1))


    #lines(x=xv, y=yv, type = 'l')

    #dev.off()
}
