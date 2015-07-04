source('src/gno.R')
source('src/utl.R')
source('src/hwu.R')
library('data.table')

## randomly pick encoded image data from a folder
img.pck <- function(src)
{
    if(!file.exists(src))
        stop(paste(src, "not exists"))

    imgs <- dir(src, '*.enc')
    img <- imgs[sample.int(n = length(imgs), size = 1)]
    vtx <- sub('.enc', '', x = img)
    idv <- sub('enc', 'sbj', x = img)

    img.file <- paste(src, img, sep = '/')
    sbj.file <- paste(src, idv, sep = '/')
    img <- scan(img.file, what = 0.0, quiet = T)
    sbj <- scan(sbj.file, what = " ", quiet = T)

    N = length(sbj)                     # subject count - col
    M = length(img)/N                   # feature count - row
    dim(img) = c(M, N)
    list(src=src, vtx=vtx, imx=img, sbj=sbj, M=M, N=N)
}

type1 <- function(N = NULL)
{
    ## pick encoded image data
    img = img.pck('hpc/encoded_wms/09_0078')
    imx <- img$imx
    
    ## pick genomic segment, and only pick those with image available
    gno <- GNO$pck(vcf = 'dat/wgs/gno', seg='dat/wgs/gen', sbj=img$sbj, wnd=5000, n = 1)[[1]]
    gmx <- gno$gmx
    
    ## randomly select common subject
    sbj = intersect(gno$sbj, img$sbj)
    if(!is.null(N))
    {
        N <- min(N, length(sbj))
    }
    else
    {
        N <- length(sbj)
    }
    sbj = sort(sbj[sample.int(length(sbj), size = N)])
    imx <- imx[, match(sbj, img$sbj)]
    gmx <- gmx[, match(sbj, gno$sbj)]
    
    ## clean up genotype
    gmx = GNO$clr.dgr(gmx)
    ## sometimes cleaup will uncheck all variants, we have
    ## to skip this type I trial
    if(length(gmx) == 0)
        return(NA)

    ## guess missing values
    gmx = GNO$imp(gmx)
    
    ## tranpose data matrices to design matrices
    gmx = t(gmx)
    imx = t(imx)
    stopifnot(nrow(gmx) == nrow(imx))
    
    ## make one binary phenotype and normal covariate
    y <- rbinom(n = N, size = 1, prob = 0.5)
    x <- rnorm(n = N, mean = 0, sd = 1.5)
    out <- list()
    
    ## g <- hwu.collapse.burden(gmx)
    ## wg <- hwu.weight.gaussian(gmx)
    wg <- hwu.weight.IBS(gmx)
    
    ## i <- hwu.collapse.burden(imx)
    ## wi <- hwu.weight.gaussian(imx)
    wi <- hwu.weight.gaussian(imx)

    out$g1IBS.i1GSN <- hwu.dg2(y, x, wg, wi)
    
    wg <- hwu.weight.IBS(gmx, w = hwu.w.MAFsd(gmx))
    out$gWIBS.i1GSN <- hwu.dg2(y, x, wg, wi)

    out
}

main <- function(n = 1000)
{
    p = replicate(n, type1(50), simplify = F)
    p
}
