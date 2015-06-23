source('src/gno.R')
source('src/utl.R')
source('src/hwu.R')
library('data.table')

## randomly pick processed image data from a folder
img.pck <- function(src, idv)
{
    if(!file.exists(src))
        stop(paste(src, "not exists"))
    
    idv <- read.table(idv, sep = '\t', as.is = T)
    N <- nrow(idv)
    
    img <- dir(src)
    img <- img[sample.int(n = length(img), size = 1)]
    cat(img, '\t')
    img <- paste(src, '/', img, sep = '')

    cnn <- file(img, "rb");
    img <- matrix(
        scan(cnn, what=0.0, quiet=T),
        nrow = nrow(idv), byrow = T);
    close(cnn);
    list(idv = idv, img = img)
}

type1 <- function()
{
    ## pick encoded image data
    img = img.pck('dat/img/enc', 'dat/img/ssn')
    imx = img$img
    
    ## pick genomic segment
    gno <- GNO$pck(
        vcf = 'dat/wgs/gno', idv='dat/wgs/idv.cmn',
        seg='dat/wgs/gen', wnd=5000, n = 1)[[1]]

    gno = GNO$clr(gno)
    gno = GNO$imp(gno)

    ## sometimes cleaup will cross all variants
    if(!is.null(gno$idx) && length(gno$idx) < 1)
        return(NA)

    ## drop = F will keep matrix structure even if
    ## on of the dimension is reduced to one
    gmx = gno$gmx[gno$idx,, drop = F]
    

    ## pick the same number of samples
    N1 = ncol(gmx)
    N2 = nrow(imx)
    N = min(N1, N2)

    ix1 = sample.int(N1, N)
    ix2 = sample.int(N2, N)
    gmx = gmx[, ix1, drop = F]
    imx = imx[ix2, , drop = F]
    
    ## make binary phenotype and normal covariate
    y <- rbinom(n = N, size = 1, prob = 0.5)
    x <- rnorm(n = N, mean = 0, sd = 1.5)
    
    g <- HWU$collapse.burden(t(gmx))
    wg <- HWU$weight.gaussian(g)

    i <- HWU$collapse.burden(imx)
    wi <- HWU$weight.gaussian(i)
    
    out.g <- HWU$dg2(y, x, wg);
    out.g
}

main <- function(n = 1000)
{
    p = replicate(n, type1())
    p
}
