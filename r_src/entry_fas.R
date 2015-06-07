source('r_src/gno.R');
source('r_src/utl.R');
source('r_src/hwu.R')

img.pck <- function(src, idv)
{
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
    img = img.pck('dat/enc', 'dat/ssn')
    imx = img$img
    
    ## pick genomic segment
    gno <- GNO$pck(
        vcf = 'dat/wgs/c03.vcf.gz', idv='dat/wgs/idv.EUR',
        seg='dat/wgs/gen', wnd=5000, n = 1)

    ## some gene has no valid variants
    if(length(gno) < 1)
        return(NA)

    gno = gno[[1]]
    gno = GNO$clr(gno)
    gno = GNO$imp(gno)
    gmx = gno$gmx[gno$idx,]

    ## pick the same number of samples
    N1 = ncol(gmx)
    N2 = nrow(imx)
    N = min(N1, N2)

    ix1 = sample.int(N1, N)
    ix2 = sample.int(N2, N)
    gmx = gmx[, ix1]
    imx = imx[ix2, ]
    
    ## make binary phenotype and normal covariate
    y <- rbinom(n = N, size = 1, prob = 0.35)
    x <- rnorm(n = N, mean = 0, sd = 1.5)
    
    g <- HWU$collapse.burden(t(gmx))
    wg <- HWU$weight.gaussian(g)

    i <- HWU$collapse.burden(imx)
    wi <- HWU$weight.gaussian(i)
    
    out.g <- HWU$dg2(y, x, wg);
    out.g
}

set.seed(2)
main <- function()
{
    p = replicate(n = 1000, type1())
    p
}
