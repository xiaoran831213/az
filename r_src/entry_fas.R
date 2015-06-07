source('r_src/gno.R');
source('r_src/utl.R');
source('r_src/hwu.R')

## wnd --- segment window size
## cnd --- number of candidate segments to pick
## vcf --- VCF file to read. (Varient Call Format)
## idv --- IDV file to read. (Individual)
## seg --- SEG file to read. (Segment file, e.g. Genes)
type1 <- function()
{
    ## pick genomic segment
    gno <- GNO$pck(
        vcf = 'dat/wgs/c03.vcf.gz', idv='dat/wgs/idv.EUR',
        seg='dat/wgs/gen', wnd=0, n = 1)

    print(length(gno))
    gno = gno[[1]]
    
    gno = GNO$clr(gno)
    gno = GNO$imp(gno)
    
    gmx = gno$gmx[gno$idx,]
    N = ncol(gmx)
    
    ## make binary phenotype and normal covariate
    y <- rbinom(n = N, size = 1, prob = 0.35)
    x <- rnorm(n = N, mean = 0, sd = 1.5)
    
    g <- HWU$collapse.burden(t(gmx));
    wg <- HWU$weight.gaussian(g);
    
    out.g <- HWU$dg2(y, x, wg);
    out.g
}

set.seed(2)
main <- function()
{
    p = replicate(n = 1000, type1())
    p
}
