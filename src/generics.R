## pick data samples from source directory
pck <- function(
    src, size = 1, replace = FALSE, seed = NULL,
    drop = TRUE, vbs = FALSE,
    ret = c('file', 'rds', 'img', 'gno'))
{
    ## pick out filenames
    fns <- file.path(src, dir(src, '*.rds'))
    set.seed(seed)
    if(replace | size < length(fns))
        fns <- sample(fns, size, replace)
    set.seed(NULL)

    ## if only requests file nemas to be returned
    ret <- match.arg(ret)
    if(ret == 'file')
        return(fns)

    env <- findFunction(ret)
    if(length(env) > 0L)
        typ <- getFunction(ret, where = env[[1L]])
    else
        typ <- readRDS
    
    dat <- sapply(fns, typ, simplify = F, USE.NAMES = F)
    if(drop & length(dat) < 2L)
        return(dat[[1]])
    dat
}

maf <- function(obj, ...) UseMethod("maf")
maf.matrix <- function(gmx, margin = 1L, level = 2L)
{
    f <- if(margin == 1L) rowMeans else colMeans
    f(gmx, na.rm = T) / level
}
maf.vector <- function(snp, level = 2L)
{
    mean(snp, na.rm = T) / level
}
## initialize object
ini <- function(obj) UseMethod("ini")

## pck a subset of subjects from a data sample
sbj <- function(dat, IDs) UseMethod("sbj")
