#!/usr/bin/env Rscript
source('src/vwa.R')

.ini.img <- function(src, ssn, dst = NULL)
{
    img <- readRDS(file.path(src, paste(ssn, 'rds', sep='.')))

    ## append dimension names to the surface data
    names(dimnames(img$sfs)) <- c('ftr', 'vtx', 'sbj')

    ## basic decoration
    img <- within(
        img,
    {
        ssn <- ssn
        names(dimnames(cmx)) <- list('a', 'b')
        sbj <- dimnames(sfs)$sbj
        vtx <- dimnames(sfs)$vtx
        enc <- lapply(enc, function(u)
        {
            dimnames(u) <- list(sbj=sbj, vtx=sprintf('C%04X', 1L:ncol(u)))
            u
        })
    })

    ## preparation for vertex wise analysis
    img <- .vwa.ini(img)                # basic
    img <- .vwa.gsb(img)                # gaussian blur

    ## cache & return
    if(!is.null(dst))
       saveRDS(img, file.path(dst, paste(ssn, 'rds', sep='.')))
    img
}

## randomly pick encoded image data from a folder
img.pck <- function(
    src, size = 1, replace = FALSE, seed = NULL,
    drop = TRUE, vbs = FALSE, ret = c('data', 'file'), recache = FALSE)
{
    ## pick out images by file name
    fns <- file.path(src, dir(src, '*.rds'))
    set.seed(seed)
    if(replace | size < length(fns))
        fns <- sample(fns, size, replace)
    set.seed(NULL)
    
    ## if only requests file nemas to be returned
    if(match.arg(ret) == 'file')
        return(fns)

    ims <- sapply(fns, readRDS, simplify = F, USE.NAMES = F)
    if(drop & length(ims) < 2L)
        return(ims[[1]])
    ims
}

img.sbj.pck <- function(img, sbj)
{
    I <- match(sbj, img$sbj)
    img <- within(
        img,
    {
        sbj <- sbj[I]
        sfs <- sfs[, , I]
        enc <- lapply(enc, function(u)
        {
            u[I, ]
        })
        ## remove location to prevent accidental overwite
        fgs <- union(fgs, 'sbj.pck')
    })

    ## in case sb. think {sbj} means sample size instead of
    ## the IDs of wanted subject
    if(length(img$sbj) == 0L)
        warning('no subject ID matches image data source.')
    img
}

.az.ec2 <- Sys.getenv('AZ_EC2')         # 1/2 encoding
.az.ec3 <- Sys.getenv('AZ_EC3')         # super fitted encoding
.az.ec4 <- Sys.getenv('AZ_EC4')         # 3/4 encoding
.az.ec5 <- Sys.getenv('AZ_EC5')         # 2/3 encoding

.cml.img <- function(argv = NULL, ...)
{
    library(argparser)
    p <- arg_parser('AZ image processing.')
    p <- add_argument(
        p, 'src', help = 'source directory of surface archives.')
    p <- add_argument(
        p, 'ssn', help = 'name of the surface sample to be processed')
    p <- add_argument(
        p, '--dst', help = 'target directory to store processed images.',
        default = '.')
    p <- add_argument(
        p, '--act', help = 'action to be done with the sample.',
        default='ini')

    ## if any argument is supplied by code or via R console, a parsing
    ## test is scheduled, otherwise, take arguments from command line
    argv <- unlist(c(argv, ...))
    if(is.null(argv))
        argv <- commandArgs(trailingOnly = TRUE)
    opt <- parse_args(p, argv)

    opt <- within(
        opt,
    {
        act <- getFunction(paste('', act, 'img', sep='.'))
    })

    with(opt, act(src, ssn, dst))
}

## run the command line parser
opt <- .cml.img()
