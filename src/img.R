#!/usr/bin/env Rscript
source('src/vwa.R')

ini.img <- function(img)
{
    ## append dimension names to the surface data
    names(dimnames(img$sfs)) <- c('ftr', 'vtx', 'sbj')

    ## basic decoration
    img <- within(
        img,
    {
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
    img <- ini.vwa(img)
    img
}

## randomly pick encoded image data from a folder
pck.img <- function(
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

pck.sbj.img <- function(img, sbj)
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
        gsb <- gsb[, , , I]
    })

    ## in case sb. think {sbj} means sample size instead of
    ## the IDs of wanted subject
    if(length(img$sbj) == 0L)
        warning('no subject ID matches image data source.')
    img
}

.az.img.ec2 <- Sys.getenv('AZ_EC2')
.az.img <- Sys.getenv('AZ_SM2')         # 1/2 encoding
cml.img <- function()
{
    argv <- commandArgs(trailingOnly = TRUE)
    if(length(argv) < 1L)
        return(NULL)
    
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
        default='ini.img')
    
    opt <- parse_args(p, argv)
    attach(opt)
    act <- getFunction(act)
    img <- readRDS(file.path(src, paste(ssn, 'rds', sep='.')))
    img <- act(img)
    saveRDS(img, file.path(dst, paste(ssn, 'rds', sep='.')))
}

## run the command line parser
opt <- cml.img()
