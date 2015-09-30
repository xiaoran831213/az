#!/usr/bin/env Rscript
source('src/vwa.R')
source('src/generics.R')

img <- function(x)
{
    if(!is.list(x))
        x <- readRDS(x)
    structure(x, class=c('img', 'list'))
}

str.img <- function(img)
{
    dnm <- 
    ftr <- paste(dimnames(img$sfs)$ftr, collapse=' ,')
    n.v <- length(img$vtx)
    n.s <- length(img$sbj)
    fmt <-  'ftr=%s; n.v=%d; n.s=%d'
    sprintf(fmt, ftr, n.v, n.s)
}

print.img <- function(img)
{
    print(str.img(img))
}

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

sbj.img <- function(img, IDs)
{
    I <- match(IDs, img$sbj)
    img <- within(
        img,
    {
        sbj <- sbj[I]
        if('sfs' %in% ls())
            sfs <- sfs[, , I]
        if('enc' %in% ls())
            enc <- lapply(enc, function(u) u[I, ])
        if('gsb' %in% ls())
            gsb <- gsb[, , , I]
    })

    ## in case sb. think {sbj} means sample size instead of
    ## the IDs of wanted subject
    if(length(img$sbj) == 0L)
        warning('no subject ID matches image source.')
    img
}

.az.img.ec2 <- Sys.getenv('AZ_EC2')
.az.img <- Sys.getenv('AZ_SM2')         # 1/2 encoding
