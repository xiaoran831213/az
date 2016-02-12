#!/usr/bin/env Rscript
source('src/vwa.R')
source('src/hwu.R')
source('src/generics.R')

## * -------- generics -------- * ##
## pick or list subjects of a data
sbj <- function(obj, ...) UseMethod("sbj")

## pick or list vertices of a data
vtx <- function(obj, ...) UseMethod("vtx")

## number of subjects
nsb <- function(obj, ...) UseMethod("nsb")

## number of vertices
nvt <- function(obj, ...) UseMethod("nvt")

## number of features
nft <- function(obj, ...) UseMethod("nft")
## * --------  (end)  -------- * ##

## the constructor: vetex based image
vimage <- function(x)
{
    if(is.null(x))
    {
        x <- list(
            cmx = NULL,                 # connection matrix
            enc = NULL,                 # vertex encoding 
            gsb = NULL,                 # gaussian blur
            nwk = NULL,                 # network
            nwt = NULL,                 # net weight
            sbj = NULL,                 # subject ids
            sfs = NULL,                 # surfaces
            vtx = NULL)                 # vertex ids
    }
    else if(is.character(x))            # read from file
    {
        x <- readRDS(x)
    }
    structure(x, class=c('vimage', 'list'))
}

## number of subjects
nsb.vimage <- function(vmg, ...) length(dimnames(vmg$sfs)$sbj)

## number of vertices
nvt.vimage <- function(vmg, ...) length(dimnames(vmg$sfs)$vtx)

## extend default dim, also nrow and ncol would work
dim.vimage <- function(vmg) c(nvt.vimage(vmg), nsb.vimage(vmg))

## features
ftr.vimage <- function(vmg) dimnames(vmg$sfs)$ftr

str.vimage <- function(img)
{
    dnm <- 
    ftr <- paste(dimnames(img$sfs)$ftr, collapse=' ,')
    nvt <- nvt(img)
    nsb <- nsb(img)
    fmt <-  'ftr=%s; n.v=%d; n.s=%d'
    sprintf(fmt, ftr, nvt, nsb)
}

print.vimage <- function(img)
{
    print(str.vimage(img))
    print(ls(img))
}

ini.vimage <- function(img)
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

## list or extract subjects
sbj.vimage <- function(img, who = NULL)
{
    if(is.null(who))
        with(img, dimnames(img$sfs)$sbj)
    else
    {
        who <- unlist(who)
        if(is.character(who))
            I <- match(who, dimnames(img$sfs)$sbj)
        within(img,
        {
            sfs <- sfs[, , I]
            if('enc' %in% ls())
                enc <- lapply(enc, function(u) u[I, ])
            if('gsb' %in% ls())
                gsb <- gsb[, , , I]
        })
    }
}

## list or extract vertices
vtx.vimage <- function(vmg, who = NULL)
{
    if(is.null(who))
        with(vmg, dimnames(vmg$sfs)$vtx)
    else
    {
        who <- unlist(who)
        if(is.character(who))
            who <- match(who, dimnames(vmg$sfs)$sbj)
        within(vmg,
        {
            sfs <- sfs[, who, ]
            if('gsb' %in% ls())
                gsb <- gsb[, who, ,]
            ## leave enc intact
        })
    }
}

.az.img.ec2 <- Sys.getenv('AZ_EC2')
.az.img <- Sys.getenv('AZ_SM2')         # 1/2 encoding
