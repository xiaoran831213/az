library(igraph, quietly = T, warn.conflicts = F)
library(Matrix)
library(parallel)
source('src/hwu.R')

ini.vwa <- function(img)
{
    ## gaussian blur
    if(is.null(img$gsb))
    {
        ## vertex coordinate for each subject
        xyz <- aperm(img$sfs[c('x', 'y', 'z'), ,], c(2, 1, 3))
        vtx <- img$vtx

        ## get non-zero entries in the vertex connection matrix
        cnn <- subset(summary(img$cmx), select = c(i,j))
        i <- cnn$i                          # from v.index
        j <- cnn$j                          # to v.index
        
        ## the euclidean distance between adjecent vertices,
        ## will be the weight of the edges in the network.
        nwt <- apply(xyz, 3L, function(p)
        {
            sqrt(rowSums((p[i, ] - p[j, ]) ^ 2L))
        })

        ## replace vertex index with id
        i <- vtx[i]                         # from v.id
        j <- vtx[j]                         # to v.id

        ## assign names to network weight list
        dimnames(nwt) <- list(vtx=paste(i, j, sep = '.'), sbj=img$sbj)
        img$nwt <- nwt
        
        ## replace vertex index with id, create the graph
        img$nwk <- graph_from_data_frame(data.frame(i, j), directed = F)

        ## gaussian blur
        img$gsb <- .gsb.vwa(img)
    }
    img
}

.gsb.vwa <- function(img)
{
    library(abind)
    vtx <- img$vtx
    sbj <- img$sbj
    
    ## extract fatures, scale to [0, 1]
    ftr <- c('slc', 'tck')
    sfs <- apply(img$sfs[ftr, , , drop = F], 'ftr', function(v)
    {
        (v - min(v)) / (max(v) - min(v))    
    })
    dim(sfs) <- c(length(vtx), length(sbj), length(ftr))
    dimnames(sfs) <- list(vtx=vtx, sbj=sbj, ftr=ftr)
    sfs <- aperm(sfs, c('ftr', 'vtx', 'sbj'))

    ## sanity check
    stopifnot(all(t(sfs['tck', , ])==img$enc$tck.0))
    stopifnot(all(t(sfs['slc', , ])==img$enc$slc.0))
    
    ## gaussian brush sizes
    sdv <- 2^(0:5)
    names(sdv) <- paste('sd', seq_along(sdv), sep='')
    
    ## for gaussian blur, go through subjects
    gsb <- sapply(img$sbj, function(k)      # subject k
    {
        w <- img$nwt[, k]               # k th. sbject weights
        f <- sfs[, , k]                 # k th. sbject features
        
        ## shortest inter vertex geodesic distance
        d <- shortest.paths(img$nwk, weights = w)[vtx, vtx] * 0.866
        
        ## gaussian blur
        sapply(sdv, function(j)         # the j.th SDV
        {
            d <- dnorm(d, 0, j)
            diag(d) <- nrow(d) * diag(d)
            d <- d / colSums(d)
            f %*% d
        }, simplify = 'array')
    }, simplify = 'array')
    
    ## bind the original surface as sdv0 (non blur)
    gsb <- abind(sd0=sfs, gsb, along = 3L)
    names(dimnames(gsb)) <- list('ftr', 'vtx', 'sdv', 'sbj')

    gsb
}
