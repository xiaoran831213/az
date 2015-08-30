source('src/gno.R')
source('src/utl.R')
source('src/hwu.R')
source('src/hlp.R')
library(igraph)

.vwa.ini <- function(img)
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
    img
}

.vwa.gbr <- function(img)
{
    ## gaussian brush sizes
    sdv <- 2^(0:4)
    names(sdv) <- paste('sd', sdv, sep='')

    ## for gaussian blur, go through subjects
    img$gbl <- sapply(img$sbj, function(k)      # subject k
    {
        w <- img$nwt[, k]                  # weights
        f <- img$sfs[c('slc', 'tck'), , k] # features
        v <- img$vtx                       # vertices
        
        ## shortest inter vertex geodesic distance
        d <- shortest.paths(img$nwk, weights = w)[v, v]
        
        ## gaussian blur
        sapply(sdv, function(j)         # the j.th sdv, 2^(j+1)
        {
            d <- dnorm(d, 0, j)
            diag(d) <- nrow(d) * diag(d)
            d <- d / colSums(d)
            f %*% d
        }, simplify = 'array')
    }, simplify = 'array')
    names(dimnames(img$gbl)) <- list('ftr', 'vtx', 'sdv', 'sbj')

    img
}
