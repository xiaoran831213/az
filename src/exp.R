main <- function()
{
    ## read subject id, probe id, loci name
    sid <- readLines('raw/exp/0/sid.txt')
    prb <- readLines('raw/exp/0/prb.txt')
    loc <- readLines('raw/exp/0/loc.txt')
    smb <- readLines('raw/exp/0/smb.txt')
    vtm <- readLines('raw/exp/0/vtm.txt')

    ## value matrix
    val.prb <- matrix(
        scan('raw/exp/0/val.txt'),
        length(prb), length(sid), T, list(prb=prb, sid=sid))

    ## probe map
    map <- mapply(function(p, l, s)
    {
        sp <- unlist(strsplit(s, ','))
        if(length(sp) > 0L)
            s <- sp

        ## expand probe of multiple symbols to multiple rows
        cbind(PRB = p, LOC = l, SMB = s)
    }, prb, loc, smb, USE.NAMES = F)
    map <- do.call(rbind, map)
    map <- as.data.frame(map, stringsAsFactors = F)

    ## express matrix mapped to gene symbol
    val.smb <- val.prb[map$PRB,]
    dimnames(val.smb) <- list(smb=map$SMB, sid=sid)
    
    invisible(list(
        val.prb=val.prb,
        val.smb=val.smb,
        vtm=vtm,
        map=map))
}
