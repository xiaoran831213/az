source('src/utl.R')
GNO<-new.env();
gno <- function(x)
{
    if(!is.list(x))
        x <- readRDS(as.character(x))
    structure(x, class = c('gno', 'list'))
}

maf.gno <- function(obj, level = 2L)
{
    rowMeans(obj$gmx, na.rm = T) / level
}

## make missing values in genotype matrix
GNO$mkMis<-function(gmx, frq, val=NA)
{
    gmx[sample(x=c(T,F), size=length(gmx), replace=T, prob=c(frq,1-frq))]<-val;
    gmx;
}

## get genotype counts, only works for biallelic dosage data for now.
GNO$cnt <- function(gmx)
{
    if(!is.matrix(gmx))
        stop('genotype must be a matrix')
    
    M <- nrow(gmx); # number of variants (features of a subject)
    N <- ncol(gmx); # number of subjects
    ## get value statistics for all variants
    hdr <- c('N0', 'N1', 'N2', 'NN', 'NV')
    stt <- matrix(0L, nrow = M, ncol = 5L, dimnames = list(NULL, hdr))
    for(i in 1L:M)
    {
        g <- gmx[i, ]                     # the i th. variant
        n0 <- sum(g == 0L, na.rm = T)     # 1.Homo-major
        n1 <- sum(g == 1L, na.rm = T)     # 2.Hete
        n2 <- sum(g == 2L, na.rm = T)     # 3.Homo-minor
        nn <- sum(is.na(g))               # 4.missing
        nv <- N - nn - max(n0, n1, n2)    # 5.net variation
        ## write down
        stt[i, ] <- c(n0, n1, n2, nn, nv) 
    }
    stt
}

## remove degenerated variants
GNO$clr.dgr<-function(gmx)
{
    if(!is.matrix(gmx))
        stop('genotype must be a matrix')

    ## get index
    idx <- 1L:nrow(gmx);
    
    # calculate genotype statistics is necessary
    cnt <- GNO$cnt(gmx);
    
    # mask degenerated variants.
    idx <- which(cnt[idx, 'NV']>0); # variation count
    
    return(gmx[idx, , drop = F])
}

## fix variants whoes MAF is greater than 0.5 by flipping their coding
GNO$fix.maf <- function(gmx, ret.idx = F)
{
    if(!is.matrix(gmx))
        stop('genotype must be a matrix')
    
    # calculate genotype counts
    maf <- rowMeans(gmx, na.rm = T) * 0.5
    
    # mask degenerated variants.
    idx <- which(maf > 0.5)

    if(ret.idx)
        return(idx)

    gmx[idx, ] = 2L - gmx[idx, ]
    gmx
}

## guess missing values based on the frequency of know values
GNO$imp<-function(gmx)
{
    ## calculate genotype count
    cnt<-GNO$cnt(gmx);

    idx <- which(cnt[, 'NN'] > 0L)
    # guess missings for a variant, maintain type frequency
    for(i in idx)
    {
        v<-sample(x = 0L:2L, size = cnt[i,4L], replace = T, prob = cnt[i, 1L:3L])
        gmx[i, is.na(gmx[i,])] <- v;
    }
    gmx
}

## read genotype from compressed VCF file
.read.vcf <- function(chr, bp1, bp2, vcf)
{
    ## the output
    gno<-list()
    gno$dir <- dirname(vcf); #VCF pool
    gno$vcf <- basename(vcf)
    gno$chr <- chr; #CHR is decided
    gno$bp1 <- bp1; #BP1 is decided
    gno$bp2 <- bp2; #BP2 is decided
    gno$err <- NULL
    
    ## read genome map
    cmd <- "bcftools query -f '%CHROM %POS %ID %REF %ALT\\n'";
    if(chr == 23L)
        rng <- sprintf("-r X:%d-%d", bp1, bp2)
    else if (chr == 24L)
        rng <- sprintf("-r Y:%d-%d", bp1, bp2)
    else
        rng <- sprintf("-r %d:%d-%d", chr, bp1, bp2)
    cmd <- paste(cmd, rng, vcf)

    sed<-'sed';
    sed<-paste(sed, "-e 's/\\(^X\\)/23/'"); # X chromosome is coded 23
    sed<-paste(sed, "-e 's/\\(^Y\\)/24/'"); # Y chromosome is coded 23
    cmd<-paste(cmd, sed, sep= "|");
    pip<-pipe(cmd, "r");
    map <- try(read.table(file = pip, header = F, as.is = T), silent = T)
    if(inherits(map, 'try-error'))
    {
        close(pip)
        stop('null g-map.')
    }
    close(pip);
    colnames(map) <- c("CHR", "POS", "UID", "REF", "ALT")
    rownames(map) <- sprintf('v%04X', 1L:nrow(map))
    
    ## -------- get subject id from VCF header -------- #
    cmd <- paste("bcftools query -l", vcf)
    pip <- pipe(cmd, "r");
    sbj <- scan(file = pip, what = " ", quiet = T)
    close(pip);
    sbj <- toupper(sbj)
    
    ## -------- get genotype matrix -------- #
    cmd<-"bcftools query -f '[%GT ]\\n'";
    cmd<-paste(cmd, rng);
    cmd<-paste(cmd, vcf);
    
    ## construct sed filter command
    sed<-'sed';
    sed<-paste(sed, "-e 's/0[|/]0/0/g'")
    sed<-paste(sed, "-e 's/[1-9][|/]0/1/g'")
    sed<-paste(sed, "-e 's/0[|/][1-9]/1/g'")
    sed<-paste(sed, "-e 's/[1-9][|/][1-9]/2/g'")
    sed<-paste(sed, "-e 's/\\.[|/]./3/g'")
    sed<-paste(sed, "-e 's/.[|/]\\./3/g'")
    
    ## the final command
    cmd<-paste(cmd, sed, sep="|");
    pip<-pipe(cmd, "r");
    gmx<-matrix(
        scan(pip, what=0L, na.strings = '3', quiet=T),
        nrow = nrow(map), ncol=length(sbj), byrow = T,
        dimnames = list(gvr=rownames(map), sbj=sbj))
    close(pip);
    
    gno$gmx<-gmx;
    gno$map<-map;
    gno$sbj<-sbj;
    gno
}

.read.vcf1 <- function(chr, bp1, bp2, vcf)
{
    ## the output
    gno <- list()
    gno$err <- NULL
    
    ## read genome map
    cmd <- sprintf("bcftools +dosage -r %s:%d-%d %s", chr, bp1, bp2, vcf)
    pip<-pipe(cmd, "r");

    ## parse the header
    hdr <- scan(pip, what="", nlines = 1L, quiet = T)
    sbj <- sub('^[[].*[]]', '', hdr[-(1:4)])
 
    ## read rest of the data
    dsg <- scan(pip, what="", na.strings='-1.0', quiet = T)
    close(pip);
    if(length(dsg) == 0L)
        stop('null genome.')

    ## the first 4 rows are genome metadata, onwards are subjects
    dim(dsg) <- c(length(hdr), length(dsg) / length(hdr))

    ## makeup variant name
    gvr <- sprintf('v%04X', 1L:ncol(dsg))

    ## variant metadata
    map <- data.frame(
        CHR=as.integer(dsg[1L,]),
        POS=as.integer(dsg[2L,]),
        REF=dsg[3L,],
        ALT=dsg[4L,],
        row.names=gvr,
        stringsAsFactors = FALSE)

    ## dosage data
    dsg <- matrix(
        as.integer(dsg[-(1L:4L), ]),
        nrow = length(gvr), byrow = TRUE,
        dimnames = list(gvr=gvr, sbj=sbj))
        
    within(
        gno,
    {
        dir=dirname(vcf); vcf=basename(vcf)
        chr=chr; bp1=bp1; bp2=bp2
        sbj=sbj; map=map; gmx=dsg
     })
}

## pick out subject from genotype data
sbj.gno <- function(gno, IDs)
{
    I <- match(IDs, gno$sbj)
    within(
        gno,
    {
        sbj <- sbj[I, drop = F];
        gmx <- gmx[, I, drop = F]
    })
}

## default parametennnnnrs
.hkg <- Sys.getenv('HG38_1KG')          # vcf pool
.hkg.bin <- paste(.hkg, 'bin', sep='.')
.hgn <- Sys.getenv('HG38_GEN')          # gene list
.bp0 <- 0L                              # lowest base pair
.bpN <- .Machine$integer.max            # highest base pair
.wnd <- 5000L                           # sampling window

## read genome segmentagion from file
.ls.seg <- function(seg.asc = .hgn, re.cache = F)
{
    dat <- read.table(seg.asc, sep = "\t", header = T, as.is = T)
    dat <- dat[with(dat, order(CHR, BP1, BP2)), ]
    rownames(dat) <- sprintf('G%04X', 1L:nrow(dat))
    dat
}

.ls.vcf <- function(vcf.dir = .hkg, ret.url = F, ret.chr = F)
{
    ## check and fetch vcf files
    if(!file.exists(vcf.dir))
        stop(paste(vcf.dir, 'does not exists.'))
    if(!file.info(vcf.dir)$isdir)
        stop(paste(vcf.dir, 'is not a directory.'))
    vcf.dir <- dir(vcf.dir, '.vcf.gz$', full.names = T)

    ## pick out chromosome files
    chr.ptn <- 'chr([[:digit:]]{1,2}).*vcf.gz$'
    chr.dir <- grep(chr.ptn, vcf.dir, value = T)
    
    ## sort chromosomes
    chr.mat = regexec(chr.ptn, chr.dir)
    chr.num <- lapply(regmatches(chr.dir, chr.mat), function(u)
    {
        as.integer(u[[2L]])
    })
    chr.srt <- sort.int(unlist(chr.num), index.return = T)
    chr.dir <- chr.dir[chr.srt$ix]
    chr.num <- chr.srt$x

    ## manage return values
    if(ret.chr)
    {
        if(ret.url)                     # both lists
            ret = list(chr=chr.num, url=chr.dir)
        else                            # numbers only
            ret = chr.num
    }
    else
    {
        if(ret.url)
            ret <- chr.dir              # only vcf files
        else                            # (num, vcf) tuples
            ret <- mapply(list, chr=chr.num, url=chr.dir, SIMPLIFY = F)
    }
    ret
}

str.gno <- function(gno)
{
    if(is.null(gno$ssn))
        with(gno, sprintf('G???? %2s:%-9d - %9d', chr, bp1, bp2))
    else
        with(gno, sprintf('%s %2s:%-9d - %9d', ssn, chr, bp1, bp2))
}

print.gno <- function(gno)
{
    print(str.gno(gno))
}

## segment genome
seg.gno <- function(
    vcf.dir = .hkg, seg.asc = .hgn,
    tgt.dir = paste(vcf.dir, 'bin', sep='.'),
    wnd = .wnd, ovr = FALSE)
{
    ## list chromosomes shared by vcf files and segmentation table
    vcf <- .ls.vcf(vcf.dir, ret.url=T, ret.chr=T)
    seg <- within(
        .ls.seg(seg.asc),
    {
        BP1 <- BP1 - wnd;
        BP2 <- BP2 + wnd;
    })
    
    chr <- intersect(vcf$chr, unique(seg$CHR))
    vcf.url <- vcf$url[vcf$chr %in% chr]
    vcf.chr <- vcf$chr[vcf$chr %in% chr]
    seg <- subset(seg, CHR %in% chr)
    SSN <- rownames(seg)

    dir.create(tgt.dir, showWarnings = F, recursive = T)    
    ret <- with(seg, mapply(CHR, BP1, BP2, SSN, FUN = function(chr, bp1, bp2, ssn)
    {
        fnm <- file.path(tgt.dir, paste(ssn, 'rds', sep='.'))
        if(file.exists(fnm) && !ovr)
        {
            cat(ssn, 'binery exists.\n')
            return(TRUE)
        }
        
        vcf <- vcf.url[match(chr, vcf.chr)]
        gno <- try(.read.vcf(chr, bp1, bp2, vcf), silent = T)
        if(inherits(gno, 'try-error'))
        {
            cat(ssn, geterrmessage())
            return(FALSE)
        }

        gno$ssn <- ssn
        saveRDS(gno, fnm)
        cat(ssn, gno.str(gno), '\n')
        return(TRUE)
    }))
    saveRDS(seg[ret, ], file.path(tgt.dir, '.seg.rds'))

    ## number of extracted segments
    sum(ret)
}

GNO$pic<-function(gfx, mrg=1L, pos=NULL, xlim=NULL, ylim=NULL, out=NULL, ...)
{
    pardefault <- par(no.readonly = T) # save plot settings
    if(!is.null(out))
    {
        png(out);
    }

    if(is.null(pos))
    {
        pos<-GNO$sdp(1L:dim(gfx)[mrg]);
    }

    # position limit
    if(is.null(xlim))
        xlim=c(pos[1L], pos[length(pos)]);
    # genotype value limit
    if(is.null(ylim))
        ylim=range(gfx);

    # dummy plot to lay the background
    plot(x=0, y=0, xlim=xlim, ylim=ylim);
    par(new=T);
    apply(X = gfx, MARGIN = mrg, FUN = function(gvr)
    {
        plot(x=pos, y=gvr, type="l", xlim=xlim, ylim=ylim, xaxt="n", yaxt="n", ylab="", xlab="", ...)
        par(new=T);
    });
    
    if(!is.null(out))
    {
        dev.off();
    }
    pardefault <- par(pardefault)    # restore plot settings
}
