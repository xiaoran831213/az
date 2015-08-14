source('src/utl.R')
GNO<-new.env();

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
GNO$clr.dgr<-function(gmx, ret.idx = F)
{
    stopifnot(is.matrix(gmx))
    ## get index
    idx <- 1L:nrow(gmx);
    
    # calculate genotype statistics is necessary
    cnt <- GNO$cnt(gmx);
    
    # mask degenerated variants.
    idx <- which(cnt[idx, 'NV']>0); # variation count

    if(ret.idx)
    {
        return(idx)
    }
    else
    {
        return(gmx[idx, , drop = F])
    }
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
GNO$vcf<-function(vcf, chr, bp1=0, bp2=Inf, sbj=NULL)
{
    ## the output
    gno<-list()
    gno$vcf<-vcf; #VCF file path is decided
    gno$chr<-chr; #CHR is decided
    gno$bp1<-bp1; #BP1 is decided
    gno$bp2<-bp2; #BP2 is decided
    
    ## read genome map
    cmd<-"bcftools query -f '%CHROM %POS %REF %ALT\\n'";
    rng<-sprintf("-r %s:%d-%d", chr, bp1, bp2);
    cmd<-paste(cmd, rng);
    cmd<-paste(cmd, vcf);

    sed<-'sed';
    sed<-paste(sed, "-e 's/\\(^X\\)/22/'"); # X chromosome is coded 22
    sed<-paste(sed, "-e 's/\\(^Y\\)/23/'"); # Y chromosome is coded 23
    cmd<-paste(cmd, sed, sep= "|");
    pip<-pipe(cmd, "r");
    ret <- try(map <- read.table(file = pip, header = F))
    if(inherits(ret, "try-error"))
    {
        gno$gmx = NULL
        gno$err = ret
        close(pip);
        return(gno)
    }
    close(pip);
    colnames(map) <- c("CHR", "POS", "REF", "ALT")
    
    ## -------- get subject id from VCF header -------- #
    cmd <- paste("bcftools query -l", vcf)
    pip <- pipe(cmd, "r");
    sbj.vcf <- scan(file = pip, what = " ", quiet = T)
    sbj.usr <- sbj
    close(pip);
    
    ## consider user specified subjects
    sbj.tmpfile <- NULL
    if(!is.null(sbj.usr))
    {
        sbj <- intersect(sbj.usr, sbj.vcf)
        if(length(sbj) < length(sbj.vcf))
        {
            sbj.tmpfile <- tempfile()      # write to tempfile
            write(sbj, sbj.tmpfile)
        }
    }
    else
    {
        sbj <- sbj.vcf
    }
    
    ## -------- get genotype matrix -------- #
    cmd<-"bcftools query -f '[%GT ]\\n'";
    if(!is.null(rng)) #genome range exist
        cmd<-paste(cmd, rng);
    if(!is.null(sbj.tmpfile)) #sample list file exists
        cmd<-paste(cmd, "-S", sbj.tmpfile);
    cmd<-paste(cmd, vcf);
    
    ## construct sed filter command
    sed<-'sed';
    sed<-paste(sed, "-e 's/0[|/]0/0/g'");
    sed<-paste(sed, "-e 's/[1-9][|/]0/1/g'");
    sed<-paste(sed, "-e 's/0[|/][1-9]/1/g'");
    sed<-paste(sed, "-e 's/[1-9][|/][1-9]/2/g'");
    sed<-paste(sed, "-e 's/\\.[|/]./3/g'");
    sed<-paste(sed, "-e 's/.[|/]\\./3/g'");
    
    ## the final command
    cmd<-paste(cmd, sed, sep="|");
    pip<-pipe(cmd, "r");
    gmx<-matrix(
        scan(pip, what=0L, na.strings = '3', quiet=T),
        nrow = nrow(map), ncol=length(sbj), byrow = T);
    close(pip);
    
    gno$gmx<-gmx;
    gno$map<-map;
    gno$sbj<-sbj;
    gno
}

## read genome segmentagion from file
.ls.seg <- function(seg.txt = NULL, chr = NULL, bp1 = NULL, bp2 = NULL)
{
    ## default parameters to ease typing during test
    if(is.null(seg.txt)) seg.txt <- Sys.getenv('HG38_GEN')
    if(is.null(chr)) chr <- 1L:26L
    if(is.null(bp1)) bp1 <- 0L
    if(is.null(bp2)) bp2 <- .Machine$integer.max

    ## sanity check!
    stopifnot(file.exists(seg.txt))
    
    ## try cache
    cache.dir <- file.path(tempdir(), dirname(seg.txt))
    cache.rds <- sub('[.][^.]*$', '.rds', basename(seg.txt))
    cache <- file.path(cache.dir, cache.rds)
    if(file.exists(cache))
        dat <- readRDS(cache)
    else
    {
        dat<-read.table(seg.txt, sep = "\t", header = T, as.is = T)
        dir.create(cache.dir, showWarnings = F, recursive = T)
        saveRDS(dat, cache)
    }

    ## filter and return
    subset(dat, CHR %in% chr & bp1 < BP1 & BP2 < bp2)
}

.pk.seg <- function(seg, wnd = 0L, size = 1L, replace = FALSE)
{
    idx <- sort(sample.int(nrow(seg), size, replace))
    within(seg[idx, ], {BP1 <- BP1 - wnd; BP2 <- BP2 + wnd;})
}

## randomly pick segments of genotype variant from genome
## wnd --- segment window size
## vcf --- VCF file to read. (Varient Call Format)
## sbj --- list of subjects to read
## seg --- segment table to read.
## size   --- number of segments to pick
GNO$pck<-function(
    vcf.dir = NULL, seg.txt = NULL, sbj.txt = NULL, chr = NULL,
    wnd = 5000L, size=1L, replace = FALSE)
{
    ## limit chromosome to available vcf files
    vcf <- .ls.vcf(vcf.dir, chr = chr, ret.num = T, ret.vcf = T)
    seg <- .ls.seg(seg.txt, chr = vcf$num)

    
    ret <- list()
    while(length(ret) < size)
    {
        ## pick some segments
        seg.pck <- .pk.seg(seg, wnd, size - length(ret), replace)
        
        ## choose vcf file
        vcf.pck <- vcf$vcf[match(seg.pck$CHR, vcf$num)]

        ## sample segments
        smp <- with(
            data.frame(seg.pck, VCF=vcf.pck, stringsAsFactors = F),
            mapply(GNO$vcf, VCF, CHR, BP1, BP2, SIMPLIFY = F))

        ## weed out failed sample (e.g., empty region)
        msk <- unlist(lapply(smp, function(u) is.null(u$err)))
        smp <- smp[msk]

        ## assign segment numbers to the samples
        names(smp) <- paste('G', seg.pck$SSN[msk], sep='')

        ## append sample pool
        ret <- c(ret, smp)
    }
    ret
}

.ls.vcf <- function(
    vcf.dir = NULL, chr = NULL,
    ret.vcf = FALSE, ret.num = FALSE, drop = TRUE)
{
    
    ## default parameters
    if(is.null(vcf.dir))
        vcf.dir <- Sys.getenv('HG38_1KG')
    if(is.null(chr))
        chr <- 1L:24L
    
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

    ## pick one or more chromosomes
    chr <- try(as.integer(unlist(chr)))
    if(inherits(chr, 'try-error'))
        stop(paste('chromosome must be integers.'))

    chr.idx <- na.omit(match(chr, chr.num))
    if(length(chr.idx) == 0)
        stop(paste('chr', chr, 'is unavailable in', vcf.dir))
    chr.num <- chr.num[chr.idx]
    chr.dir <- chr.dir[chr.idx]
    
    ## manage return values
    if(ret.num)
    {
        if(ret.vcf)                     # both lists
            ret = list(num=chr.num, vcf=chr.dir)
        else                            # numbers only
            ret = chr.num
    }
    else
    {
        if(ret.vcf)
            ret <- chr.dir              # only vcf files
        else                            # (num, vcf) tuples
            ret <- mapply(list, num=chr.num, vcf=chr.dir, SIMPLIFY = F)
    }

    ## try drop list structure
    if(drop & is.list(ret) & length(ret) == 1L)
         ret <- ret[[1L]]
    ret
}

#standarize genome position
GNO$sdp<-function(pos)
{
    if(length(pos) > 0L)
    {
        pos = (pos - pos[1L]) / (pos[length(pos)] - pos[1L]);
    }
    pos;
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
