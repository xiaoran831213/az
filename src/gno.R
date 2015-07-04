source('src/utl.R')
require(data.table)

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
GNO$vcf<-function(vcf, chr=NULL, bp1=NULL, bp2=NULL, sbj=NULL)
{
    ## the output
    gno<-new.env();
    cmd<-NULL;
    rng<-NULL;
    chr<-as.integer(chr);
    bp1<-as.integer(bp1);
    bp2<-as.integer(bp2);
    
    gno$vcf<-vcf; #VCF file path is decided
    gno$chr<-chr; #CHR is decided
    gno$bp1<-bp1; #BP1 is decided
    gno$bp2<-bp2; #BP2 is decided
    
    ## read genome map
    cmd<-"bcftools query -f '%CHROM %POS %REF %ALT\\n'";
    if(!is.null(chr))       # consider range
    {
        if(is.null(bp1))	bp1<-1L;
        if(is.null(bp2))	bp2<-0x7FFFFFFFL
        if(chr==22L)		chr<-'X';
        if(chr==23L)		chr<-'Y';
        rng<-sprintf("-r %s:%d-%d", chr, bp1, bp2);
        cmd<-paste(cmd, rng);
    }
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
    map<-data.table(map);
    setnames(map, c("chr","pos","ref","alt"));
    setkey(map, chr, pos);
    
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
    
    gno<-new.env();
    gno$gmx<-gmx;
    gno$map<-map;
    gno$sbj<-sbj;
    gno$chr<-chr;
    gno$bp1<-bp1;
    gno$bp2<-bp2;
    gno
}

## read genome segmentagion from file
GNO$seg <- function(whr, chr=NULL, bp1=NULL, bp2=NULL)
{
    dat<-read.table(file=whr, sep = "\t", header = T, as.is = T);
    dat<-data.table(
        dat[, c('SSN', 'CHR', 'BP1', 'BP2', 'GEN')],
        key = c('CHR', 'BP1', 'BP2'));
    
    if(is.null(chr))
    {
        c1=0;
        c2=0x7FFFFFFFL;
    }
    else
    {
        c1=chr;
        c2=chr;
    }
    
    if(is.null(bp1))
    {
        b1=0L;
    }
    else
    {
        b1=bp1;
    }
    
    if(is.null(bp2))
    {
        b2=0x7FFFFFFFL;
    }
    else
    {
        b2=bp2;
    }
    dat[c1<=CHR & CHR<=c2 & b1<=BP1 & BP2<=b2,];
}

## randomly pick segments of genotype variant from genome
## wnd --- segment window size
## vcf --- VCF file to read. (Varient Call Format)
## sbj --- list of subjects to read
## seg --- segment table to read.
## n   --- number of segments to pick
GNO$pck<-function(vcf, seg, sbj=NULL, wnd=5000L, n=20L)
{
    ## list available files
    if(file.info(vcf)$isdir)
    {
        vcf = dir(vcf, full.names = T)[grep('vcf.gz$', dir(vcf))]
    }
    
    ## load segment list & genotype data.
    seg <- GNO$seg(whr = seg)

    ## pick segments
    ret <- list()
    while(length(ret) < n)
    {
        ## choose file, and decide chromosome
        f = sample(vcf, size = 1)
        pt = '^.*chr([0-9]{1,2}).*'
        mt = grep(pt, f)
        if(length(mt) > 0)
            ch = as.integer(sub(pt, "\\1", f))
        else
            ch = sample.int(n=24, size = 1)

        ## decide segment
        s = seg[CHR == ch,][sample.int(n = .N, size = 1), drop = T]

        ## extract genotype
        g<-GNO$vcf(
            f,
            chr = s$CHR, bp1 = s$BP1 - wnd, bp2 = s$BP2 + wnd,
            sbj = sbj);
        if(!is.null(g$err))
        {
            cat(g$err, '\n', file = stderr());
        }
        else
        {
            ## record segment index
            key <- sprintf('G%s.%s', s$SSN, s$GEN);
            ret[[key]] <- g;
            cat(sprintf(
                '%-20s%8d%4d%12d%12d\n',
                key, nrow(g$map), g$chr, g$bp1, g$bp2));
        }
    }
    ret;
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
