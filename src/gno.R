source('r_src/utl.R')
require(data.table)

GNO<-new.env();

## make missing values in genotype matrix
GNO$mkMis<-function(gmx, frq, val=NA)
{
    gmx[sample(x=c(T,F), size=length(gmx), replace=T, prob=c(frq,1-frq))]<-val;
    gmx;
}

## get genotype matrix statistics, only works for biallelic dosage data.
GNO$stt <- function(gmx)
{
    if(length(gmx) == 0)
        stop("empty genotype matrix.")
    
    if(!is.matrix(gmx))                 # only one variant
        gmx = t(gmx)
    
    ngv<-nrow(gmx); # number of genome variants
    ndv<-ncol(gmx); # number of individuals
    # get value statistics for all variants
    stt<-sapply(X=1L:ngv, simplify = F, FUN = function(i)
    {
        g<-gmx[i,];
        n0<-sum(g==0L);            # 1.Homo-major
        n1<-sum(g==1L);            # 2.Hete
        n2<-sum(g==2L);            # 3.Homo-minor
        nn<-ndv-sum(n0,n1,n2);     # 4.missing
        vr<-ndv-nn-max(n0,n1,n2);  # 5.variation
        c(n0, n1, n2, nn, vr);
    });
    stt<-do.call(rbind, stt);
    colnames(stt)<-c('N0', 'N1', 'N2', 'NN', 'NV');
    stt
}

## non-physically remove degenerated variants by index.
GNO$clr<-function(gno)
{
    # get gnomic matrix
    gmx<-gno$gmx;
    
    # get index
    idx<-gno$idx;
    if(is.null(idx))
        idx <- 1L:nrow(gmx);
    
    # calculate genotype statistics is necessary
    stt<-gno$stt;
    if(is.null(stt))
        stt <- GNO$stt(gmx);
    
    # mask degenerated variants.
    msk <- which(stt[idx, 'NV']>0); # variation count
    idx <- idx[msk]
    ret<-new.env();
    for(n in ls(gno))
        assign(n, get(n, gno), ret);
    ret$stt<-stt;
    ret$idx<-idx;
    ret
}

## impute missing genotype.
GNO$imp<-function(gno)
{
    # get gnomic matrix
    gmx<-gno$gmx;
    
    # create gnomic index if necessary
    idx <- gno$idx;
    if(is.null(idx))
        idx <- 1L:nrow(gmx);
    
    # calculate genotype statistics if necessary
    stt <- gno$stt;
    if(is.null(stt))
        stt<-GNO$stt(gmx);

    msk <- stt[idx, 'NN']>0L;
    imp <- idx[msk];
    # guess missings for a variant, maintain type frequency
    for(i in imp)
    {
        v<-sample(x=0L:2L, size=stt[i,4L], replace=T, prob=stt[i,1L:3L]);
        gmx[i,gmx[i,]==3L] <- v;
        stt[i,] <- GNO$stt(gmx[i,]);
    }

    ret<-new.env();
    for(n in ls(gno))
        assign(n, get(n, gno), ret);
    ret$gmx<-gmx;
    ret$stt<-stt;
    ret;
}

#read genotype from compressed VCF file
GNO$vcf<-function(vcf, chr=NULL, bp1=NULL, bp2=NULL, idv=NULL, map=vcf)
{
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
    
    # read genome map
    cmd<-"bcftools query -f '%CHROM %POS %REF %ALT\\n'";
    if(!is.null(chr))
    {
        if(is.null(bp1))	bp1<-1L;
        if(is.null(bp2))	bp2<-0x7FFFFFFFL
        if(chr==22L)		chr<-'X';
        if(chr==23L)		chr<-'Y';
        rng<-sprintf("-r %s:%d-%d", chr, bp1, bp2);
        cmd<-paste(cmd, rng);
    }
    cmd<-paste(cmd, map);
    sed<-'sed';
    sed<-paste(sed, "-e 's/\\(^X\\)/22/'"); # X chromosome is coded 22
    sed<-paste(sed, "-e 's/\\(^Y\\)/23/'"); # Y chromosome is coded 23
    cmd<-paste(cmd, sed, sep= "|");
    pip<-pipe(cmd, "r");
    map<-try(read.table(file = pip, header = F), silent = T);
    close(pip);
    if (inherits(map, 'try-error'))
    {
        gno$err<-UTL$err_msg(try_error = map)
        return(gno);
    }
    map<-data.table(map);
    setnames(map, c("chr","pos","al1","al2"));
    setkey(map, chr, pos);
    
    # -------- get subject information -------- #
    iid<-'/tmp/iid';
    if(is.null(idv))
    {
        pip<-pipe(paste("bin/bcftools query -l", vcf), "r");
        idv<-read.table(file = pip, as.is = T);
        close(pip);
        iid<-NULL;
    }
    else if(length(idv)==1 & is.character(idv))
    {
        # write single column individual ID list file.
        cmd<-sprintf("sed -n -e '1,1 d' -e 's/^\\([^\t]*\\).*/\\1/p' <%s >%s", idv, iid);
        system(command = cmd, intern = F);
        idv<-read.table(file = idv, sep="\t", header = T, stringsAsFactors = F);
    }
    else
    {
        # write single column individual ID list file.
        write.table(x = idv[,1L], file = iid, quote = F, row.names = F, col.names = F);
    }
    
    # -------- get genotype matrix -------- #
    cmd<-"bin/bcftools query -f '[%GT ]\\n'";
    if(!is.null(rng)) #genome range exist
        cmd<-paste(cmd, rng);
    if(!is.null(iid)) #sample list file exists
        cmd<-paste(cmd, "-S", iid);
    cmd<-paste(cmd, vcf);
    
    # construct sed filter command
    sed<-'sed';
    sed<-paste(sed, "-e 's/0[|/]0/0/g'");
    sed<-paste(sed, "-e 's/[1-9][|/]0/1/g'");
    sed<-paste(sed, "-e 's/0[|/][1-9]/1/g'");
    sed<-paste(sed, "-e 's/[1-9][|/][1-9]/2/g'");
    sed<-paste(sed, "-e 's/\\.[|/]./3/g'");
    sed<-paste(sed, "-e 's/.[|/]\\./3/g'");
    
    # the final command
    cmd<-paste(cmd, sed, sep="|");
    pip<-pipe(cmd, "r");
    gmx<-matrix(scan(pip, what=0L, quiet=T), nrow = nrow(map), ncol=nrow(idv), byrow = T);
    close(pip);
    
    gno<-new.env();
    gno$gmx<-gmx;
    gno$map<-map;
    gno$idv<-idv;
    gno$chr<-chr;
    gno$bp1<-bp1;
    gno$bp2<-bp2;
    gno
}

#read genome segmentagion from file
GNO$seg<-function(whr, chr=NULL, bp1=NULL, bp2=NULL)
{
    dat<-read.table(file=whr, sep = "\t", header = T, as.is = T);
    dat<-data.table(
        dat[, c('sn', 'chr', 'bp1', 'bp2', 'id')],
        key = c("chr", "bp1", "bp2"));
    
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
    dat[c1<=chr & chr<=c2 & b1<=bp1 & bp2<=b2,];
}

## randomly pick segments of genotype variant from genome
## wnd --- segment window size
## vcf --- VCF file to read. (Varient Call Format)
## idv --- individual list to read.
## seg --- segment table to read.
## n   --- number of segments to pick
GNO$pck<-function(vcf, seg, idv, wnd=5000L, n=20L)
{
    ## list available chromosomes
    if(file.info(vcf)$isdir)
        vcf = dir(vcf)
    mt = grep('^.*chr([0-9]{1,2}).*', cs)
    ch = sub('^.*chr([0-9]{1,2}).*', '\\1', cs)[mt]
    
    ## load segment list & genotype data.
    seg<-GNO$seg(whr = seg, chr = 3);   # for now, only chr03

    ## load individual list, they are shared among all candicate segments
    idv<-read.table(file=idv, header = F, as.is = T);

    ## pick candidate segments
    sel <- sample.int(n = nrow(seg), size = n, replace = F)
    sel <- sort(sel);
    
    ## extract genome segments
    ret <- list()
    for(i in sel)
    {
        ## the range
        r<-seg[i,, drop=T];
        
        ## extract genotype
        g<-GNO$vcf(
            vcf,
            chr = r$chr, bp1 = r$bp1 - wnd, bp2 = r$bp2 + wnd,
            idv = idv);
        if(!is.null(g$err))
        {
            cat(g$err, '\n', file = stderr());
            next;
        }
        
        ## record segment index
        key <- sprintf('G%s.%s', r$sn, r$id);
        ret[[key]] <- g;
        cat(sprintf(
            '%-20s%8d%4d%12d%12d\n',
            key, nrow(g$map), g$chr, g$bp1, g$bp2));
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
