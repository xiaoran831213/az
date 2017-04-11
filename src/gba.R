#!/usr/bin/env Rscript
source('src/hlp.R')
source('src/dsg.R')
source('src/GGRF.R')
source('src/VarCompScoreTest.R')
library(SKAT)

get.phe <- function()
{
    # use baseline case and control
    d <- subset(
        read.csv('raw/ADNIMERGE.csv', as.is = T),
        VISCODE == 'bl' & PTETHCAT == 'Not Hisp/Latino' &
        PTRACCAT == 'White' & PTMARRY != 'Unknow',
        c(PTID, AGE, PTGENDER, PTEDUCAT, PTMARRY, APOE4,
          CDRSB, ADAS11, ADAS13, MMSE, RAVLT_immediate, RAVLT_learning,
          RAVLT_forgetting, RAVLT_perc_forgetting, FAQ,
          Ventricles, Hippocampus, WholeBrain, Entorhinal, Fusiform, MidTemp, ICV))
    d <- na.omit(d)
    rownames(d) <- d$PTID

    ## normalize ventricals
    d <- within(d,
    {
        VentricRNQ <- qnorm((rank(Ventricles)-0.5)/length(Ventricles))
    })
    d
}

histos <- function(p)
{
    nm <- c('Ventricles', 'Hippocampus', 'WholeBrain', 'Entorhinal', 'Fusiform', 'MidTemp', 'ICV')
    sapply(nm, function(x)
    {
        dat = p[, x]
        png(paste(x, '_hist', '.png', sep = ''))
        ## par(mar=c(2,2,2,2))
        hist(dat, breaks = 20, main=paste('Histogram of', x), xlab = NULL, ylab = NULL)
        dev.off()

        ## png(paste(x, '_log_hist', '.png', sep = ''))
        ## hist(log(dat))
        ## dev.off()
    })
    invisible(NULL)
}

## genomic group analysis (GBA)
gba <- function(gno, ...)
{
    set.seed(1234)
    
    ## load dosage data
    gno <- dosage(gno)
    if(nrow(gno) == 0L)
        return('NG')                    # null genotype
    gno <- impute(gno)

    ## phenotypes and covariants
    phe <- get.phe()[sbj(gno), ]

    sid <- intersect(rownames(phe), sbj(gno))
    phe <- phe[sid, ]
    
    ## minner allele frequency
    maf <- maf(gno)
    
    y <- phe[,  'Hippocampus']
    # y <- qnorm((rank(y)-0.5)/length(y))

    dim(y) <- c(length(y), 1)
    x <- within(
        phe[, c('AGE', 'PTGENDER', 'PTEDUCAT', 'PTMARRY')],
    {
        SEX = c(Male = 0, Female = 1)[PTGENDER]
        MAR_DIV = as.numeric(PTMARRY == 'Divorced')
        MAR_MAR = as.numeric(PTMARRY == 'Married')
        MAR_NVR = as.numeric(PTMARRY == 'Never married')
        rm(PTMARRY, PTGENDER)
    })
    x <- as.matrix(x)
    
    ## genomic variant weights
    g <- rmDgr(gno$gmx[, sid])
    g <- t(g)
    
    EQU = rep(1, ncol(g))               # unweighted
    ##:ess-bp-start::browser@nil:##
browser(expr=is.null(.ESSBP.[["@10@"]]))##:ess-bp-end:##
    BTA = dbeta(maf, 1, 25)             # beta(1,25)
    MAF = sqrt(1/(maf * (1-maf)))       # weight sum statistics
    LOG = sqrt(-log10(maf))             # logged weight

    ## SKAT kernels
    LNR = "linear.weighted"
    IBS = "IBS.weighted"

    ## SKAT methods
    OPT <- "optimal.adj"
    DVS <- "davies"

    ## covariant
    CVR <- "5CV"

    ## SKAT
    ## SN <- SKAT_Null_Model(y ~ 1 + x, out_type = 'C')
    ## rt <- rbind(
        ## c('LNR', 'EQU', 'OPT', CVR, SKAT(g, SN, LNR, OPT, weights = EQU)$p.value))
        ## c('LNR', 'BTA', 'OPT', CVR, SKAT(g, SN, LNR, OPT, weights = BTA)$p.value))
        ## c('LNR', 'MAF', 'OPT', CVR, SKAT(g, SN, LNR, OPT, weights = MAF)$p.value))
        ## c('LNR', 'LOG', 'OPT', CVR, SKAT(g, SN, LNR, OPT, weights = LOG)$p.value))
        ##c('IBS', 'EQU', 'DVS', CVR, SKAT(g, SN, IBS, DVS, weights = EQU)$p.value),
        ##c('IBS', 'BTA', 'DVS', CVR, SKAT(g, SN, IBS, DVS, weights = BTA)$p.value),
        ##c('IBS', 'MAF', 'DVS', CVR, SKAT(g, SN, IBS, DVS, weights = MAF)$p.value)),
        ##c('IBS', 'LOG', 'DVS', CVR, SKAT(g, SN, IBS, DVS, weights = LOG)$p.value))
    ## rt.SKT <- cbind('SKT', rt)

    ## GGRF
    ## rt.GRF <- c('GRF', 'IBS', 'BTA', 'DVS', CVR, GGRF(y, g, x, weights = BTA^2)$pvalue)
    
    ## VarScoreTest
    rt.CAR <- rbind(
        ## c('CAR', 'IBS', 'EQU', 'LMT', CVR, VarScoreTest(y, g, cbind(1, x), weights = EQU^2)),
        c('CAR', 'IBS', 'BTA', 'LMT', CVR, VarScoreTest(y, g, cbind(1, x), weights = BTA^2)),
        c('CAR', 'IBS', 'MAF', 'LMT', CVR, VarScoreTest(y, g, cbind(1, x), weights = MAF^2)),
        c('CAR', 'IBS', 'LOG', 'LMT', CVR, VarScoreTest(y, g, cbind(1, x), weights = LOG^2)))

    ## final report
    ## rt <- rbind(rt.SKT, rt.GRF, rt.CAR)
    rt <- rt.CAR
    
    ## compile and return
    set.seed(NULL)
    colnames(rt) <- c('alg', 'krn', 'wgt', 'mtd', 'cvr', 'pvl')
    rt <- with(gno,
    {
        cbind(
            seq=sub('.rds$', '', basename(gno$rds)), # serial number
            chr=chr, bp1=bp1, bp2=bp2,               # location
            smb=smb,                                 # symble
            wnd=wnd,                                 # flanking window
            ngv=nrow(gmx),                           # No. of variants
            rt)
    })
    
    rt <- data.frame(rt, stringsAsFactors = FALSE)
    rownames(rt) <- NULL
    rt$pvl <- as.numeric(rt$pvl)
    rt
}

## the command line parser
.argparser <- function()
{
    library(argparser, quietly = T)
    p <- arg_parser('VCF to DSG converter, RDS version.')
    aa <- add_argument

    p <- aa(p, 'gno', help = 'the input genotype .')
    p <- aa(p, 'out', help = 'the output.')
    p <- aa(p, '--ovr', help = 'overwrite existing output.', default = FALSE)
    p <- aa(p, '--vbs', help = 'verbose', default = 1L)
    p
}

main <- function(...)
{
    p <- .argparser()
    argv <- c(commandArgs(trailingOnly = TRUE), ...)
    
    if(length(argv) < 1)
    {
        print(p)
        return()
    }

    ## parse the arguments
    opt <- parse_args(p, argv)
    msg <- 'xt:'
    vbs <- opt$vbs
    
    ## initialization at command line
    opt <- within(
        opt,
    {
        msg <<- paste(msg, gno, out)
        
        ## verbosity option
        if(vbs > 1L)
            print(opt)

        ## check overwrite
        if(!ovr && file.exists(out))
            msg <<- paste(msg, 'EX')

        rm(vbs, ovr, help, opts)
    })
    ## drop NULL to facilitate function call
    opt <- opt[!sapply(opt, is.null)]

    ## skip existing file if OVR is off
    if(!grepl('EX$', msg))
    {
        rt <- do.call(gba, args = opt)
        if(is.data.frame(rt))
        {
            msg <- paste(msg, 'OK')
            saveRDS(rt, opt$out)
        }
        else
            msg <- paste(msg, rt)
    }
    else
        rt <- readRDS(opt$out)

    cat(msg); cat('\n')
    invisible(rt)
}

test <- function()
{
    a <- unlist(strsplit("dat/gs1/G4E70.rds aa.rds --ovr T", ' '))
}

main()
