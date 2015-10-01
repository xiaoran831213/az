options(warn = 2);
UTL<-new.env();

                                        # environment copy
UTL$ecp<-function(env)
{
    cpy<-new.env();
    for(x in ls(env))
        assign(x, get(x, env), cpy);
    cpy;
}

                                        # get error message from try-error object
UTL$err_msg<-function(try_error)
{
    sub('\\(converted from warning\\) ', '', attr(try_error, 'condition')[['message']])
}

                                        # clean the environmente
UTL$clr<-function()
{
    lsr<-ls(envir = parent.frame());
    rm(list=lsr, envir = parent.frame(), inherits = F);
}

UTL$binGet<-function(x, root='bin')
{
    name<-as.character(substitute(x));
    path<-sprintf('%s/%s.%s', root, name, 'bin');
    envr<-parent.frame();
    if(file.exists(path))
    {
        load(file = path, envir = envr);
        TRUE;
    }
    else
    {
        FALSE;
    }
}

UTL$binPut<-function(x, overwrite=F, root='bin')
{
    name<-as.character(substitute(x));
    path<-sprintf('%s/%s.%s', root, name, 'bin');

    if(file.exists(path))
    {
        if(overwrite)
        {
            ret <- 1L;
        }
        else
        {
            ret <- -1L;
        }
    }
    else
        ret <- 0L;
    if(ret > -1L)
        save(file = path, x);
    ret;
}

## check is an object is a scalar
.scalar <- function(obj)
{
    if(is.null(obj))
        return(TRUE)
    if(is.list(obj))
        return(FALSE)
    if(is.vector(obj) & length(obj) < 2L)
        return(TRUE)
    if(is.factor(obj) & length(obj) < 2L)
        return(TRUE)
    FALSE
}

## collect object in a function environment, by default only
## visible scalars are collected
.record <- function(pass=.scalar)
{
    ret <- list()
    env <- parent.frame()
    for(nm in ls(env))
    {
        obj <- env[[nm]]
        if(!pass(obj))
            next
        if(is.null(obj))
           obj <- 'NL'
        ret[[nm]] <- obj
    }
    ret
}

.sc1 <- function(x)
{
    (x - min(x)) / (max(x) - min(x))
}

## lower triangle of a matrix
.lwt <- function(x, ret.idx=0L)
{
    n <- nrow(x)
    z <- sequence(n)
    idx <- cbind(
        i = unlist(lapply(2L:n, seq.int, to=n), use.names = FALSE),
        j = rep.int(z[-n], rev(z)[-1L]))
        
    if(ret.idx == 1L)
        with(idx, i * n + j)
    else if(ret.idx == 2L)
        idx
    else
        x[idx]
}

.rds.rpt <- function(src, ...)
{
    ## pick out images by file name
    dirs <- c(src, ...)
    
    fns <- unlist(lapply(dirs, function(d) file.path(d, dir(d, '*.rds'))))
    dat <- lapply(fns, readRDS)
    do.call(rbind, dat)
}

