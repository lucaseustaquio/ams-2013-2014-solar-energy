library("data.table")
library("compiler")
library("Rcpp")


enableJIT(3) 
setCompilerOptions(suppressUndefined = T)
options(stringsAsFactors = FALSE)
options(max.print = 1000)

path.wd <- getwd()

sourceCpp("fn.base.cpp")
#############################################################
# tic toc
#############################################################
tic <- function(gcFirst = TRUE, type=c("elapsed", "user.self", "sys.self")) {
  type <- match.arg(type)
  assign(".type", type, envir=baseenv())
  if(gcFirst) gc(FALSE)
  tic <- proc.time()[type]         
  assign(".tic", tic, envir=baseenv())
  invisible(tic)
}

toc <- function() {
  type <- get(".type", envir=baseenv())
  toc <- proc.time()[type]
  tic <- get(".tic", envir=baseenv())
  print(toc - tic)
  invisible(toc)
}

##############################################################
## Registers parallel workers
##############################################################
fn.register.wk <- function(n.proc = NULL) {
  if (file.exists("data/cluster.csv")) {
    cluster.conf <- read.csv("data/cluster.csv", 
                             stringsAsFactors = F,
                             comment.char = "#")
    n.proc <- NULL
    for (i in 1:nrow(cluster.conf)) {
      n.proc <- c(n.proc, 
                  rep(cluster.conf$host[i], 
                      cluster.conf$cores[i]))
    }
  }
  if (is.null(n.proc)) {
    n.proc = as.integer(Sys.getenv("NUMBER_OF_PROCESSORS"))
    if (is.na(n.proc)) {
      library(parallel)
      n.proc <-detectCores()
    }
  }
  workers <- mget(".pworkers", envir=baseenv(), ifnotfound=list(NULL));
  if (!exists(".pworkers", envir=baseenv()) || length(workers$.pworkers) == 0) {
    
    library(doSNOW)
    library(foreach)
    workers<-suppressWarnings(makeSOCKcluster(n.proc));
    suppressWarnings(registerDoSNOW(workers))
    clusterSetupRNG(workers, seed=5478557)
    assign(".pworkers", workers, envir=baseenv());
    
    tic()
    cat("Workers start time: ", format(Sys.time(), 
                                       format = "%Y-%m-%d %H:%M:%S"), "\n")
  }
  invisible(workers);
}

##############################################################
## Kill parallel workers
##############################################################
fn.kill.wk <- function() {
  library("doSNOW")
  library("foreach")
  workers <- mget(".pworkers", envir=baseenv(), ifnotfound=list(NULL));
  if (exists(".pworkers", envir=baseenv()) && length(workers$.pworkers) != 0) {
    stopCluster(workers$.pworkers);
    assign(".pworkers", NULL, envir=baseenv());
    cat("Workers finish time: ", format(Sys.time(), 
                                        format = "%Y-%m-%d %H:%M:%S"), "\n")
    toc()
  }
  invisible(workers);
}

##############################################################
## init worker setting work dir and doing path redirect
##############################################################
fn.init.worker <- function(log = NULL, add.date = F) {
  setwd(path.wd)
  
  if (!is.null(log)) {
    date.str <- format(Sys.time(), format = "%Y-%m-%d_%H-%M-%S")
    
    if (add.date) {
      output.file <- fn.log.file(paste(log, "_",date.str,
                                       ".log", sep=""))
    } else {
      output.file <- fn.log.file(paste(log,".log", sep=""))
    }
    output.file <- file(output.file, open = "wt")
    sink(output.file)
    sink(output.file, type = "message")
    
    cat("Start:", date.str, "\n")
  }
  
  tic()
}

##############################################################
## clean worker resources
##############################################################
fn.clean.worker <- function() {
  gc()
  
  try(toc(), silent=T)
  suppressWarnings(sink())
  suppressWarnings(sink(type = "message"))
}

##############################################################
## wait clean
##############################################################
fn.gc.wait <- function() {
  invisible(gc())
  Sys.sleep(1)
  invisible(gc())
  Sys.sleep(1)
  invisible(gc())
}

#############################################################
# log file path
#############################################################
fn.base.dir <- function(extra) {
  paste0(path.wd, "/../data/", extra)
}

#############################################################
# log file path
#############################################################
fn.log.file <- function(name) {
  fn.base.dir(paste0("log/", name))
}

#############################################################
# input file path
#############################################################
fn.in.file <- function(name) {
  fn.base.dir(paste0("input/", name))
}

#############################################################
# r output file path
#############################################################
fn.out.file <- function(name) {
  fn.base.dir(paste0("output-R/", name))
}

#############################################################
# python output file path
#############################################################
fn.py.file <- function(name) {
  fn.base.dir(paste0("output-py/", name))
}

#############################################################
# submission file path
#############################################################
fn.submission.file <- function(name) {
  fn.base.dir(paste0("submission/", name, ".csv"))
}

#############################################################
# data file path
#############################################################
fn.data.file <- function(name) {
  fn.out.file(name)
}

#############################################################
# save data file
#############################################################
fn.save.data <- function(dt.name, envir = parent.frame()) {
  save(list = dt.name, 
       file = fn.data.file(paste0(dt.name, ".RData")), envir = envir)
}

#############################################################
# load saved file
#############################################################
fn.load.data <- function(dt.name, envir = parent.frame()) {
  load(fn.data.file(paste0(dt.name, ".RData")), envir = envir)
}

#############################################################
# error evaluation
#############################################################
fn.print.err <- function(actual, pred, pred.col = "pred", do.print = T) { 
  library("data.table")
  if (is.data.table(actual)) {
    actual <- fn.join.dt(pred, actual)
    actual <- actual$power 
  }
  
  if (is.data.frame(actual)) {
    if (!is.null(actual$power)) {
      actual <- actual$power
    }  
  }
  
  if (is.data.frame(pred)) {
    pred <- pred[[pred.col]]
  }
  
  pred.na <- is.na(pred) | is.na(actual)
  if (any(pred.na)) {
    pred <- pred[!pred.na]
    actual <- actual[!pred.na]
  }
  
  if (is.null(actual) || all(is.na(actual))) {
    df <- summary(pred)
  } else {
    library("Metrics")
    df <- data.frame(size = length(pred), mae = mae(actual,pred))
  }
  
  
  if (do.print) {
    print(df)
  }
  
  invisible(df)
}
# debug(fn.print.err)


#############################################################
# join data tables - overwrite existing cols
#############################################################
fn.join.dt <- function(df1, df2) {
  cols.key <- key(df2)
  if (is.data.table(df1)) {
    data.key <- df1[,cols.key,with=F]
    df1 <- data.table(df1, key = key(df1))
  } else {
    data.key <- data.table(df1[,cols.key,drop=F])
    df1 <- data.frame(df1)
  }
  df2.join <- df2[data.key,allow.cartesian=TRUE]
  df2.join <- df2.join[
    ,!(colnames(df2.join) %in% cols.key),with=F]
  cols.overwrite <- colnames(df2.join)[
    colnames(df2.join) %in% colnames(df1)]
  for (col.name in cols.overwrite) {
    df1[[col.name]] <- NULL
  }
  cbind(df1, df2.join)
}



##################################################
# write libsvm data
##################################################
fn.write.libsvm <- function(
  data.tr, 
  data.test, 
  name,
  fn.y.transf = NULL, 
  dir = "../data/output",
  col.y = "stars",
  col.x = colnames(data.tr)[!(colnames(data.tr) %in% c(col.y, col.qid))],
  col.qid = NULL,
  feat.start = 1,
  vw.mode = F,
  data.val = NULL,
  y.def = min(data.tr[[col.y]]))
{
  options(scipen=999)
  col.x <- unique(col.x)
  library("data.table")
  cat("Building feature map ...")
  tic()
  model.dts <- list()
  val.xi <- feat.start
  col.x.groups <- NULL
  feat.size <- 0
  for (i in (1:length(col.x))) {
    col <- col.x[i]
      if (is.character(data.tr[[col]]) | is.factor(data.tr[[col]])) {
      
      col.ids.all   <- unique(as.character(data.tr[[col]]))
      if (!is.null(data.val)) {
        col.ids.all <- unique(c(col.ids.all, as.character(data.val[[col]])))
      }
      if (!is.null(data.test)) {
        col.ids.all <- unique(c(col.ids.all, as.character(data.test[[col]])))
      }
      col.ids.all <- col.ids.all[!is.na(col.ids.all)]
      model.dts[[col]] <- data.table(ID = col.ids.all, key = "ID")
      feat.size <- feat.size + nrow(model.dts[[col]])
      
      model.dts[[col]]$X.Idx <- val.xi:(val.xi+nrow(model.dts[[col]])-1)
      
      val.xi <- val.xi+nrow(model.dts[[col]])
      
    } else {
      model.dts[[col]] <- val.xi
      val.xi <- val.xi + 1
    }
  }
  map.name <- paste(name, ".map", sep="")
  assign(map.name, model.dts)
  save(list = map.name, file = paste(dir, "/", map.name, ".RData", sep=""))
  cat("done \n")
  toc()
  
  if (!exists("cmpfun")) {
    cmpfun <- identity
  }
  write.file <- cmpfun(function (data, file) { 
    tic()
    col.chunk <- col.x
    if (!is.null(col.qid)) {
      col.chunk <- c(col.chunk, col.qid)
    }
    if (!is.null(data[[col.y]])) {
      col.chunk <- c(col.y, col.chunk)
    }
    unlink(file)
    cat("Saving ", file, "...")
    fileConn <- file(file, open="at")
    
    data.chunk <- data[, col.chunk]
    if (is.null(data.chunk[[col.y]])) {
      data.chunk[[col.y]] <- y.def
    }
    data.chunk[[col.y]][is.na(data.chunk[[col.y]])] <- y.def
    if (!is.null(fn.y.transf)) {
      data.chunk[[col.y]] <- fn.y.transf(data.chunk[[col.y]])
    }
    data.chunk[[col.y]][data.chunk[[col.y]] == Inf] <- y.def
    
    for (col in col.x) {
      if (is.numeric(data.chunk[[col]])) {
        data.chunk[[col]][data.chunk[[col]] == 0] <- NA
      }
    }
    
    for (col in col.x) {
#       print(col)
      col.is.na <- is.na(data.chunk[[col]])
      if (is.factor(data.chunk[[col]]) || is.character(data.chunk[[col]])) {
        data.chunk[[col]] <-  paste(
          model.dts[[col]][J(as.character(data.chunk[[col]]))]$X.Idx,
          c(1), sep = ":")
      } else {
        data.chunk[[col]] <- paste(
          rep(model.dts[[col]], nrow(data.chunk)),
          data.chunk[[col]], sep = ":")
      }
      data.chunk[col.is.na,col] <- ""
    }
    
    if (!is.null(col.qid)) {
      data.chunk[[col.qid]] <- paste("qid", data.chunk[[col.qid]], sep = ":")
    }
    
    data.chunk <- do.call(paste, data.chunk[, c(col.y, col.qid, col.x)])
    chunk.size <- as.numeric(object.size(data.chunk))
    chunk.size.ch <- T
    while (chunk.size.ch) {
      data.chunk <- gsub(" [0-9]+\\:?NA", "", data.chunk)
      data.chunk <- gsub(" NA\\:-?[0-9]+", "", data.chunk)
      chunk.size.ch <- chunk.size != as.numeric(object.size(data.chunk))
      chunk.size <- as.numeric(object.size(data.chunk))
    }
    data.chunk <- gsub("\\s+", " ", data.chunk)
    data.chunk <- gsub("^([0-9]+(\\.[0-9]+)?)\\s*$", 
                       paste0("\\1 ", val.xi, ":1"), data.chunk)
    
    if (vw.mode) {
      data.chunk <- gsub("^([-]?[0-9]+(\\.[0-9]+)?)\\s+", "\\1 | ", data.chunk)
    }
    
    writeLines(c(data.chunk), fileConn)
    
    close(fileConn)
    cat("done.\n")
    toc()
  })
  #     debug(write.file)
  write.file(data.tr, paste(dir, "/", name, ".tr.libsvm", sep=""))
  if (!is.null(data.val)) {
    write.file(data.val, paste(dir, "/", name, ".val.libsvm", sep=""))
  }
  if (!is.null(data.test)) {
#     debug(write.file)
    write.file(data.test, paste(dir, "/", name, ".test.libsvm", sep=""))
  }
}
# debug(fn.write.libsvm)

##################################################
# write libfm relation files
##################################################
fn.write.libfm.relation <- function(tr.ids,
                                    test.ids,
                                    out.file,
                                    data.file,
                                    id.file = paste0(data.file, ".ids"),
                                    dir = "../data/libfm/") {
  
  tr.ids <- as.character(tr.ids)
  test.ids <- as.character(test.ids)
  
  id.file <- paste0(dir, id.file)
  out.file <- paste0(dir, out.file)
  data.file <- paste0(dir, data.file)
  
  library("data.table")
  dt.ids <- data.table(id=readLines(id.file), key="id")
  dt.ids$idx <- 1:nrow(dt.ids)
  data.rows <- readLines(data.file)
  writeLines(data.rows, out.file)
  
  writeLines(as.character(dt.ids[J(tr.ids)]$idx-1), paste0(out.file, ".train"))
  writeLines(as.character(dt.ids[J(test.ids)]$idx-1), paste0(out.file, ".test"))
  invisible(NULL)
}

##################################################
# expand factors
##################################################
fn.expand.factors <- function(data.df) {
  data.frame(model.matrix( ~ . - 1, data=data.df))
}

##############################################################
## get with default value
##############################################################
fn.get <- function(var.name, def.val = NULL) {
  if ( exists(var.name)) {
    return(get(var.name))
  } else {
    return(def.val)
  }
}


#############################################################
# extract prediction
#############################################################
fn.extract.pred <- function(data.all, data.type) {
  
  library("data.table")
  
  data.extracted <- data.table(data.all[data.all$datatype == data.type, ])
  data.extracted$datatype <- NULL
  if (is.null(data.extracted$modeltype)) {
    data.extracted <- data.frame(data.extracted[
      ,list(pred = mean(pred)), by="test.idx"])
  } else {
    data.extracted <- data.frame(data.extracted[
      ,list(
        pred1 = mean(pred[modeltype==1]),
        pred2 = mean(pred[modeltype==2]),
        pred3 = mean(pred[modeltype==3]),
        pred4 = mean(pred[modeltype==4])
      ), by="test.idx"])
  }
  data.extracted <- data.frame(data.extracted)
  
  data.extracted[data.extracted$test.idx,] <- data.extracted
  rownames(data.extracted) <- 1:nrow(data.extracted)
  
  cols.out <- colnames(data.extracted)[!(
    colnames(data.extracted) %in% c("test.idx", "datatype"))]
  
  data.extracted[, cols.out, drop = F]
}

# debug(fn.extract.pred)

fn.extract.tr <- function(data.all, ...) {
  fn.extract.pred(data.all, "tr", ...)
}

fn.extract.test <- function(data.all, ...) {
  fn.extract.pred(data.all, "test", ...)
}

#############################################################
# write submission
#############################################################
fn.write.submission <- function(pred, file.name) { 
  
  library("data.table")
  
  pred <- pred[pred$date >= 20080101,]
  pred <- fn.mean.pred(pred)
  
  to.list <- function(stid, pred) {
    list <- lapply(pred, identity)
    names(list) <- stid
    list
  }
  pred <- pred[order(pred$date, as.character(pred$stid)),]
  pred <- data.table(pred)[
    ,to.list(stid=stid, pred=pred), by="date"]
  setnames(pred, "date", "Date")
  
#   print(summary(pred))
  pred <- data.frame(pred)
  write.csv(pred,
            file = fn.submission.file(paste0(file.name)), 
            row.names = F, quote = F)
  invisible(pred)
}

#############################################################
# read submission
#############################################################
fn.read.submission <- function (sub.file) {
  data.ens <- read.csv(fn.submission.file(sub.file))
  data.ens <- data.ens[order(data.ens$Date),
                       c("Date", sort(colnames(data.ens)[-1]))]
  data.ens
}

##############################################################
## cross val folds
##############################################################
fn.cv.folds <- function(dates, K, seed = NULL, type = "consecutive") {
  n <- length(dates)
  library("cvTools")
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  data.cv.folds <- cvFolds(n, K = K, type = type)
  if (!is.null(seed)) {
    set.seed(Sys.time())
  }
  
  data.cv.folds <- list(
    n = data.cv.folds$n,
    K = data.cv.folds$K,
    which = data.cv.folds$which,
    subsets = data.cv.folds$subsets,
    dates = sort(dates))
  data.cv.folds
}

##############################################################
## cross val selection
##############################################################
fn.cv.which <- function(data.dates, cv.data, k, type) {
  if (type == "tr") {
    tr.dates <- cv.data$dates[cv.data$subsets[!cv.data$which %in% k]]
    return (which(data.dates %in% tr.dates))
  } else if (type == "test") {
    if (k <= cv.data$K) {
      test.dates <- cv.data$dates[cv.data$subsets[cv.data$which %in% k]]
      return (which(data.dates %in% test.dates))
    } else {
      return (which(!data.dates %in% cv.data$dates))
    }
  }
}
# debug(fn.write.submission)

##############################################################
## dor xor in prediction
##############################################################
fn.find.ntile <- function(val, ntile) {
  val <- sort(val)
  val[round(length(val)*ntile)]
}

#############################################################
# print rf importance
#############################################################
fn.rf.print.imp <- function(rf) {
  imp <- try(print(data.frame(
    importance=rf$importance[order(-rf$importance[,1]),]), silent = T))
}

#############################################################
# print rf importance
#############################################################
fn.haversine.km <- function(lat1,long1,lat2,long2) {
  K_PI <- 3.141592653589793;
  d2r <- (K_PI / 180.0)
  dlong <- (long2 - long1) * d2r;
  dlat  <- (lat2 - lat1) * d2r;
  a     <- sin(dlat/2.0)^2 + cos(lat1*d2r) * cos(lat2*d2r) * sin(dlong/2.0)^2;
  c     <- 2 * atan2(sqrt(a), sqrt(1-a));
  d     <- 6371.19645 * c; # mean of mean of earth radius in km
  d
}

#############################################################
# data.transformations
#############################################################
fn.to.transf.df <- function(transf.fun, df, 
                            cols=colnames(df)[grepl("^fc_", colnames(df))]) {
  df[,cols] <- transf.fun(df[,cols])
  return (df)
}

#############################################################
# no tranformation
#############################################################
fns.no.transf <- list(
  to = identity,
  from = identity,
  to.df = identity)

#############################################################
# sqrt tranformation
#############################################################
fn.to.sqrt <- function(val) {
  sign(val)*sqrt(abs(val))
}

fn.from.sqrt <- function(val) {
  sign(val)*(val^2)
}

fn.to.sqrt.df <- function(...) {
  fn.to.transf.df(transf.fun=fn.to.sqrt, ...)
}

fns.sqrt.transf <- list(
  to = fn.to.sqrt,
  from = fn.from.sqrt,
  to.df = fn.to.sqrt.df)

#############################################################
# log tranformation
#############################################################
fn.to.log <- function(val) {
  sign(val)*log(1+abs(val))
}

fn.from.log <- function(val) {
  sign(val)*(exp(abs(val))-1)
}

fn.to.log.df <- function(...) {
  fn.to.transf.df(transf.fun=fn.to.log, ...)
}

fns.log.transf <- list(
  to = fn.to.log,
  from = fn.from.log,
  to.df = fn.to.log.df)

#############################################################
# simple average of predictions
#############################################################
fn.mean.pred <- function(data.pred) {
  library("data.table")
  df <- data.frame(data.table(data.pred)[
    ,list(pred= mean(pred)),
    by=c("date", "stid")])
  df[order(df$date),]
}

#############################################################
# default sorting
#############################################################
fn.sort.def <- function(data) {
  data[order(data$date, as.character(data$stid)),]
}


##############################################################
## train using optmin
##############################################################
fn.opt.predict <- cmpfun(function(pars, data.x) {
  pars.m <- matrix(rep(pars,each=nrow(data.x)),nrow=nrow(data.x))
  rowSums(data.x*pars.m)
})

fn.opt.train <- function(x, y, pars = rep(1/ncol(x),ncol(x)), 
                         use.cpp = T, do.trace = F,
                         ...) {
  
  tic()
  
  not.na <- !is.na(y)
  x <- x[not.na,]
  y <- y[not.na]
  
  
  if (use.cpp) {
    x.mat <- as.matrix(x)
    fn.opt <- cmpfun(function(pars) {
      fn_opt_mae(x.mat, y, pars)
    })
  } else {
    library("Metrics")
    fn.opt <- cmpfun(function(pars) {
      mae(y, fn.opt.predict(pars, x))
    })
  }

  opt.model <- optim(pars, fn.opt, 
                     control = list(trace = do.trace),
                     ...)
  opt.model$summary <- data.frame(m=NA,weight=NA)
  col.x <- colnames(x)
  for (i in 1:ncol(x)) {
    opt.model$summary[i,] <- c(col.x[i],  opt.model$par[i])
  }
  print(opt.model$summary)
  toc()
  opt.model
}

##############################################################
## train using optmin
##############################################################
fn.sub.to.vector <- function(data.sub) {
  vec <- NULL
  for (col in sort(colnames(data.sub)[-1])) {
    vec <- c(vec, data.sub[[col]])
  }
  vec
}
