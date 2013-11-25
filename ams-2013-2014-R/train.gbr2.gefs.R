rm(list = ls())
source("fn.base.R")

fn.load.data("data.tr.out")
fn.load.data("data.cv.folds")
fn.load.data("data.tr.out.exp")

#############################################################
# training
#############################################################
grid.train <- expand.grid(k = 1:(data.cv.folds$K+1), gefs = 0:12)
model.overwrite.input <- T
model.load.if.exists <- F

fn.register.wk(1)
tmp <- foreach(k=grid.train$k, gefs = grid.train$gefs, 
                        .combine=rbind) %dopar% {
  
  source("fn.base.R")
  
  data.fold <- list()
  data.fold$name <- paste0("gbr2_gefs_", gefs, "_k_", k)
  data.fold$log.fname <- paste0(data.fold$name, ".log")
  data.fold$model.fname <- fn.out.file(paste0(data.fold$name, ".RData"))
  
  data.fold$data.name <- paste0("data.all.gefs.", gefs)
  fn.load.data(data.fold$data.name)
  data.all <- data.table(get(data.fold$data.name))
  rm(list=data.fold$data.name)
  
  cols.fc <- colnames(data.all)[grepl("fc_", colnames(data.all))]
  cols.fc <- c("m_elev_dist", "m_harv_dist", cols.fc)
  
  data.m.ord <- unique(data.all[,c("stid", "m", "m_harv_dist"),with=F])
  data.m.ord <- data.m.ord[
    order(as.character(data.m.ord$stid), data.m.ord$m_harv_dist),]
  data.m.ord <- data.m.ord[,lapply(m,identity),by="stid"]
  setkeyv(data.m.ord, "stid")
  data.m.ord <- data.frame(data.m.ord[J(data.all$stid)])
  
  fn.build.df <- function(m_ord) {
  
    cols.first <- c("date", "power", "stid", "month", "meso_loc", cols.fc)
    data.all.1 <- data.all[data.all$m == "1", cols.first, with=F]
#     data.all.1$m_ord <- factor(paste(m_ord, collapse=","))
    setnames(data.all.1, cols.fc, paste0("m1_", cols.fc))
    
    data.all.2 <- data.all[data.all$m == "2", cols.fc, with=F]
    setnames(data.all.2, cols.fc, paste0("m2_", cols.fc))
    
    data.all.3 <- data.all[data.all$m == "3", cols.fc, with=F]
    setnames(data.all.3, cols.fc, paste0("m3_", cols.fc))
    
    data.all.4 <- data.all[data.all$m == "4", cols.fc, with=F]
    setnames(data.all.4, cols.fc, paste0("m4_", cols.fc))
    
    fn.gc.wait()
    
    data.frame(cbind(data.all.1,
                      data.all.2,
                      data.all.3, 
                      data.all.4))
  }
  data.all <- fn.build.df(data.m.ord)
  
  fn.gc.wait()
  
  data.fold$transf <- fns.no.transf
  cols.exclude <- c("year")
  cols.in <- setdiff(colnames(data.all), c("date", "power", cols.exclude))
  
  data.all.series.date <- data.all$date
  data.all.series.x <- data.fold$transf$to.df(data.all[,cols.in])
  data.all.series.y <- data.fold$transf$to(data.all$power)
  rm(data.all)
  fn.gc.wait()
    
  tr.idx <- fn.cv.which(data.all.series.date, data.cv.folds, k, "tr")
  data.fold$tr.file <- fn.py.file(paste0(data.fold$name,"_tr.csv"))
  if (model.overwrite.input || !file.exists(data.fold$tr.file)) {
    write.csv(
      data.frame(power=data.all.series.y[tr.idx],
                 data.all.series.x[tr.idx,]),
      file=data.fold$tr.file,
      row.names = F
    )
  }
  
  test.idx <- fn.cv.which(data.all.series.date, data.cv.folds, k, "test")
  data.fold$test.file <- fn.py.file(paste0(data.fold$name,"_test.csv"))
  data.fold$test.pred.file <- fn.py.file(paste0(data.fold$name,"_test_pred.csv"))
  
  data.fold$test.y <- data.all.series.y[test.idx]
  data.fold$test.stid <- data.all.series.x$stid[test.idx]
  data.fold$test.dates <- data.all.series.date[test.idx]
  
  if (model.overwrite.input || !file.exists(data.fold$test.file)) {
    write.csv(
      data.frame(power=data.fold$test.y,
                 data.all.series.x[test.idx,]),
      file=data.fold$test.file,
      row.names = F
    )
  }
  
  
  rm(data.all.series.x, data.all.series.y)
  
  fn.gc.wait()

  data.fold$test.pred <- data.frame(
    date = data.fold$test.dates,
    stid = data.fold$test.stid,
    gefs = gefs)
  
  save(data.fold, file=data.fold$model.fname)
  
  NULL
}
fn.kill.wk()

fn.register.wk()
gbr2.gefs.pred.all.tmp <- foreach(k=grid.train$k, gefs = grid.train$gefs, 
                        .combine=rbind) %dopar% {
  
  source("fn.base.R")
  load(fn.out.file(paste0("gbr2_gefs_", gefs, "_k_", k, ".RData")))
  
  if (!model.load.if.exists || !file.exists(data.fold$test.pred.file)) {
    
    fn.init.worker(data.fold$name)
  
    fit.args <- paste("'{\"loss\": \"lad\", \"n_estimators\": 3000, \"learning_rate\": 0.035, ",
        "\"max_features\": 80, \"max_depth\": 7, \"random_state\": 67899, \"subsample\": 0.5, \"verbose\": 2}'")
    
    if (gefs %in% 11:12) {
      fit.args <- paste("'{\"loss\": \"lad\", \"n_estimators\": 3000, \"learning_rate\": 0.035, ",
        "\"max_depth\": 7, \"random_state\": 67899, \"subsample\": 0.5, \"verbose\": 2}'")
    }
    
    if (gefs %in% 13) {
      fit.args <- paste("'{\"loss\": \"lad\", \"n_estimators\": 2000, \"learning_rate\": 0.035, ",
        "\"max_features\": 100, \"max_depth\": 7, \"random_state\": 67899, \"subsample\": 0.5, \"verbose\": 2}'")
    }
    
    system(
      paste(
        "cd ../ams-2013-2014-py && python -u sci_learn_train.py",
        "-train_data_file", data.fold$tr.file,
        "-test_data_file", data.fold$test.file,
        "-test_pred_file", data.fold$test.pred.file,
        "-test_metric mae",
        "-target_col power", 
        "-model_type GradientBoostingRegressor",
        "-model_file ", paste0(data.fold$tr.file, ".pkl"),
        "-fit_args ", fit.args,
        " >> ", paste0("../data/log/", data.fold$log.fname)
      )
    )
    
  }
    
  data.fold$test.pred$pred <- data.fold$transf$from(
    read.csv(file=data.fold$test.pred.file)$pred)
  
  fn.print.err(data.fold$transf$from(data.fold$test.y), 
               data.fold$test.pred$pred)
  
  fn.clean.worker()
  
  data.frame(data.fold$test.pred)
}
fn.kill.wk()


fn.save.data("gbr2.gefs.pred.all.tmp")

gbr2.gefs.pred.all.tmp <- fn.sort.def(gbr2.gefs.pred.all.tmp)
fn.print.err(data.tr.out.exp, gbr2.gefs.pred.all.tmp) # 2063427

for (gefs.cur in sort(unique(gbr2.gefs.pred.all.tmp$gefs))) {
  cat("saving gefs", gefs.cur, "data.frame \n")
  data.pred.cur <- gbr2.gefs.pred.all.tmp[gbr2.gefs.pred.all.tmp$gefs == gefs.cur,]
  pred.cur.name <- paste0("data.pred.gbr2.gefs.", gefs.cur)
  assign(pred.cur.name, data.pred.cur)
  fn.print.err(data.tr.out.exp, data.pred.cur)
  fn.save.data(pred.cur.name)
}

# Load saved pieces and put them togheter
while (T) {
  gbr2.gefs.pred.all <- NULL
  gefs.cur <- 0
  while (T) {
    pred.cur.name <- paste0("data.pred.gbr2.gefs.", gefs.cur)
    if (!file.exists(fn.data.file(paste0(pred.cur.name, ".RData")))) break
    fn.load.data(pred.cur.name)
    cat("Loading", pred.cur.name, "...\n")
    pred.cur <- get(pred.cur.name)
    rm(list=pred.cur.name)
    invisible(gc())
    pred.cur <- data.table(pred.cur)[
      ,mean(pred),
      by=c("date", "stid", "gefs")]
    gbr2.gefs.pred.all <- rbind(gbr2.gefs.pred.all, pred.cur)
    gefs.cur <- gefs.cur + 1
  }
  
  gbr2.gefs.pred.all.dt <- data.table(gbr2.gefs.pred.all,
                                          key=c("date", "stid", "gefs"))
  
  if (nrow(unique(gbr2.gefs.pred.all.dt)) == nrow(gbr2.gefs.pred.all.dt)) {
    break
  }
  rm(gbr2.gefs.pred.all.dt)
  cat("Rbind bug, re-reading\n")
}
gbr2.gefs.pred.all <- data.frame(gbr2.gefs.pred.all)
gbr2.gefs.pred.all <- fn.sort.def(gbr2.gefs.pred.all)

# for (gefs.cur in sort(unique(gbr2.gefs.pred.all.V1$gefs))) {
#   fn.print.err(data.tr.out.exp, 
#                gbr2.gefs.pred.all.V1[gbr2.gefs.pred.all.V1$gefs == gefs.cur,],
#                "V1")
# }

fn.save.data("gbr2.gefs.pred.all")

# fn.load.data("gbr2.gefs.pred.all")
# gbr2.gefs.pred.all.best <- gbr2.gefs.pred.all
# fn.save.data("gbr2.gefs.pred.all.best")

#############################################################
# ensenbling all predictions together
#############################################################

if (T) {

  rm(list = ls())
  source("fn.base.R")
  
  fn.load.data("data.tr.out.exp")
  fn.load.data("gbr2.gefs.pred.all")
  fn.load.data("data.cv.folds.ens")
  
#   fn.load.data("gbr2.gefs.pred.all.best")
  
#   gbr2.gefs.pred.all$V1 <- (gbr2.gefs.pred.all$V1 +
#     gbr2.gefs.pred.all.best$V1)/2
  
  fn.register.wk()
  gbr2.gefs.pred.ens.all <- foreach(k=1, .combine=rbind, .noexport=c("fn_opt_mae")) %dopar% {
    
    source("fn.base.R")
    
    fn.init.worker(paste0("gbr2_gefs_ens"))
    
    data.ens.gbr2.all.df <- data.table(gbr2.gefs.pred.all)[
      ,lapply(V1, identity)
      ,by=c("date", "stid")]
    
    data.ens.gbr2.all.df <- fn.join.dt(data.frame(data.ens.gbr2.all.df), 
                                       data.tr.out.exp)
    
    data.ens.gbr2.all.ext.df <- data.frame(data.table(gbr2.gefs.pred.all)[
      ,list(
        VMax = max(V1),
        VMed = median(V1),
        VMean = mean(V1)
      ),
      by=c("date", "stid")])
    
    data.ens.gbr2.all.df$VMax <- data.ens.gbr2.all.ext.df$VMax
    data.ens.gbr2.all.df$VMed <- data.ens.gbr2.all.ext.df$VMed
#     data.ens.gbr2.all.df$VMean <- data.ens.gbr2.all.ext.df$VMean
    
    cols.ens <- colnames(data.ens.gbr2.all.df)[
      grepl("V[0-9A-Za-z]+", colnames(data.ens.gbr2.all.df))]
    
    data.pred.ens.gbr2 <- data.ens.gbr2.all.df[, c("date", "stid")]
    data.pred.ens.gbr2$pred <- NA
    
    lb.idx <- fn.cv.which(data.ens.gbr2.all.df$date, data.cv.folds.ens, 
                          data.cv.folds.ens$K+1, "test")
    lb.idx.mult <- 1/data.cv.folds.ens$K
    data.pred.ens.gbr2$pred[lb.idx] <- 0
    
    for (k in 1:(data.cv.folds.ens$K)) {
      cat("Fold", k, "\n")
      tr.idx <- fn.cv.which(data.ens.gbr2.all.df$date, 
                            data.cv.folds.ens, k, "tr")
      test.idx <- fn.cv.which(data.ens.gbr2.all.df$date, 
                              data.cv.folds.ens, k, "test")
      
      opt.model <- fn.opt.train(
        y = data.ens.gbr2.all.df$power[tr.idx],
        x = data.ens.gbr2.all.df[tr.idx, cols.ens]
      )
    
      data.pred.ens.gbr2$pred[test.idx] <- 
        fn.opt.predict(pars=opt.model$par, 
                       data.ens.gbr2.all.df[test.idx,cols.ens])
      
      data.pred.ens.gbr2$pred[lb.idx] <- 
        data.pred.ens.gbr2$pred[lb.idx] + 
        fn.opt.predict(pars=opt.model$par, 
                       data.ens.gbr2.all.df[lb.idx,cols.ens])*lb.idx.mult
    }
    data.pred.ens.gbr2$pred <- 
                   (data.pred.ens.gbr2$pred +  
                   data.ens.gbr2.all.ext.df$VMean)/2
    
    fn.print.err(data.tr.out.exp, 
                 data.frame(
                   data.ens.gbr2.all.df[, c("date", "stid")],
                   pred = data.ens.gbr2.all.ext.df$VMean))
                 
    
    fn.print.err(data.tr.out.exp, data.pred.ens.gbr2)
    
    fn.clean.worker()
    
    data.pred.ens.gbr2
  }
  fn.kill.wk()
  
  fn.print.err(data.tr.out.exp, gbr2.gefs.pred.ens.all)
  #     size     mae
  # 1 5501074 2012959
  fn.save.data("gbr2.gefs.pred.ens.all")
  
  #   summary(gbr2.gefs.pred.ens.all[gbr2.gefs.pred.ens.all$date >= 20080101,])
  
}

#############################################################
# smoothing using neighbooring stations
#############################################################

if (T) {

  rm(list = ls())
  source("fn.base.R")
  fn.load.data("data.tr.out.exp")
  fn.load.data("data.station.dist")
  fn.load.data("data.station.dist.diff.loc")
  fn.load.data("data.station.dist.mae")
  fn.load.data("data.station.dist.mape")
  fn.load.data("gbr2.gefs.pred.all")
  fn.load.data("gbr2.gefs.pred.ens.all")
  fn.load.data("data.cv.folds.ens")
  
  gbr2.gefs.pred.all$pred <- gbr2.gefs.pred.all$V1
  
  ens.gefs <-  -1
  gbr2.gefs.pred.ens.all$gefs <- ens.gefs
  
  gbr2.gefs.pred.all <- gbr2.gefs.pred.all[
    ,colnames(gbr2.gefs.pred.all) %in% colnames(gbr2.gefs.pred.ens.all)]
  gbr2.gefs.pred.ens.all <- gbr2.gefs.pred.ens.all[,colnames(gbr2.gefs.pred.all)]
  
  data.ens.stid.dt <- data.table(rbind(gbr2.gefs.pred.all, 
                                       gbr2.gefs.pred.ens.all),
                                 key=c("date", "stid", "gefs"))
 
  fn.register.wk()
  gbr2.gefs.pred.stid.smooth.all <- foreach(stid = unique(data.tr.out.exp$stid), .combine=rbind, 
                                           .noexport=c("fn_opt_mae")) %dopar% {
  
    source("fn.base.R")
    
    fn.init.worker(paste0("gbr2_gefs_stid_neigh_", stid))
    
    data.ens.stid.df <- fn.join.dt(
      data.frame(data.ens.stid.dt)[data.ens.stid.dt$stid == stid,], 
      data.tr.out.exp)
#     setnames(data.ens.stid.df, "pred", "VPred")
    
    # Select neartest stations
    fn.select.nearest <- function(data.station.dist, n) {
      stations.nearest <- data.frame(data.station.dist[stid])
      stations.nearest <- stations.nearest[,as.character(data.station.dist$stid)]
      stations.dist <- sort(t(stations.nearest))
      vec.dist <- t(stations.nearest)
      stations.dist <- stations.nearest[
        ,vec.dist <= stations.dist[n]
        & t(stations.nearest) < Inf]
      stations.dist <- stations.dist[,order(t(stations.dist))]
      colnames(stations.dist)
    }
    
    data.ens.stid.list <- list(
      list(dist = data.station.dist.diff.loc, n=6, df = data.ens.stid.df, 
           n.back = 0),
      list(dist = data.station.dist, n=7, df = data.ens.stid.df, 
           n.back = 0),
      list(dist = data.station.dist.mae, n=7, df = data.ens.stid.df, 
           n.back = 0),
      list(dist = data.station.dist.mape, n=7, df = data.ens.stid.df, 
           n.back = 0)
    )
    
    data.pred.ens.gbr2 <- data.ens.stid.df[, c("date", "stid", "gefs")]
    lb.idx <- fn.cv.which(data.ens.stid.df$date, data.cv.folds.ens, 
                          data.cv.folds.ens$K+1, "test")
    lb.idx.mult <- 1/data.cv.folds.ens$K
    
    for (li in 1:length(data.ens.stid.list)) {
      cur.lst <- data.ens.stid.list[[li]]
      for (stid.neigh in fn.select.nearest(cur.lst$dist, cur.lst$n)) {
        join.key <- data.table(
          date = cur.lst$df$date,
          stid = factor(stid.neigh),
          gefs = cur.lst$df$gefs
        )
        cur.lst$df[[paste0("V", stid.neigh)]] <- 
          data.ens.stid.dt[join.key]$pred
      }
      cur.lst$cols.ens <- c("pred", colnames(cur.lst$df)[
        grepl("V[0-9A-Za-z]+", colnames(cur.lst$df))])
      cur.lst$col.pred <- paste0("pred",li)
      data.ens.stid.list[[li]] <- cur.lst
      data.pred.ens.gbr2[[cur.lst$col.pred]] <- 0
#       data.pred.ens.gbr2[[cur.lst$col.pred]][lb.idx] <- 0
    }
    
    # cols.ens <- setdiff(cols.ens, c("V1", "V2", "V3"))

    for (k in 1:(data.cv.folds.ens$K)) {
      cat("Fold", k, "\n")
      tr.idx <- fn.cv.which(data.ens.stid.df$date, 
                            data.cv.folds.ens, k, "tr")
      test.idx <- fn.cv.which(data.ens.stid.df$date, 
                              data.cv.folds.ens, k, "test")
      
      for (li in 1:length(data.ens.stid.list)) {
        cur.lst <- data.ens.stid.list[[li]]
        
        all.n.back <- 0:cur.lst$n.back
        n.back.mult <- 1/length(all.n.back)
        
        for (n.back in all.n.back) {
          cur.cols <- head(cur.lst$cols.ens, n=length(cur.lst$cols.ens)-n.back)
#           cur.cols <- cur.lst$cols.ens
          opt.model <- fn.opt.train(
            y = cur.lst$df$power[tr.idx],
            x = cur.lst$df[tr.idx, cur.cols]
          )
          
          data.pred.ens.gbr2[[cur.lst$col.pred]][test.idx] <- 
            data.pred.ens.gbr2[[cur.lst$col.pred]][test.idx] + 
            fn.opt.predict(pars=opt.model$par, 
                           cur.lst$df[test.idx, cur.cols])*n.back.mult
          
          data.pred.ens.gbr2[[cur.lst$col.pred]][lb.idx] <- 
            data.pred.ens.gbr2[[cur.lst$col.pred]][lb.idx] + 
            fn.opt.predict(pars=opt.model$par, 
                           cur.lst$df[lb.idx, cur.cols])*
            lb.idx.mult*n.back.mult
        }
      }
    }
    
    data.pred.ens.gbr2 <- data.pred.ens.gbr2[
      data.pred.ens.gbr2$gefs == ens.gefs,]
    data.pred.ens.gbr2$gefs <- NULL
    
    fn.print.err(data.tr.out.exp, 
                 data.ens.stid.df[data.ens.stid.df$gefs == ens.gefs,])
    
    cols.pred <- colnames(data.pred.ens.gbr2)[
      grepl("pred[0-9]+", colnames(data.pred.ens.gbr2))]
    for (col.pred in cols.pred) {
      fn.print.err(data.tr.out.exp, data.pred.ens.gbr2, col.pred)
    }
    
    fn.clean.worker()
    
    data.pred.ens.gbr2
  }
  fn.kill.wk()
  

  gbr2.gefs.pred.stid.smooth.all <- fn.sort.def(gbr2.gefs.pred.stid.smooth.all)
  cols.pred <- colnames(gbr2.gefs.pred.stid.smooth.all)[
    grepl("pred[0-9]+", colnames(gbr2.gefs.pred.stid.smooth.all))]
  for (col.pred in cols.pred) {
    fn.print.err(data.tr.out.exp, gbr2.gefs.pred.stid.smooth.all, col.pred)
  }
  #     size     mae
  # 1 501074 1989999
  #     size     mae
  # 1 501074 1986884
  #     size     mae
  # 1 501074 1987069
  #     size     mae
  # 1 501074 1992468

  #   summary(gbr2.gefs.pred.stid.smooth.all[gbr2.gefs.pred.stid.smooth.all$date >= 20080101,])
  fn.save.data("gbr2.gefs.pred.stid.smooth.all")
  
}
  
#############################################################
# final ensenble
#############################################################

if (T) {

  rm(list = ls())
  source("fn.base.R")
  
  fn.load.data("data.tr.out.exp")
  fn.load.data("gbr2.gefs.pred.ens.all")
  fn.load.data("gbr2.gefs.pred.stid.smooth.all")
  fn.load.data("data.cv.folds.ens")
  
  fn.register.wk()
  gbr2.gefs.pred.ens.final <- foreach(k=1, .combine=rbind, .noexport=c("fn_opt_mae")) %dopar% {
    
    source("fn.base.R")
    
    fn.init.worker(paste0("gbr2_gefs_ens_final"))
    
    data.ens.gbr2.final.df <- gbr2.gefs.pred.stid.smooth.all
    data.ens.gbr2.final.df$pred0 <- gbr2.gefs.pred.ens.all$pred
    data.ens.gbr2.final.df <- fn.join.dt(data.ens.gbr2.final.df, data.tr.out.exp)
    
    cols.ens.m <- colnames(data.ens.gbr2.final.df)[
      grepl("pred[0-9A-Za-z]+", colnames(data.ens.gbr2.final.df))]
    
    data.ens.gbr2.final.df$predMed <- apply(data.ens.gbr2.final.df[,cols.ens.m], 
                                           1, median)
    data.ens.gbr2.final.df$predMax <- apply(data.ens.gbr2.final.df[,cols.ens.m], 
                                           1, max)
    data.ens.gbr2.final.df$predMean <- rowMeans(
      data.ens.gbr2.final.df[,cols.ens.m])
    
    cols.ens <- colnames(data.ens.gbr2.final.df)[
      grepl("pred[0-9A-Za-z]+", colnames(data.ens.gbr2.final.df))]
    
    cols.ens <- cols.ens[!cols.ens %in% c("predMean", "predMax")]
    
    # cols.ens <- setdiff(cols.ens, c("V1", "V2", "V3"))
    data.pred.ens.gbr2 <- data.ens.gbr2.final.df[, c("date", "stid")]
    data.pred.ens.gbr2$pred2 <- NA    
    
    
    all.k <- (data.cv.folds.ens$K+1)
    
    lb.idx <- fn.cv.which(data.ens.gbr2.final.df$date, 
                          data.cv.folds.ens, 
                          data.cv.folds.ens$K+1, "test")
    lb.idx.mult <- 1/length(all.k)
    
    for (k in all.k) {
      cat("Fold", k, "\n")
      tr.idx <- fn.cv.which(data.ens.gbr2.final.df$date, 
                            data.cv.folds.ens, k, "tr")
      test.idx <- fn.cv.which(data.ens.gbr2.final.df$date, 
                              data.cv.folds.ens, k, "test")
      test.idx <- tr.idx
      
      opt.model <- fn.opt.train(
        y = data.ens.gbr2.final.df$power[tr.idx],
        x = data.ens.gbr2.final.df[tr.idx, cols.ens]
      )
    
      data.pred.ens.gbr2$pred2[test.idx] <- 
        fn.opt.predict(pars=opt.model$par, 
                       data.ens.gbr2.final.df[test.idx,cols.ens])
      
      data.pred.ens.gbr2$pred2[lb.idx] <- 
        data.pred.ens.gbr2$pred2[lb.idx] + 
        fn.opt.predict(pars=opt.model$par, 
                       data.ens.gbr2.final.df[lb.idx,cols.ens])*lb.idx.mult
    }
    
    data.pred.ens.gbr2$pred <- 
      (data.ens.gbr2.final.df$predMean + 
         data.ens.gbr2.final.df$pred2)/2
    
    fn.print.err(data.tr.out.exp, data.pred.ens.gbr2)
    
    fn.print.err(data.tr.out.exp, data.pred.ens.gbr2, "pred2")
    
    data.pred.ens.gbr2[, c("date", "stid", "pred")]
  }
  fn.kill.wk()
  
  fn.print.err(data.tr.out.exp, gbr2.gefs.pred.ens.final)
  #     size     mae
  # 1 501074 1984926
  fn.save.data("gbr2.gefs.pred.ens.final")
  
  summary(gbr2.gefs.pred.ens.final[gbr2.gefs.pred.ens.final$date < 20080101,])
  summary(gbr2.gefs.pred.ens.final[gbr2.gefs.pred.ens.final$date >= 20080101,])
  
  gbr2.gefs.pred.ens.final.sub <- gbr2.gefs.pred.ens.final
  gbr2.gefs.pred.ens.final.sub$pred <- gbr2.gefs.pred.ens.final.sub$pred*1.02
  fn.write.submission(gbr2.gefs.pred.ens.final.sub, "gbr2.gefs.pred.ens.final")
  
#   fn.load.data("gbr2.gefs.pred.ens.final")
#   gbr2.gefs.pred.ens.final.best <- gbr2.gefs.pred.ens.final
#   fn.save.data("gbr2.gefs.pred.ens.final.best")
}