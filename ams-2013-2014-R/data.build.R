rm(list = ls())

##############################################################
## load data
##############################################################
source("fn.base.R")
library("ncdf4")
library("data.table")

tic()
cat("Loading csv data... ")

data.station <- read.csv(fn.in.file("station_info.csv"))
data.station$stid <- factor(data.station$stid, sort(unique(data.station$stid)))
data.station$elon <- data.station$elon+360

data.station$m1_nlat <- ceiling(data.station$nlat)
data.station$m1_elon <- ceiling(data.station$elon)
data.station$m1_elev <- NA
data.station$m1_harv_dist <- fn.haversine.km(
  data.station$nlat, data.station$elon,
  data.station$m1_nlat, data.station$m1_elon)
data.station$m1_elev_dist <- NA

data.station$m2_nlat <- ceiling(data.station$nlat)
data.station$m2_elon <- floor(data.station$elon)
data.station$m2_elev <- NA
data.station$m2_harv_dist <- fn.haversine.km(
  data.station$nlat, data.station$elon,
  data.station$m2_nlat, data.station$m2_elon)
data.station$m2_elev_dist <- NA

data.station$m3_nlat <- floor(data.station$nlat)
data.station$m3_elon <- ceiling(data.station$elon)
data.station$m3_elev <- NA
data.station$m3_harv_dist <- fn.haversine.km(
  data.station$nlat, data.station$elon,
  data.station$m3_nlat, data.station$m3_elon)
data.station$m3_elev_dist <- NA

data.station$m4_nlat <- floor(data.station$nlat)
data.station$m4_elon <- floor(data.station$elon)
data.station$m4_elev <- NA
data.station$m4_harv_dist <- fn.haversine.km(
  data.station$nlat, data.station$elon,
  data.station$m4_nlat, data.station$m4_elon)
data.station$m4_elev_dist <- NA

data.tr.out <- read.csv(fn.in.file("train.csv"))
data.test.out  <- read.csv(fn.in.file("sampleSubmission.csv"))
data.test.out[,-1] <- NA


data.gefs <- unique(data.frame(
  nlat = c(data.station$m1_nlat, data.station$m2_nlat, 
           data.station$m3_nlat, data.station$m4_nlat),
  elon = c(data.station$m1_elon, data.station$m2_elon, 
           data.station$m3_elon, data.station$m4_elon),
  elev = NA))

nc.elev = nc_open(fn.in.file("gefs_elevations.nc"), 
                      write=F, readunlim=T,  verbose=F)
gefs.lats <- 31:39
gefs.longs <- 254:269
for (idx in 1:nrow(data.gefs)) {
  idx.lon <- which(gefs.longs == data.gefs$elon[idx])
  idx.lat <- which(gefs.lats == data.gefs$nlat[idx])
  data.gefs$elev[idx] <- 
    (ncvar_get(nc.elev, "elevation_control")[idx.lon, idx.lat] +
    ncvar_get(nc.elev, "elevation_perturbation")[idx.lon, idx.lat])/2

}

nc_close(nc.elev)

data.gefs <- data.table(data.gefs, key=c("nlat", "elon"))

data.station$m1_elev <- data.gefs[J(data.station$m1_nlat,
                                    data.station$m1_elon)]$elev
data.station$m1_elev_dist <- data.station$elev - data.station$m1_elev

data.station$m2_elev <- data.gefs[J(data.station$m2_nlat,
                                    data.station$m2_elon)]$elev
data.station$m2_elev_dist <- data.station$elev - data.station$m2_elev

data.station$m3_elev <- data.gefs[J(data.station$m3_nlat,
                                    data.station$m3_elon)]$elev
data.station$m3_elev_dist <- data.station$elev - data.station$m3_elev

data.station$m4_elev <- data.gefs[J(data.station$m4_nlat,
                                    data.station$m4_elon)]$elev
data.station$m4_elev_dist <- data.station$elev - data.station$m4_elev

data.station$meso_loc <- factor(paste(
  paste(data.station$m1_nlat, data.station$m1_elon, sep=","),
  paste(data.station$m2_nlat, data.station$m2_elon, sep=","),
  paste(data.station$m3_nlat, data.station$m3_elon, sep=","),
  paste(data.station$m4_nlat, data.station$m4_elon, sep=","),
  sep = ";"))

data.station.meso <- NULL
for (m in 1:4) {
  cols.preffix <- paste0("m", m, "_") 
  cols.target <- c("nlat",  "elon", "elev", "harv_dist", "elev_dist")
  cols.station <- paste0(cols.preffix, cols.target)
  data.station.meso.cur <- data.frame(
    stid = data.station$stid,
    m = factor(m, levels=1:4))
  data.station.meso.cur[,paste0("m_", cols.target)] <- data.station[,cols.station]
  
  for (col.rm in cols.station) {
    data.station[[col.rm]] <- NULL
  }
  data.station.meso <- rbind(data.station.meso, data.station.meso.cur)
}
data.station.meso <- data.table(data.station.meso, key = c("stid", "m"))
data.station <- data.table(data.station, key = c("stid"))

data.season <- data.table(
  date = c(data.tr.out$Date, data.test.out$Date))
data.season$month <- format(
  as.POSIXct(strptime(data.season$date, "%Y%m%d")), format = "%m")
data.season$month <- factor(data.season$month)
data.season$year <- as.integer(format(
  as.POSIXct(strptime(data.season$date, "%Y%m%d")), format = "%m"))
data.season$year <- data.season$year - min(data.season$year)
setkeyv(data.season, "date")

data.tr.dates <- data.frame(
  date = rep(data.tr.out$Date, each=nrow(data.station)),
  stid = rep(data.station$stid, nrow(data.tr.out))
  )
data.test.dates <- data.frame(
  date = rep(data.test.out$Date, each=nrow(data.station)),
  stid = rep(data.station$stid, nrow(data.test.out))
  )

data.dates <- rbind(data.tr.dates, data.test.dates)
data.dates <- data.table(data.dates[, c("stid", "date")], key="stid")

fn.save.data("data.station")
fn.save.data("data.station.meso")
fn.save.data("data.tr.out")
fn.save.data("data.test.out")
fn.save.data("data.gefs")
fn.save.data("data.season")
fn.save.data("data.dates")

toc()


##############################################################
## creating distance data
##############################################################
source("fn.base.R")

fn.load.data("data.station")

tic()

data.station.dist <- data.table(data.station[,"stid", with=F])
for (stid.neigh in as.character(data.station$stid)) {
    data.station.dist[[stid.neigh]] <- Inf
}
for (r.cur in 1:nrow(data.station)) {
  stid.cur <- as.character(data.station$stid[r.cur])
  for (r.neigh in 1:nrow(data.station)) {
    stid.neigh <- as.character(data.station$stid[r.neigh])
    if (stid.cur != stid.neigh) {
      data.station.dist[[stid.neigh]][r.cur] <- 
        fn.haversine.km(data.station$nlat[r.cur], data.station$elon[r.cur],
                        data.station$nlat[r.neigh], data.station$elon[r.neigh])
    }
  }
}

data.station.dist.same.loc <- data.table(data.station.dist)
data.station.dist.diff.loc <- data.table(data.station.dist)

same.cnt <- 0
diff.cnt <- 0
for (r.cur in 1:nrow(data.station)) {
  stid.cur <- as.character(data.station$stid[r.cur])
  loc.cur <- as.character(data.station$meso_loc[r.cur])
  for (r.neigh in 1:nrow(data.station)) {
    stid.neigh <- as.character(data.station$stid[r.neigh])
    loc.neigh <- as.character(data.station$meso_loc[r.neigh])
    if (stid.cur != stid.neigh) {
      if (loc.cur == loc.neigh) {
        same.cnt <- same.cnt + 1
        data.station.dist.diff.loc[[stid.neigh]][r.cur] <- Inf
      } else {
        diff.cnt <- diff.cnt + 1
        data.station.dist.same.loc[[stid.neigh]][r.cur] <- Inf
      }
    }
  }
}
cat("Same loc", same.cnt, "diff loc", diff.cnt, "\n")

setkeyv(data.station.dist, key(data.station))
setkeyv(data.station.dist.same.loc, key(data.station))
setkeyv(data.station.dist.diff.loc, key(data.station))

fn.save.data("data.station.dist")
fn.save.data("data.station.dist.same.loc")
fn.save.data("data.station.dist.diff.loc")

toc()

##############################################################
## loading forecast data
##############################################################
source("fn.base.R")
fn.load.data("data.station")
fn.load.data("data.tr.out")
fn.load.data("data.test.out")
fn.load.data("data.gefs")

train.dir <- fn.in.file("train")
train.suf <- "_latlon_subset_19940101_20071231.nc"
test.dir <- fn.in.file("test")
test.suf <- "_latlon_subset_20080101_20121130.nc"
forecast.files <- list.files(train.dir)
forecast.files <- gsub(train.suf, "", forecast.files)

cat("Loading forecast data... ")

fn.register.wk()
forecasts.names <- foreach(fname=forecast.files,.combine=c) %dopar% {
  
  library("ncdf4")
  library("data.table")
  
  fn.init.worker(paste("forecast_",fname,sep=""))
  
  nc.tr.cur = nc_open(paste0(train.dir, "/", fname, train.suf), 
                      write=F, readunlim=T,  verbose=F)
  nc.test.cur = nc_open(paste0(test.dir, "/", fname, test.suf), 
                      write=F, readunlim=T,  verbose=F)

  gefs.name <- tolower(names(nc.tr.cur$var)[3])
  gefs.longs <- ncvar_get(nc.tr.cur, "lon") # 16
  gefs.lats <- ncvar_get(nc.tr.cur, "lat") # 9
  gefs.fhours <- ncvar_get(nc.tr.cur, "fhour") # 5
  gefs.ids <- ncvar_get(nc.tr.cur, "ens") # 11
  # nc.cur$var$Total_precipitation$varsize = 16    9    5   11 5113
  data.nc <- NULL

  for (r in 1:nrow(data.gefs)){
    gefs.lat <- data.gefs$nlat[r]
    gefs.long <- data.gefs$elon[r]
    idx.lat <- which(gefs.lats == gefs.lat)
    idx.lon <- which(gefs.longs == gefs.long)
    
    for (gefs.id in gefs.ids) {
      idx.gefs <- which(gefs.ids == gefs.id)
      data.nc.cur <- data.frame(
          date = c(data.tr.out$Date, data.test.out$Date),
          gefs = gefs.id,
          nlat = gefs.lat,
          elon = gefs.long)
      
      cat(gefs.lat, gefs.long, gefs.id, "\n")
      gefs.name.cur <- gsub("[^a-zA-Z0-9]", "_", 
                            paste0("fc_", gefs.name,"_",gefs.fhours))
      
      nc.tr.vec <- ncvar_get(nc.tr.cur)[idx.lon, idx.lat,, idx.gefs,]
      nc.test.vec <- ncvar_get(nc.test.cur)[idx.lon, idx.lat,, idx.gefs,]
      mc.mat <- matrix(c(nc.tr.vec,nc.test.vec),
                       nrow=nrow(data.nc.cur),
                       ncol=length(gefs.fhours),
                       byrow=T)
      
      
      data.nc.cur[,gefs.name.cur] <- mc.mat
      data.nc <- rbind(data.nc, data.nc.cur)
    }
  }
  
  data.nc <- data.table(data.nc, key = c("date", "gefs", "nlat", "elon"))
  cur.name <- paste0("data.forecast.", sub("\\_", ".", fname))
  assign(cur.name, data.nc)
  fn.save.data(cur.name)
  print(summary(data.nc))
  
  nc_close(nc.tr.cur)
  nc_close(nc.test.cur)

  fn.clean.worker()
  
  cur.name
}
fn.kill.wk()
fn.save.data("forecasts.names")

library("data.table")
data.forecast.all <- NULL
for (fc.name in forecasts.names) {
  fn.load.data(fc.name)
  data.forecast.cur <- get(fc.name)
  if (is.null(data.forecast.all)) {
    data.forecast.all <- data.forecast.cur
  } else {
    cols.key <- key(data.forecast.cur)
    if (!all(data.forecast.all[,cols.key, with=F] == 
               data.forecast.cur[,cols.key, with=F])) {
      stop("Error matching dimensions")
    }
    cols.cur <- setdiff(colnames(data.forecast.cur), 
                        key(data.forecast.cur))
    data.forecast.all[,cols.cur] <- data.forecast.cur[,cols.cur,with=F]
  }
}
fn.save.data("data.forecast.all")

##############################################################
## build training data
##############################################################
source("fn.base.R")
fn.load.data("data.station")
fn.load.data("data.station.meso")
fn.load.data("data.tr.out")
fn.load.data("data.test.out")
fn.load.data("data.season")
fn.load.data("data.forecast.all")

tic()
cat("Building training data... \n")

data.all.join <- data.frame(data.station.meso)
data.all.join <- fn.join.dt(data.all.join, data.station)

data.dates <- sort(c(data.tr.out$Date, data.test.out$Date))
data.all.join <- data.all.join[rep(1:nrow(data.all.join), 
                                   each = length(data.dates)),]
data.all.join$date <- data.dates
data.all.join <- data.all.join[order(data.all.join$date, 
                                     as.character(data.all.join$stid)),]

data.all.join <- fn.join.dt(data.all.join, data.season)

data.tr.out.exp <- NULL
for (stid in colnames(data.tr.out)[-1]) {
  data.tr.out.cur <- data.tr.out[, c("Date", stid)]
  setnames(data.tr.out.cur, c("Date", stid) , c("date", "power"))
  data.tr.out.cur$stid <- stid
  data.tr.out.cur <- data.tr.out.cur[,c("date", "stid", "power")]
  data.tr.out.exp <- rbind(data.tr.out.exp, data.tr.out.cur)
}
data.tr.out.exp <- data.table(data.tr.out.exp, key=c("date", "stid"))
data.all.join <- fn.join.dt(data.all.join, data.tr.out.exp)

gefs.all <- c(sort(unique(data.forecast.all$gefs)),11,12)
for (gefs.cur in gefs.all) {
  
  cat("creating gefs", gefs.cur, "data.frame \n")
  data.forecast.all.cur <- data.forecast.all
  if (gefs.cur %in% c(11:13)) {
    data.forecast.all.cur <- data.table(data.forecast.all.cur)
    if (gefs.cur != 13) {
      data.forecast.all.cur$gefs <- as.integer(gefs.cur)
    }
    cols.fc <- colnames(data.forecast.all.cur)[
      grepl("^fc_", colnames(data.forecast.all.cur))]
    for (col.fc in cols.fc) {
      col.fc.target <- sub("\\_[0-9]+$", "", col.fc)
      if (is.null(data.forecast.all.cur[[col.fc.target]])) {
        data.forecast.all.cur[[col.fc.target]] <- 0
      }
      data.forecast.all.cur[[col.fc.target]] <- 
        data.forecast.all.cur[[col.fc.target]] + 
        data.forecast.all.cur[[col.fc]]
      data.forecast.all.cur[[col.fc]] <- NULL
    }
  }
  if (gefs.cur == 11) {
    data.forecast.all.cur <- data.forecast.all.cur[
      ,lapply(sapply(.SD, median), identity),by=c("date", "gefs", "nlat", "elon")]
  } else if (gefs.cur == 12) {
    data.forecast.all.cur <- data.forecast.all.cur[
      ,lapply(sapply(.SD, max), identity),by=c("date", "gefs", "nlat", "elon")]
  }
  
  data.forecast.cur <- data.forecast.all.cur[gefs == gefs.cur,]
  data.forecast.cur$gefs <- NULL
  setnames(data.forecast.cur, c("nlat", "elon"), c("m_nlat", "m_elon"))
  setkeyv(data.forecast.cur, c("date", "m_nlat", "m_elon"))
  data.all.gefs.cur <- fn.join.dt(data.all.join, data.forecast.cur)
  cols.first <- c("date", "power", "stid")
  data.all.gefs.cur <- data.all.gefs.cur[
    ,c(cols.first, setdiff(colnames(data.all.gefs.cur), cols.first))]
  gefs.cur.name <- paste0("data.all.gefs.", gefs.cur)
  assign(gefs.cur.name, data.all.gefs.cur)
  fn.save.data(gefs.cur.name)
}

fn.save.data("data.tr.out.exp")

cat("done \n")
toc()

##############################################################
## cross validation indexes
##############################################################
source("fn.base.R")
fn.load.data("data.tr.out")

tic()
cat("Building cv... ")

data.cv.folds <- fn.cv.folds(data.tr.out$Date, 
                             K=3, 
                             type="consecutive")
fn.save.data("data.cv.folds")

data.cv.folds.ens <- fn.cv.folds(data.tr.out$Date, 
                             K=10, 
                             type="consecutive")
fn.save.data("data.cv.folds.ens")

cat("done \n")
toc()

##############################################################
## cross validation indexes
##############################################################
source("fn.base.R")
fn.load.data("data.tr.out")
fn.load.data("data.cv.folds.ens")


fn.register.wk()
data.station.dist.cor <- foreach(stid = colnames(data.tr.out)[-1], .combine=rbind, 
  .noexport=c("fn_opt_mae")) %dopar% {
  
  source("fn.base.R")
  library("Metrics")
  
  fn.init.worker(paste0("dist_cor_", stid))
  
  data.station.dist.cur.mae <- data.frame(stid=stid, 
                                          type="mae",
                                          data.tr.out[1,-1])
  data.station.dist.cur.mae[,-c(1:2)] <- Inf
  
  data.station.dist.cur.mape <- data.station.dist.cur.mae
  data.station.dist.cur.mape$type <- "mape"
  
  cur.y <- data.tr.out[,stid]
  
  
  for (stid.neigh in colnames(data.tr.out)[-1]) {
    if (stid.neigh == stid) next
    cur.pred <- rep(NA, nrow(data.tr.out))
    cur.x = data.tr.out[,stid.neigh,drop=F]
    cur.x$VINTERCEPT <- mean(cur.y)
    for (k in 1:(data.cv.folds.ens$K)) {
        cat(stid.neigh, "fold", k, "\n")
        tr.idx <- fn.cv.which(data.tr.out$Date, 
                          data.cv.folds.ens, k, "tr")
        test.idx <- fn.cv.which(data.tr.out$Date, 
                            data.cv.folds.ens, k, "test")
    
        opt.model <- fn.opt.train(
          y = cur.y[tr.idx],
          x = cur.x[tr.idx,,drop=F]
        )
        cur.pred[test.idx] <- fn.opt.predict(
          pars=opt.model$par,
          cur.x[test.idx,,drop=F])
    }
    data.station.dist.cur.mae[,stid.neigh] <- mae(cur.y, cur.pred)
    data.station.dist.cur.mape[,stid.neigh] <- 
      mean(ae(cur.y, cur.pred)/cur.y)
  }
  
  fn.clean.worker()
  
  rbind(data.station.dist.cur.mae,
        data.station.dist.cur.mape)
}
fn.kill.wk()

data.station.dist.mae <- data.station.dist.cor[
  data.station.dist.cor$type == "mae",]
data.station.dist.mae$type <- NULL
data.station.dist.mae <- data.table(data.station.dist.mae, key="stid")

data.station.dist.mape <- data.station.dist.cor[
  data.station.dist.cor$type == "mape",]
data.station.dist.mape$type <- NULL
data.station.dist.mape <- data.table(data.station.dist.mape, key="stid")

fn.save.data("data.station.dist.mae")
fn.save.data("data.station.dist.mape")

