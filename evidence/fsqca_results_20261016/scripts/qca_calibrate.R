args <- commandArgs(trailingOnly=TRUE)
root <- args[1]; out <- file.path(root,'fsqca_results_20261016')
.libPaths(c(file.path(out,'tools/R_library'),.libPaths()))
library(QCA)
options(digits=17,warn=1)
emit <- function(name,x) {cat(paste0('\nBEGIN_EXPORT:',name,'\n'));write.csv(x,stdout(),row.names=FALSE,na='');cat('\nEND_EXPORT\n')}
d <- read.csv(file.path(root,'fsqca_preparation_20261016/fsqca_case_dataset.csv'),check.names=FALSE)
stopifnot(nrow(d)==111,!anyDuplicated(d$country_code))
raw <- c('bank_depth_mean_2018_2020','RQ_mean_2018_2020','real_GDP_pc_mean_2018_2020','financing_mean_per_capita_2021_2023','financing_mean_musd2022_2021_2023')
anchors <- list(c(20,50,100),c(-1,0,1),c(2000,8000,20000),c(1,10,50),c(1,25,100)); sets<-c('FD','RQ','GDP','Ypc','Yabs')
direct <- function(x,a) plogis(log(19)*(x-a[2])/ifelse(x<a[2],a[2]-a[1],a[3]-a[2]))
checks <- list()
for(j in seq_along(sets)) {
 x<-d[[raw[j]]];a<-anchors[[j]];stopifnot(!anyNA(x))
 d[[sets[j]]]<-direct(x,a)
 q<-QCA::calibrate(x,type='fuzzy',method='direct',thresholds=a,logistic=TRUE,idm=.95)
 checks[[j]]<-data.frame(set=sets[j],QCA_max_difference=max(abs(q-d[[sets[j]]])),anchor_max_difference=max(abs(direct(a,a)-c(.05,.5,.95))),exact_crossovers=sum(d[[sets[j]]]==.5))
}
d$executed_calibration_status <- 'Approved piecewise logistic; original raw case file preserved'
emit('fsqca_calibrated_dataset.csv',d);emit('R_calibration_checks.csv',do.call(rbind,checks))
cat('\nSESSION_INFO\n');print(sessionInfo())
