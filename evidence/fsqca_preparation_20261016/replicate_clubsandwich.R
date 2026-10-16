options(timeout=120,warn=1)
root <- Sys.getenv('REPLICATION_EVIDENCE_ROOT'); if (!nzchar(root)) stop('Set REPLICATION_EVIDENCE_ROOT to the absolute evidence directory')
out <- file.path(root,'fsqca_preparation_20261016');ev <- file.path(out,'audit_evidence');lib <- file.path(out,'tools/R_library')
dir.create(lib,recursive=TRUE,showWarnings=FALSE);.libPaths(c(lib,.libPaths()))
if(!requireNamespace('clubSandwich',quietly=TRUE)) install.packages('clubSandwich',lib=lib,repos='https://cloud.r-project.org',type='binary',dependencies=NA)
library(clubSandwich)
raw <- read.csv(file.path(root,'reports/final_model_selection_20261016/four_income_group_final_dataset.csv'),check.names=FALSE,stringsAsFactors=FALSE,na.strings=c('','NA'))
groups <- c('Low income','Lower middle income','Upper middle income','High income')
x <- c('FD.AST.PRVT.GD.ZS','regulatory_quality_estimate','gdp_per_capita_constant_2015_usd','population','trade_pct_gdp','inflation_cpi_pct','irena2025_pdf')
ok <- complete.cases(raw[,x]) & raw$group %in% groups & raw$population>0 & raw$gdp_per_capita_constant_2015_usd>0 & raw$irena2025_pdf>=0
dat <- raw[which(ok),]
repeat { before<-nrow(dat);ct<-table(dat$country_code);yt<-table(dat$year);dat<-dat[ct[dat$country_code]>=2 & yt[as.character(dat$year)]>=2,];if(nrow(dat)==before) break }
stopifnot(nrow(dat)==1104,length(unique(dat$country_code))==137)
write.csv(dat[,c('country_code','year','group')],file.path(ev,'R_independent_sample.csv'),row.names=FALSE)
dat$FD <- dat[['FD.AST.PRVT.GD.ZS']]-mean(dat[['FD.AST.PRVT.GD.ZS']]);dat$Governance<-dat$regulatory_quality_estimate-mean(dat$regulatory_quality_estimate)
dat$y<-log1p(dat$irena2025_pdf);dat$log_real_GDP_pc<-log(dat$gdp_per_capita_constant_2015_usd);dat$log_population<-log(dat$population);dat$trade_GDP<-dat$trade_pct_gdp;dat$inflation<-dat$inflation_cpi_pct
dat$FD_x_Governance<-dat$FD*dat$Governance
for(i in 2:4) { short<-c('LMI','UMI','HI')[i-1];dat[[paste0('group_',short)]]<-as.numeric(dat$group==groups[i]);dat[[paste0('FD_x_',short)]]<-dat$FD*dat[[paste0('group_',short)]] }
base <- 'y ~ FD + Governance + log_real_GDP_pc + log_population + trade_GDP + inflation + factor(country_code) + factor(year)'
forms <- list('Primary Baseline'=base,'Primary RQ interaction'=paste(base,'+ FD_x_Governance'),'Primary Income interaction'=paste(base,'+ group_LMI + group_UMI + group_HI + FD_x_LMI + FD_x_UMI + FD_x_HI'))
map<-c(group_LMI='group_Lower middle income',group_UMI='group_Upper middle income',group_HI='group_High income',FD_x_LMI='FD_x_Lower middle income',FD_x_UMI='FD_x_Upper middle income',FD_x_HI='FD_x_High income')
cc<-list();cis<-list();jj<-list();mm<-list()
for(label in names(forms)) {
 fit<-lm(as.formula(forms[[label]]),data=dat,x=TRUE,y=TRUE)
 V<-vcovCR(fit,cluster=dat$country_code,type='CR2',target=rep(1,nrow(dat)))
 ct<-as.data.frame(coef_test(fit,vcov=V,test='Satterthwaite'));ct$R_term<-rownames(ct);ct$term<-ifelse(ct$R_term %in% names(map),unname(map[ct$R_term]),ct$R_term);ct$model<-label;cc[[label]]<-ct
 ci<-as.data.frame(conf_int(fit,vcov=V,test='Satterthwaite',level=.95));ci$R_term<-rownames(ci);ci$term<-ifelse(ci$R_term %in% names(map),unname(map[ci$R_term]),ci$R_term);ci$model<-label;cis[[label]]<-ci
 if(label!='Primary Baseline') {
  terms<-if(label=='Primary RQ interaction') 'FD_x_Governance' else c('FD_x_LMI','FD_x_UMI','FD_x_HI')
  w<-as.data.frame(Wald_test(fit,constraints=constrain_zero(which(names(coef(fit)) %in% terms)),vcov=V,test='HTZ'));w$model<-label;jj[[label]]<-w
 }
 mm[[label]]<-data.frame(model=label,n=nobs(fit),clusters=length(unique(dat$country_code)),rank=fit$rank,residual_df=df.residual(fit),FD_center=mean(dat[['FD.AST.PRVT.GD.ZS']]),RQ_center=mean(dat$regulatory_quality_estimate))
 cat(label,'complete; rank',fit$rank,'\n')
}
write.csv(do.call(rbind,cc),file.path(ev,'R_clubSandwich_coefficients.csv'),row.names=FALSE)
write.csv(do.call(rbind,cis),file.path(ev,'R_clubSandwich_confidence_intervals.csv'),row.names=FALSE)
write.csv(do.call(rbind,jj),file.path(ev,'R_clubSandwich_joint_tests.csv'),row.names=FALSE)
write.csv(do.call(rbind,mm),file.path(ev,'R_model_metadata.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(ev,'R_session_info.txt'))
writeLines(paste('clubSandwich',as.character(packageVersion('clubSandwich'))),file.path(ev,'clubSandwich_version.txt'))
