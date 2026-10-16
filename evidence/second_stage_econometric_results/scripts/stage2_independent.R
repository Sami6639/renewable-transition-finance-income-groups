root <- Sys.getenv('REPLICATION_EVIDENCE_ROOT'); if (!nzchar(root)) stop('Set REPLICATION_EVIDENCE_ROOT to the absolute evidence directory')
.libPaths(c(file.path(root,'fsqca_preparation_20261016/tools/R_library'),.libPaths()))
library(clubSandwich)
raw <- read.csv(file.path(root,'second_stage_investigation_20261016/processed_data/second_stage_candidate_grid.csv'),check.names=FALSE,stringsAsFactors=FALSE,na.strings=c('','NA'))
groups<-c('Low income','Lower middle income','Upper middle income','High income')
basefields<-c('FD.AST.PRVT.GD.ZS','regulatory_quality_estimate','gdp_per_capita_constant_2015_usd','population','trade_pct_gdp','inflation_cpi_pct','irena2025_pdf')
eligible<-complete.cases(raw[,basefields])&raw$income_group %in% groups & raw$population>0 & raw$gdp_per_capita_constant_2015_usd>0 & raw$irena2025_pdf>=0
core<-function(dat){repeat{before<-nrow(dat);ct<-table(dat$country_code);yt<-table(dat$year);dat<-dat[ct[dat$country_code]>=2 & yt[as.character(dat$year)]>=2,];if(nrow(dat)==before)break};dat}
benchmark<-core(raw[which(eligible),]);stopifnot(nrow(benchmark)==1104,length(unique(benchmark$country_code))==137)
coefs<-list();intervals<-list();meta<-list();members<-list();joint<-list()
for(label in c('A primary','B primary','C primary','A readiness interaction')){
 family<-substr(label,1,1);fields<-switch(family,A='L1_RISE_RE',B=c('L1_FDI_GDP','gdp_growth_pct'),C='L1_access_deficit')
 dat<-core(benchmark[complete.cases(benchmark[,fields,drop=FALSE]),]);nexpected<-switch(family,A=914,B=1094,C=1094);stopifnot(nrow(dat)==nexpected)
 dat$FD<-dat[['FD.AST.PRVT.GD.ZS']]-mean(dat[['FD.AST.PRVT.GD.ZS']]);dat$Governance<-dat$regulatory_quality_estimate-mean(dat$regulatory_quality_estimate)
 dat$y<-log1p(dat$irena2025_pdf);dat$log_real_GDP_pc<-log(dat$gdp_per_capita_constant_2015_usd);dat$log_population<-log(dat$population);dat$trade_GDP<-dat$trade_pct_gdp;dat$inflation<-dat$inflation_cpi_pct
 form<-'y ~ FD + Governance + log_real_GDP_pc + log_population + trade_GDP + inflation + factor(country_code) + factor(year)'
 if(family=='A'){dat$RISE<-dat$L1_RISE_RE-mean(dat$L1_RISE_RE);form<-paste(form,'+ RISE');target<-'RISE'}
 if(family=='B'){dat$FDI<-dat$L1_FDI_GDP;dat$growth<-dat$gdp_growth_pct;form<-paste(form,'+ FDI + growth');target<-'FDI'}
 if(family=='C'){dat$Deficit<-dat$L1_access_deficit;form<-paste(form,'+ Deficit');target<-'Deficit'}
 if(label=='A readiness interaction'){dat$RQ_x_RISE<-dat$Governance*dat$RISE;form<-paste(form,'+ RQ_x_RISE');target<-'RQ_x_RISE'}
 fit<-lm(as.formula(form),data=dat,x=TRUE,y=TRUE);V<-vcovCR(fit,cluster=dat$country_code,type='CR2',target=rep(1,nrow(dat)))
 ct<-as.data.frame(coef_test(fit,vcov=V,test='Satterthwaite'));ct$term<-rownames(ct);ct$model<-label;ct<-ct[!grepl('^factor|^\\(Intercept',ct$term),];coefs[[label]]<-ct
 ci<-as.data.frame(conf_int(fit,vcov=V,test='Satterthwaite',level=.95));ci$term<-rownames(ci);ci$model<-label;ci<-ci[ci$term %in% ct$term,];intervals[[label]]<-ci
 jt<-as.data.frame(Wald_test(fit,constraints=constrain_zero(which(names(coef(fit))==target)),vcov=V,test='HTZ'));jt$model<-label;joint[[label]]<-jt
 meta[[label]]<-data.frame(model=label,n=nobs(fit),clusters=length(unique(dat$country_code)),rank=fit$rank,residual_df=df.residual(fit),FD_center=mean(dat[['FD.AST.PRVT.GD.ZS']]),RQ_center=mean(dat$regulatory_quality_estimate));members[[label]]<-data.frame(model=label,country_code=dat$country_code,year=dat$year)
 cat('Completed independent R model: ',label,'\n',sep='')
}
emit<-function(x,name){cat('\nBEGIN_EXPORT:',name,'\n',sep='');write.csv(x,stdout(),row.names=FALSE);cat('\nEND_EXPORT\n')}
emit(do.call(rbind,coefs),'R_coefficients.csv');emit(do.call(rbind,intervals),'R_confidence_intervals.csv');emit(do.call(rbind,joint),'R_joint_tests.csv');emit(do.call(rbind,meta),'R_metadata.csv');emit(do.call(rbind,members),'R_sample_membership.csv')
cat('\nBEGIN_EXPORT:R_session_info.txt\n');print(sessionInfo());cat('\nEND_EXPORT\n')
