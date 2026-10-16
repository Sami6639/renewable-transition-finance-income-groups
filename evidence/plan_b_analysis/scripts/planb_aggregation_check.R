root<-commandArgs(trailingOnly=TRUE)[1];ev<-file.path(root,'plan_b_analysis/audit_evidence');d<-read.csv(file.path(root,'reports/final_model_selection_20261016/four_income_group_final_dataset.csv'),check.names=FALSE)
out<-list();i<-0;G<-c('Low income','Lower middle income','Upper middle income','High income');cols<-c('FD.AST.PRVT.GD.ZS','regulatory_quality_estimate','gdp_per_capita_constant_2015_usd','population','trade_pct_gdp','inflation_cpi_pct','irena2025_pdf')
ok<-complete.cases(d[,cols]) & d$group %in% G & d$population>0 & d$gdp_per_capita_constant_2015_usd>0 & d$irena2025_pdf>=0
a<-d[which(ok),];repeat{n<-nrow(a);ct<-table(a$country_code);tt<-table(a$year);a<-a[ct[a$country_code]>=2 & tt[as.character(a$year)]>=2,];if(n==nrow(a))break}
p<-read.csv(file.path(ev,'annual_TWFE_absolute_data.csv'));key<-function(z)paste(z$country_code,z$year);i<-i+1;out[[i]]<-data.frame(check='Fresh raw-source annual eligibility',absolute_error=length(setdiff(key(a),key(p)))+length(setdiff(key(p),key(a))),tolerance=0)
for(h in 2:3) {
 w<-read.csv(file.path(ev,paste0('complete_',h,'year_windows.csv')));errs<-numeric();complete<-TRUE
 for(j in seq_len(nrow(w))) {
  r<-w[j,];z<-d[d$country_code==r$country_code & d$year>=r$t & d$year<=r$end_year,];complete<-complete && nrow(z)==h && all(complete.cases(z[,cols]))
  calc<-c(FD=mean(z[['FD.AST.PRVT.GD.ZS']]),RQ=mean(z$regulatory_quality_estimate),lnGDP=mean(log(z$gdp_per_capita_constant_2015_usd)),lnPop=mean(log(z$population)),Trade=mean(z$trade_pct_gdp),Inflation=mean(z$inflation_cpi_pct),Yabs=mean(z$irena2025_pdf),Ypc=mean(z$irena2025_pdf*1e6/z$population))
  errs<-c(errs,abs(as.numeric(r[names(calc)])-calc))
 }
 i<-i+1;out[[i]]<-data.frame(check=paste0(h,'year independent R window aggregation'),absolute_error=max(errs),tolerance=1e-9)
 i<-i+1;out[[i]]<-data.frame(check=paste0(h,'year complete observed-year requirement'),absolute_error=as.numeric(!complete),tolerance=0)
}
cat('\nBEGIN_EXPORT:R_independent_aggregation_checks.csv\n');write.csv(do.call(rbind,out),stdout(),row.names=FALSE);cat('\nEND_EXPORT\n')
