root<-commandArgs(trailingOnly=TRUE)[1];out<-file.path(root,'plan_b_analysis');ev<-file.path(out,'audit_evidence')
.libPaths(c(file.path(root,'fsqca_preparation_20261016/tools/R_library'),.libPaths()));library(clubSandwich);library(nlme);options(digits=17,warn=1)
emit<-function(name,x){cat(paste0('\nBEGIN_EXPORT:',name,'\n'));write.csv(x,stdout(),row.names=FALSE,na='');cat('\nEND_EXPORT\n')}
plan<-read.csv(file.path(ev,'model_plan.csv'),stringsAsFactors=FALSE);co<-list();md<-list();jo<-list();fails<-list();covs<-list();fi<-0
for(i in seq_len(nrow(plan))) {
 m<-plan[i,];z<-read.csv(file.path(ev,paste0(m$model,'_data.csv')),check.names=FALSE);X<-read.csv(file.path(ev,paste0(m$model,'_design.csv')),check.names=FALSE);dat<-cbind(z[,c('country_code','y')],X);namesX<-names(X)
 form<-as.formula(paste('y ~ 0 +',paste(namesX,collapse=' + ')))
 fit<-tryCatch(if(m$method=='CRE')nlme::lme(form,random=~1|country_code,data=dat,method='REML',control=lmeControl(maxIter=200,msMaxIter=300,niterEM=50,opt='nlminb')) else lm(form,data=dat,x=TRUE,y=TRUE),error=function(e)e)
 if(inherits(fit,'error')) {fi<-fi+1;fails[[fi]]<-data.frame(model=m$model,procedure='fit',error=conditionMessage(fit));next}
 V<-if(m$method=='CRE') vcovCR(fit,cluster=dat$country_code,type='CR2') else vcovCR(fit,cluster=dat$country_code,type='CR2',target=rep(1,nrow(dat)))
 ct<-as.data.frame(coef_test(fit,vcov=V,test='Satterthwaite'))
 focus<-if(m$method=='TWFE')c('FD','RQ','lnGDP','lnPop','Trade','Inflation') else c(paste0('W_',c('FD','RQ','lnGDP','lnPop','Trade','Inflation')),paste0('B_',c('FD','RQ','lnGDP','lnPop','Trade','Inflation')))
 ct$term<-rownames(ct);ct<-ct[ct$term %in% focus,];co[[m$model]]<-data.frame(model=m$model,term=ct$term,estimate=ct$beta,SE_CR2=ct$SE,Satterthwaite_df=ct$df_Satt,p_value=ct$p_Satt,CI95_low=ct$beta-qt(.975,ct$df_Satt)*ct$SE,CI95_high=ct$beta+qt(.975,ct$df_Satt)*ct$SE)
 tau<-if(m$method=='CRE')as.numeric(VarCorr(fit)[1,'Variance']) else 0;sigma<-if(m$method=='CRE')fit$sigma^2 else sum(resid(fit)^2)/(nrow(dat)-ncol(X))
 md[[m$model]]<-data.frame(model=m$model,n=nrow(dat),clusters=length(unique(dat$country_code)),rank=ncol(X),residual_df=nrow(dat)-ncol(X),tau2=tau,sigma2=sigma,ICC=tau/(tau+sigma),logLik=as.numeric(logLik(fit)),CR2_target=if(m$method=='CRE')'Fitted random-intercept covariance' else 'Identity')
 mat<-as.matrix(V);focusidx<-match(focus,namesX);q<-expand.grid(term1=focus,term2=focus,stringsAsFactors=FALSE);q$covariance_CR2<-as.vector(mat[focusidx,focusidx]);q$model<-m$model;covs[[m$model]]<-q
 if(m$model=='annual_CRE_absolute') {
  for(vars in list(c('FD','RQ','lnGDP','lnPop','Trade','Inflation'),'FD','RQ')) {
   C<-matrix(0,nrow=length(vars),ncol=length(namesX));for(j in seq_along(vars)){C[j,match(paste0('W_',vars[j]),namesX)]<-1;C[j,match(paste0('B_',vars[j]),namesX)]<--1}
   w<-as.data.frame(Wald_test(fit,constraints=C,vcov=V,test='HTZ'));lab<-paste(vars,collapse=';');jo[[lab]]<-data.frame(model=m$model,test=paste0('within_equals_between_',lab),restrictions=length(vars),HTZ_F=w$Fstat,HTZ_df1=w$df_num,HTZ_df2=w$df_denom,p_value=w$p_val)
  }
 }
 cat('MODEL_DONE ',m$model,' n=',nrow(dat),' clusters=',length(unique(dat$country_code)),' tau2=',tau,'\n',sep='');flush.console()
}
emit('R_coefficients.csv',do.call(rbind,co));emit('R_model_metadata.csv',do.call(rbind,md));emit('R_joint_tests.csv',do.call(rbind,jo));emit('R_CR2_covariances.csv',do.call(rbind,covs));if(length(fails))emit('R_failures.csv',do.call(rbind,fails));cat('\nSESSION_INFO\n');print(sessionInfo())
