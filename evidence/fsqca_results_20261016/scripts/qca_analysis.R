args<-commandArgs(trailingOnly=TRUE);root<-args[1];out<-file.path(root,'fsqca_results_20261016');ev<-file.path(out,'audit_evidence')
.libPaths(c(file.path(out,'tools/R_library'),.libPaths()));library(QCA);options(digits=17,warn=1)
emit<-function(name,x){cat(paste0('\nBEGIN_EXPORT:',name,'\n'));write.csv(x,stdout(),row.names=FALSE,na='');cat('\nEND_EXPORT\n')}
all<-read.csv(file.path(ev,'scenario_memberships.csv'),check.names=FALSE);meta<-read.csv(file.path(ev,'scenario_metadata.csv'),check.names=FALSE)
tts<-list();necess<-list();sols<-list();paths<-list();checks<-list();errors<-list();ii<-jj<-kk<-ll<-mm<-ee<-0
formula_text<-function(p) if(length(p)) paste(p,collapse=' + ') else ''
for(si in seq_len(nrow(meta))) {
 m<-meta[si,];d<-all[all$scenario==m$scenario,];rownames(d)<-d$country_code
 dd<-d[,c('FD','RQ','GDP','Y')];rownames(dd)<-rownames(d)
 # Fix both inclusion cutpoints to the same value: only supported 1/0, never C.
 tt<-truthTable(dd,outcome='Y',conditions=c('FD','RQ','GDP'),incl.cut=c(m$consistency,m$consistency),n.cut=m$frequency,pri.cut=m$PRI,complete=TRUE,show.cases=TRUE,sort.by='')
 t<-as.data.frame(tt$tt,stringsAsFactors=FALSE);t$corner<-apply(t[,c('FD','RQ','GDP')],1,paste0,collapse='');t$scenario<-m$scenario
 ii<-ii+1;tts[[ii]]<-t[,c('scenario','corner','FD','RQ','GDP','OUT','n','incl','PRI','cases')]
 for(c in c('FD','RQ','GDP'))for(neg in c(FALSE,TRUE)) {
  x<-if(neg)1-dd[[c]] else dd[[c]];p<-pof(setms=data.frame(X=x),outcome=dd$Y,relation='necessity')$incl.cov
  jj<-jj+1;necess[[jj]]<-data.frame(scenario=m$scenario,condition=paste0(if(neg)'~' else '',c),consistency=p[1,'inclN'],coverage=p[1,'covN'],RoN=p[1,'RoN'])
 }
 if(!any(t$OUT=='1')) {kk<-kk+1;sols[[kk]]<-data.frame(scenario=m$scenario,type='none',model=0,expression='',consistency=NA,PRI=NA,coverage=NA,status='No configuration meets approved criteria');next}
 for(typ in c('complex','parsimonious')) {
  z<-tryCatch(minimize(tt,include=if(typ=='parsimonious')'?' else '',details=TRUE,show.cases=TRUE,row.dom=FALSE,all.sol=TRUE),error=function(e)e)
  if(inherits(z,'error')) {ee<-ee+1;errors[[ee]]<-data.frame(scenario=m$scenario,type=typ,error=conditionMessage(z));next}
  if(m$scenario %in% c('primary_Ypc','primary_Yabs')) {cat('\nPRIMARY_PACKAGE_MINIMIZATION ',m$scenario,typ,'\n');print(z);str(z,max.level=2)}
  for(mi in seq_along(z$solution)) {
   terms<-z$solution[[mi]]
   mem<-sapply(terms,function(term) {
    if(term=='1')return(rep(1,nrow(dd)))
    lits<-strsplit(term,'*',fixed=TRUE)[[1]]
    Reduce(pmin,lapply(lits,function(lit)if(substr(lit,1,1)=='~')1-dd[[substring(lit,2)]] else dd[[lit]]))
   });if(is.null(dim(mem)))mem<-matrix(mem,ncol=1);colnames(mem)<-terms
   fit<-pof(setms=as.data.frame(mem),outcome=dd$Y,relation='sufficiency')$incl.cov
   s<-apply(mem,1,max);un<-pof(setms=data.frame(solution=s),outcome=dd$Y,relation='sufficiency')$incl.cov
   kk<-kk+1;sols[[kk]]<-data.frame(scenario=m$scenario,type=typ,model=mi,expression=formula_text(terms),consistency=un[1,'inclS'],PRI=un[1,'PRI'],coverage=un[1,'covS'],status='Boolean solution; evaluate union fit separately')
   for(pi in seq_along(terms)) {ll<-ll+1;paths[[ll]]<-data.frame(scenario=m$scenario,type=typ,model=mi,pathway=terms[pi],consistency=fit[pi,'inclS'],PRI=fit[pi,'PRI'],raw_coverage=fit[pi,'covS'],unique_coverage=fit[pi,'covU'],countries=paste(d$country_code[mem[,pi]>.5],collapse=';'))}
  }
 }
}
emit('R_truth_tables.csv',do.call(rbind,tts));emit('R_necessity.csv',do.call(rbind,necess));emit('R_solutions.csv',do.call(rbind,sols));emit('R_pathways.csv',do.call(rbind,paths));if(length(errors))emit('R_errors.csv',do.call(rbind,errors))
cat('\nSESSION_INFO\n');print(sessionInfo())
