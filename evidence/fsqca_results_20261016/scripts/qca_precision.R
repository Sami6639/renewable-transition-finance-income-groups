root<-commandArgs(trailingOnly=TRUE)[1];out<-file.path(root,'fsqca_results_20261016');.libPaths(c(file.path(out,'tools/R_library'),.libPaths()));library(QCA);options(digits=17)
d<-read.csv(file.path(out,'audit_evidence/scenario_memberships.csv'));b<-as.matrix(expand.grid(GDP=0:1,RQ=0:1,FD=0:1)[,c('FD','RQ','GDP')]);res<-list();i<-0
for(s in unique(d$scenario)) {z<-d[d$scenario==s,];p<-.Call('C_truthTable',as.matrix(z[,c('FD','RQ','GDP')]),z$Y,b,rep(1,3),PACKAGE='QCA');i<-i+1;res[[i]]<-data.frame(scenario=s,corner=apply(b,1,paste0,collapse=''),n=p[1,],incl=p[2,],PRI=p[3,])}
cat('\nBEGIN_EXPORT:R_truth_tables_unmasked.csv\n');write.csv(do.call(rbind,res),stdout(),row.names=FALSE);cat('\nEND_EXPORT\n')
