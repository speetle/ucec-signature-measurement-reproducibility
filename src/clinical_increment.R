# Fixed48 incremental prognostic assessment. Revision-stage exploratory; no claim of clinical utility.
library(survival)
args<-commandArgs(TRUE);O<-if(length(args))args[1] else 'UCEC_关键缺口追补_v13_20261020/analysis'
d<-read.csv(file.path(O,'TCGA48_input.csv'),stringsAsFactors=FALSE)
d$age10<-d$age/10;d$stage<-factor(d$stage,levels=c('I','II','III','IV'));d$grade<-factor(d$grade,levels=c('G1','G2','High'));d$histology<-factor(d$histology,levels=c('Endometrioid','Nonendometrioid'));d$subtype<-factor(d$subtype,levels=c('CN_LOW','CN_HIGH','MSI','POLE'))
stopifnot(nrow(d)==503,sum(d$event)==82,all(complete.cases(d[,c('time','event','age10','stage','grade','histology','subtype','response_score')])))
X0<-model.matrix(~age10+stage+grade+histology+subtype,d)[,-1,drop=FALSE];XL<-list(clinical_molecular=X0,plus_source48=cbind(X0,response48=d$response_score));tau<-36;pts<-seq(.05,.30,.05)
g_at<-function(sf,t,left=FALSE) vapply(t,function(a){ix<-which(if(left) sf$time<a else sf$time<=a);if(length(ix))sf$surv[max(ix)] else 1},numeric(1))
fit_model<-function(X,idx){
 mu<-colMeans(X[idx,,drop=FALSE]);ss<-apply(X[idx,,drop=FALSE],2,sd);ss[ss==0]<-1
 Xtr<-sweep(sweep(X[idx,,drop=FALSE],2,mu),2,ss,'/');tt<-d$time[idx];ee<-d$event[idx]
 f<-coxph(Surv(tt,ee)~ridge(Xtr,theta=1,scale=FALSE),x=TRUE,control=coxph.control(iter.max=100));bb<-coef(f);stopifnot(all(is.finite(bb)))
 Xall<-sweep(sweep(X,2,mu),2,ss,'/');lp<-as.vector(Xall%*%bb);h<-basehaz(f,centered=FALSE);h0<-if(any(h$time<=tau))tail(h$hazard[h$time<=tau],1) else 0
 list(lp=lp,risk=1-exp(-h0*exp(lp)),mu=mu,sd=ss,beta=bb,fit=f,G=survfit(Surv(tt,1-ee)~1))
}
measure<-function(idx,f){
 tt<-d$time[idx];ee<-d$event[idx];risk<-f$risk[idx];lp<-f$lp[idx];case<-tt<=tau & ee==1;control<-tt>tau;w<-rep(0,length(idx));raw<-rep(NA_real_,length(idx));raw[case]<-g_at(f$G,tt[case],TRUE);raw[control]<-g_at(f$G,rep(tau,sum(control)));w[case|control]<-1/pmax(raw[case|control],.05)
 cc<-concordance(Surv(tt,ee)~lp,reverse=TRUE)$concordance
 bs<-mean(w*(as.numeric(case)-risk)^2)
 slope<-tryCatch(as.numeric(coef(coxph(Surv(tt,ee)~lp))[1]),error=function(e)NA_real_)
 nb<-vapply(pts,function(p)mean(w*(risk>=p)*(as.numeric(case)-(1-as.numeric(case))*p/(1-p))),numeric(1));all<-vapply(pts,function(p)mean(w*(as.numeric(case)-(1-as.numeric(case))*p/(1-p))),numeric(1))
 list(metrics=c(C_index=cc,Brier36=bs,calibration_slope=slope),NB=nb,all=all,truncated_G=sum(raw[case|control]<.05))
}
full<-lapply(XL,fit_model,idx=seq_len(nrow(d)));app<-lapply(full,function(f)measure(seq_len(nrow(d)),f));rows<-list();dc<-list();bad<-list();set.seed(2026102007)
for(b in seq_len(500)){
 idx<-sample.int(nrow(d),replace=TRUE);oob<-setdiff(seq_len(nrow(d)),unique(idx));stopifnot(length(oob)>30)
 for(n in names(XL)){
  ans<-tryCatch({f<-fit_model(XL[[n]],idx);tr<-measure(idx,f);orig<-measure(seq_len(nrow(d)),f);te<-measure(oob,f);list(f=f,tr=tr,orig=orig,te=te)},error=function(e)e)
  if(inherits(ans,'error')){bad[[length(bad)+1]]<-data.frame(iteration=b,model=n,error=conditionMessage(ans));next}
  rows[[length(rows)+1]]<-data.frame(iteration=b,model=n,metric=names(ans$tr$metrics),training=ans$tr$metrics,original=ans$orig$metrics,oob=ans$te$metrics,optimism=ans$tr$metrics-ans$orig$metrics,oob_n=length(oob),oob_events=sum(d$event[oob]),truncated_G=ans$te$truncated_G)
  dc[[length(dc)+1]]<-data.frame(iteration=b,model=n,threshold=pts,oob_NB=ans$te$NB,oob_treat_all_NB=ans$te$all,oob_n=length(oob),oob_events=sum(d$event[oob]),truncated_G=ans$te$truncated_G)
 }
 if(b%%50==0)cat('completed',b,flush=TRUE,'\n')
}
bt<-do.call(rbind,rows);nb<-do.call(rbind,dc);write.csv(bt,file.path(O,'source48_increment_bootstrap.csv'),row.names=FALSE);write.csv(nb,file.path(O,'source48_decision_bootstrap.csv'),row.names=FALSE)
perf<-do.call(rbind,lapply(names(XL),function(n)do.call(rbind,lapply(names(app[[n]]$metrics),function(m){z<-bt[bt$model==n & bt$metric==m,];data.frame(model=n,metric=m,apparent=app[[n]]$metrics[m],corrected=app[[n]]$metrics[m]-mean(z$optimism,na.rm=TRUE),successful_draws=nrow(z))}))))
write.csv(perf,file.path(O,'source48_performance.csv'),row.names=FALSE)
pair<-merge(bt[bt$model=='plus_source48',],bt[bt$model=='clinical_molecular',],by=c('iteration','metric'),suffixes=c('_plus','_base'))
inc<-do.call(rbind,lapply(unique(pair$metric),function(m){z<-pair[pair$metric==m,];v<-z$oob_plus-z$oob_base;ff<-perf[perf$metric==m,];data.frame(metric=m,corrected_delta=ff$corrected[ff$model=='plus_source48']-ff$corrected[ff$model=='clinical_molecular'],oob_median_delta=median(v,na.rm=TRUE),stability_low=quantile(v,.025,na.rm=TRUE),stability_high=quantile(v,.975,na.rm=TRUE),paired_draws=sum(is.finite(v)))}))
write.csv(inc,file.path(O,'source48_increment_summary.csv'),row.names=FALSE)
np<-merge(nb[nb$model=='plus_source48',],nb[nb$model=='clinical_molecular',],by=c('iteration','threshold'),suffixes=c('_plus','_base'))
ns<-do.call(rbind,lapply(pts,function(p){z<-np[np$threshold==p,];v<-z$oob_NB_plus-z$oob_NB_base;data.frame(threshold=p,median_incremental_NB=median(v),stability_low=quantile(v,.025),stability_high=quantile(v,.975),positive_fraction=mean(v>0),paired_draws=length(v),interpretation='Illustrative thresholds, no specified clinical action; OOB stability range, not CI or clinical utility proof')}))
write.csv(ns,file.path(O,'source48_decision_summary.csv'),row.names=FALSE)
pr<-data.frame(patient=d$patient_id,time=d$time,event=d$event,source48=d$response_score,base_risk36=full$clinical_molecular$risk,plus48_risk36=full$plus_source48$risk);write.csv(pr,file.path(O,'source48_full_predictions.csv'),row.names=FALSE)
for(n in names(full))write.csv(data.frame(term=colnames(XL[[n]]),coefficient=full[[n]]$beta,mean=full[[n]]$mu,sd=full[[n]]$sd),file.path(O,paste0('source48_',n,'_coefficients.csv')),row.names=FALSE)
if(length(bad))write.csv(do.call(rbind,bad),file.path(O,'source48_failures.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(O,'clinical48_sessionInfo.txt'));print(perf);print(inc);print(ns);cat('DONE\n')
