library(survival)
args<-commandArgs(TRUE);root<-args[1];out<-file.path(root,'results/clinical');dir.create(out,recursive=TRUE,showWarnings=FALSE)
d<-read.csv(file.path(root,'input/clinical/TCGA48_input.csv'),stringsAsFactors=FALSE)
d$age10<-d$age/10;d$stage<-factor(d$stage,levels=c('I','II','III','IV'));d$grade<-factor(d$grade,levels=c('G1','G2','High'));d$histology<-factor(d$histology);d$subtype<-factor(d$subtype);d$z<-as.numeric(scale(d$response_score))
f<-coxph(Surv(time,event)~age10+stage+grade+histology+strata(subtype)+z,data=d,x=TRUE);s<-summary(f);a<-data.frame(cohort='TCGA',model='age, stage, grade, histology adjusted; subtype-stratified Cox',score='source48',n=nrow(d),events=sum(d$event),HR=s$conf.int['z','exp(coef)'],low=s$conf.int['z','lower .95'],high=s$conf.int['z','upper .95'],P=s$coefficients['z','Pr(>|z|)'],PH_global=cox.zph(f)$table['GLOBAL','p']);write.csv(a,file.path(out,'TCGA_source48_association.csv'),row.names=FALSE);write.csv(as.data.frame(cox.zph(f)$table),file.path(out,'TCGA_source48_PH.csv'))

sens<-list()
for(model in c('subtype_stratified','stage_subtype_stratified','restricted_36m')){
 dd<-d
 if(model=='restricted_36m'){dd$event<-as.integer(dd$event==1 & dd$time<=36);dd$time<-pmin(dd$time,36)}
 ff<-if(model=='stage_subtype_stratified')Surv(time,event)~age10+grade+histology+strata(stage,subtype)+z else Surv(time,event)~age10+stage+grade+histology+strata(subtype)+z
 fit<-coxph(ff,data=dd,x=TRUE);ss<-summary(fit);ph<-cox.zph(fit)
 sens[[model]]<-data.frame(model=model,n=nrow(dd),events=sum(dd$event),HR=ss$conf.int['z','exp(coef)'],low=ss$conf.int['z','lower .95'],high=ss$conf.int['z','upper .95'],P=ss$coefficients['z','Pr(>|z|)'],score_PH_P=ph$table['z','p'],global_PH_P=ph$table['GLOBAL','p'])
 write.csv(as.data.frame(ph$table),file.path(out,paste0('PH_',model,'.csv')))
}
write.csv(do.call(rbind,sens),file.path(out,'TCGA_source48_PH_sensitivity.csv'),row.names=FALSE)

x<-read.csv(file.path(root,'input/clinical/external_survival_endpoint_audit.csv'),check.names=FALSE);x<-x[tolower(x$included)=='true',];res<-list()
for(sc in c('full48','common29')){
 x$score<-as.numeric(scale(x[[sc]]));fit<-coxph(Surv(time,event)~score,data=x,x=TRUE);s<-summary(fit);ph<-cox.zph(fit);res[[sc]]<-data.frame(score=sc,n=nrow(x),events=sum(x$event),HR=s$conf.int[1,'exp(coef)'],low=s$conf.int[1,'lower .95'],high=s$conf.int[1,'upper .95'],P=s$coefficients[1,'Pr(>|z|)'],PH_P=ph$table[1,'p'])
}
a<-do.call(rbind,res);a$BH2<-p.adjust(a$P,'BH');write.csv(a,file.path(out,'external_survival_exploratory.csv'),row.names=FALSE);print(a)
writeLines(capture.output(sessionInfo()),file.path(out,'associations_sessionInfo.txt'))
