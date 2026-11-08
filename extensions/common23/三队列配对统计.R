args <- commandArgs(trailingOnly=TRUE); if(length(args)!=2) stop('provide input and output dirs'); p <- args[1]; out <- args[2]
labels <- c('发现队列','独立患者队列','宫颈癌队列');res<-list()
for (j in seq_along(labels)) {
 d<-read.csv(file.path(p,paste0(labels[j],' 统一子集评分.csv')),check.names=FALSE);a<-d$full48;b<-d$common;c<-d$protein;n<-nrow(d)
 delta<-function(ix) suppressWarnings(cor(b[ix],c[ix],method='spearman')-cor(a[ix],c[ix],method='spearman'))
 set.seed(2026110800+j); draws<-replicate(20000,delta(sample.int(n,n,replace=TRUE)));ci<-quantile(draws,c(.025,.975),na.rm=TRUE,names=FALSE)
 aa<-round(a,12);bb<-round(b,12);pairs<-combn(n,2);v<-(aa[pairs[1,]]-aa[pairs[2,]])*(bb[pairs[1,]]-bb[pairs[2,]]);k<-ceiling(n*.25);topa<-order(-aa,d$sample)[seq_len(k)];topb<-order(-bb,d$sample)[seq_len(k)]
 res[[j]]<-data.frame(cohort=labels[j],n=n,rho_full=cor(a,c,method='spearman'),rho_common=cor(b,c,method='spearman'),delta=delta(seq_len(n)),ci_low=ci[1],ci_high=ci[2],valid_boot=sum(is.finite(draws)),rna_rho=cor(a,b,method='spearman'),reverse=sum(v<0),untied=sum(v!=0),top_n=k,top_overlap=length(intersect(topa,topb)))
}
r<-do.call(rbind,res);write.csv(r,file.path(out,'统一子集配对统计结果.csv'),row.names=FALSE,fileEncoding='UTF-8');print(r)
