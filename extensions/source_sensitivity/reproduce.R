args<-commandArgs(trailingOnly=TRUE)
base<-normalizePath(if(length(args)>0) args[1] else '.')
oldroot<-base
out<-file.path(base,'results');dir.create(out,recursive=TRUE,showWarnings=FALSE)
suppressPackageStartupMessages({library(edgeR);library(limma)})
d<-read.delim(gzfile(file.path(oldroot,'input/GSE235046_Count_table.txt.gz')),check.names=FALSE)
meta<-read.csv(file.path(oldroot,'assets/sample_map.csv'),check.names=FALSE)
meta<-meta[match(names(d)[5:ncol(d)],meta$count_column),];stopifnot(!anyNA(meta$gsm),all(meta$count_column==names(d)[5:ncol(d)]))
x<-as.matrix(d[,5:ncol(d)]);storage.mode(x)<-'double';stopifnot(all(is.finite(x)),all(x>=0))
pc<-d$bioType=='protein_coding' & !is.na(d$geneSymbol) & nzchar(d$geneSymbol)
cts<-rowsum(x[pc,,drop=FALSE],d$geneSymbol[pc],reorder=FALSE)
grp<-factor(meta$treatment,levels=c('Media','LPS','CS','LPS+CS','LPS+CS+Torin'));design<-model.matrix(~0+grp);colnames(design)<-c('Media','LPS','CS','Combined','Torin')
y0<-DGEList(cts);keep<-filterByExpr(y0,design=design);cts<-cts[keep,,drop=FALSE]
contr<-cbind(Combined_vs_Media=c(-1,0,0,1,0),Interaction=c(1,-1,-1,1,0),Torin_vs_Combined=c(0,0,0,-1,1));rownames(contr)<-colnames(design)
edge<-function(x){ y<-calcNormFactors(DGEList(x)); y<-estimateDisp(y,design,robust=TRUE);fit<-glmQLFit(y,design,robust=TRUE); lapply(seq_len(ncol(contr)),function(i)topTags(glmQLFTest(fit,contrast=contr[,i]),n=Inf,sort.by='none')$table) }
a<-edge(cts);b<-edge(round(cts))
y<-calcNormFactors(DGEList(cts));v<-voom(y,design,plot=FALSE);z<-eBayes(contrasts.fit(lmFit(v,design),contr),robust=TRUE)
c<-lapply(seq_len(ncol(contr)),function(i)topTable(z,coef=i,number=Inf,sort.by='none'))
legacy<-read.csv(file.path(oldroot,'assets/legacy101.csv'),check.names=FALSE)
tr<-read.csv(file.path(base,'assets/source101_gene_traceability.csv'),check.names=FALSE,fileEncoding='UTF-8-BOM')
results<-data.frame(gene=legacy$geneSymbol,direction=legacy$方向,stringsAsFactors=FALSE)
results$fixed48<-tolower(as.character(tr$inherited_fixed48[match(results$gene,tr$mouse_symbol)]))=="true"
stopifnot(!anyNA(results$fixed48),sum(results$fixed48)==48)
for(meth in c('fractional_edgeR','rounded_edgeR','fractional_voom')){
 tabs<-switch(meth,fractional_edgeR=a,rounded_edgeR=b,fractional_voom=c)
 for(i in 1:3){ix<-match(results$gene,rownames(tabs[[i]]));results[[paste0(meth,'_',colnames(contr)[i],'_logFC')]]<-tabs[[i]]$logFC[ix];results[[paste0(meth,'_',colnames(contr)[i],'_FDR')]]<-tabs[[i]][[if(meth=='fractional_voom') 'adj.P.Val' else 'FDR']][ix]}
 f<-function(name)results[[paste0(meth,'_',name)]]
 flag<-!is.na(f('Combined_vs_Media_FDR')) & f('Combined_vs_Media_FDR')<0.05 & abs(f('Combined_vs_Media_logFC'))>=1 & f('Interaction_FDR')<0.05 & f('Torin_vs_Combined_FDR')<0.05 & f('Combined_vs_Media_logFC')*f('Torin_vs_Combined_logFC')<0
 flag[is.na(flag)]<-FALSE;results[[paste0(meth,'_retained')]]<-flag
}
orig<-read.csv(file.path(oldroot,'assets/original101_audit.csv'),check.names=FALSE)
stopifnot(identical(results$gene,orig$gene),identical(results$fractional_edgeR_retained,orig$meets_formal_response_criteria))
write.csv(results,file.path(out,'101个候选基因逐项比较.csv'),row.names=FALSE)
summary<-data.frame(method=c('fractional_edgeR','rounded_edgeR','fractional_voom'),tested=sapply(c('fractional_edgeR','rounded_edgeR','fractional_voom'),function(m)sum(!is.na(results[[paste0(m,'_Combined_vs_Media_FDR')]]))),retained=sapply(c('fractional_edgeR','rounded_edgeR','fractional_voom'),function(m)sum(results[[paste0(m,'_retained')]])),fixed48_retained=sapply(c('fractional_edgeR','rounded_edgeR','fractional_voom'),function(m)sum(results[[paste0(m,'_retained')]] & results$fixed48,na.rm=TRUE)),discordant_vs_original=sapply(c('fractional_edgeR','rounded_edgeR','fractional_voom'),function(m)sum(results[[paste0(m,'_retained')]] != results$fractional_edgeR_retained)))
write.csv(summary,file.path(out,'方法敏感性汇总.csv'),row.names=FALSE)
for(m in c('rounded_edgeR','fractional_voom')){flag<-results[[paste0(m,'_retained')]]!=results$fractional_edgeR_retained;write.csv(results[flag,,drop=FALSE],file.path(out,paste0(m,'_变动基因.csv')),row.names=FALSE)}
writeLines(capture.output(sessionInfo()),file.path(out,'运行环境.txt'))
writeLines(c(paste('source fractional entries before protein-coding aggregation',sum(abs(x-round(x))>1e-6)),paste('tested symbols common filter',sum(keep)),paste('sample group sizes',paste(names(table(grp)),as.numeric(table(grp)),collapse='; ')),'Rounding comparison is numerical sensitivity only; voom is alternative statistical model, not an independent experiment. Original gene-level expected counts lack transcript-level length offset in these analyses. Sample independence and length bias remain unresolved.'),file.path(out,'分析范围.txt'))
print(summary)
