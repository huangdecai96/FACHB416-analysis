rd_data<-rd
rd<-function(f)rd_data('Expression',sub('^Figure5/','',f))
dd<-function(f)rd_data('Phenotype_expression',f)
wt<-function(x,n)write.table(x,file.path(ROOT,'verification',n),sep='	',row.names=FALSE,quote=FALSE)
paper_style$font_size<-6.5;paper_style$axis_text_size<-6.2;paper_style$title_size<-7;paper_style$panel_label_size<-7
co<-c('LN','HN','LU','HU','LA','HA');sets<-c('E1 growth-maintenance','S0 all','S0B precursor transport','C1/C2 adaptive reserve')
set_short<-c('E1 maintenance','S0 interface','S0B transport','C1/C2 reserve');set_cols<-setNames(c('#DB8C2F','#3485AE','#8862A5','#6B7280'),sets)
metrics<-c('growth_rate_r','carrying_capacity_k','day7_cell_volume','day7_cpp','day7_single_cell_fraction','day7_fvfm','day7_phipsii')
ink<-'#20242B';rose<-'#A84F65';green<-'#1B927F'
sumval<-function(d,keys,value){do.call(rbind,lapply(split(d,interaction(d[keys],drop=TRUE)),function(q){a<-q[1,keys,drop=FALSE];v<-q[[value]];a$n<-length(v);a$mean<-mean(v);a$sd<-sd(v);a$se<-sd(v)/sqrt(length(v));a}))}
raw<-dd('plot_replicates.tsv');summ<-dd('trait_summary.tsv');tests<-dd('figure4_global_tests.tsv');cld<-dd('figure4_letters.tsv');cors<-dd('figure5_correlations.tsv')
cc<-rd_data('Morphology','colony_counts_all_days.tsv');cc<-cc[cc$day==7,]
cs<-do.call(rbind,lapply(c(1,2,4,8),function(k)do.call(rbind,lapply(split(cc,cc$condition),function(q){v<-q[[paste0('class',k,'_cell_fraction')]];data.frame(condition=q$condition[1],cell_class=k,n=length(v),mean=mean(v),sd=sd(v),se=sd(v)/sqrt(length(v)))}))))
ms<-rd('Figure5/module_summary.tsv');scores<-rd('Figure5/module_sample_scores.tsv')
E5<-new.env(parent=globalenv());E5$rd<-function(f)rd(f)
source(file.path(SCRIPT_DIR,'helpers/expression_components.R'),local=E5)
sh<-ms[ms$score_module%in%sets&ms$condition%in%c('HN','HU'),];sh$yy<-match(sh$score_module,rev(sets));sh$dy<-ifelse(sh$condition=='HN',-.08,.08)
raw5<-scores[scores$condition%in%c('HN','HU'),];raw5$yy<-match(raw5$score_module,rev(sets));raw5$yraw<-raw5$yy+ifelse(raw5$condition=='HN',-.17,.17)
plot_shift<-function(){con<-merge(sh[sh$condition=='HN',],sh[sh$condition=='HU',],by='score_module',suffixes=c('_HN','_HU'));ann<-rd('Figure5/module_contrasts.tsv');ann$yy<-match(ann$score_module,rev(sets));ann$lab<-sprintf('%+.2f %s',ann$delta_A_minus_B,ann$sig)
 ggplot()+geom_segment(data=con,aes(x=mean_score_HN,xend=mean_score_HU,y=yy_HN,yend=yy_HU),linewidth=.7,colour='#CFD6DE')+geom_point(data=raw5,aes(score,yraw,colour=score_module,shape=condition),size=.8,alpha=.6)+geom_errorbar(data=sh,aes(y=yy+dy,xmin=mean_score-se,xmax=mean_score+se),orientation='y',width=.08,linewidth=.3)+geom_point(data=sh,aes(mean_score,yy+dy,fill=condition,colour=score_module),shape=21,size=2.1,stroke=.5)+geom_text(data=ann,aes(x=9.02,y=yy,label=lab,colour=score_module),hjust=0,size=2.15,fontface='bold',family='Arial')+scale_colour_manual(values=set_cols,guide='none')+scale_shape_manual(values=c(HN=1,HU=16),guide='none')+scale_fill_manual(values=c(HN='white',HU='#697383'),name=NULL)+scale_y_continuous(breaks=4:1,labels=set_short,limits=c(.5,4.5))+scale_x_continuous(limits=c(6.5,9.85),breaks=c(7,8,9))+labs(tag='b',title='Candidate-set response',subtitle='HU - HN; means +/- SE and individual RNA samples',x='Mean VST score',y=NULL)+theme_paper()+theme(legend.position='bottom',axis.line.y=element_blank(),axis.ticks.y=element_blank())}
cors$set<-factor(cors$score_module,levels=rev(sets),labels=rev(set_short));cors$trait<-factor(cors$metric,levels=metrics,labels=c('Fitted\nr','Fitted\nK','Cell volume\nDay 7','CPP\nDay 7','Unicellular\ncell fraction','Fv/Fm\nDay 7','\u03A6PSII\nDay 7'))
plot_heat<-function(){ggplot(cors,aes(trait,set,fill=rho))+geom_tile(colour='white',linewidth=.6)+geom_text(aes(label=paste0(sprintf('%.2f',rho),ifelse(q_BH28<.05,'*',''))),size=2.15,family='Arial')+scale_fill_gradient2(low='#B86166',mid='#F6F5F1',high='#4D80AA',limits=c(-1,1),name='Spearman\nrho')+guides(fill=guide_colourbar(display='rectangles'))+labs(tag='c',title='Candidate sets and phenotype',subtitle='Six condition means; morphology and photosynthesis measured on day 7',x=NULL,y=NULL)+theme_paper()+theme(axis.line=element_blank(),axis.ticks=element_blank(),legend.position='right',legend.key.height=unit(4,'mm'),axis.text.x=element_text(size=5.9),axis.text.y=element_text(size=6))}
cand<-rd('seven_candidates_DESeq2.tsv');candexpr<-rd('seven_candidates_expression.tsv')
cand$strict<-cand$padj<.05 & abs(cand$log2FoldChange)>=1
stopifnot(sum(cand$strict)==6,!cand$strict[cand$gene_id=='g11586'])
candidate_order<-c('g8640','g43','g7401','g11251','g12581','g12816','g11586')
cand$y<-match(cand$gene_id,rev(candidate_order))
candidate_labels<-c(g8640='g8640 | PsbP-domain',g43='g43 | LHC',g7401='g7401 | AMT',
 g11251='g11251 | GOGAT',g12581='g12581 | Protein kinase',
 g12816='g12816 | Nucleotide-sugar transport (pred.)',g11586='g11586 | Nucleotide-sugar transport (pred.)')
candidate_col<-c('Photo.'='#3379A9','N transport'='#258B7B','N assimilation'='#258B7B','Reg.'='#747B87','S0B'='#8862A5')
cand$padj_label<-formatC(cand$padj,format='e',digits=2)
cand$effect_label<-sprintf('%+.2f',cand$log2FoldChange)
wt(cand,'figure5_seven_candidates.tsv')
wt(sumval(candexpr,c('gene_id','condition'),'vst'),'seven_candidates_vst_summary.tsv')

plot_candidates<-function(){
 ggplot(cand,aes(y=y))+
 geom_vline(xintercept=0,colour='#8D959E',linetype='dashed',linewidth=.3)+
 geom_segment(aes(x=0,xend=log2FoldChange,yend=y,colour=functional_class),linewidth=.6)+
 geom_point(aes(x=log2FoldChange,colour=functional_class,fill=functional_class,shape=strict),size=2.2,stroke=.5)+
 geom_text(aes(x=ifelse(log2FoldChange>0,log2FoldChange+.5,log2FoldChange-.5),
  label=effect_label,hjust=ifelse(log2FoldChange>0,0,1)),size=2.15,family='Arial')+
 geom_text(aes(x=9.3,label=padj_label),hjust=1,size=2.15,family='Arial',colour=ink)+
 annotate('text',x=9.3,y=7.9,label='DESeq2 adj. P',hjust=1,size=2.1,fontface='bold',family='Arial')+
 scale_colour_manual(values=candidate_col,guide='none')+scale_fill_manual(values=candidate_col,guide='none')+
 scale_shape_manual(values=c('TRUE'=21,'FALSE'=1),guide='none')+
 scale_y_continuous(breaks=7:1,labels=candidate_labels[candidate_order],limits=c(.4,8.15),expand=c(0,0))+
 scale_x_continuous(limits=c(-12,9.6),breaks=c(-10,-5,0,5),expand=c(0,0))+
 labs(tag='c',title='Selected gene responses',subtitle='HU relative to HN; selected examples across functional annotations',
  x=expression('HU/HN log'[2]*' fold change'),y=NULL)+theme_paper()+
 theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_text(size=6),panel.grid=element_blank())
}
resp<-summ[summ$metric%in%metrics[1:4],];ref<-setNames(resp$mean[resp$condition=='HN'],resp$metric[resp$condition=='HN']);resp$log2_ratio<-log2(resp$mean/ref[resp$metric]);resp$condition<-factor(resp$condition,levels=co);resp$label<-factor(resp$metric,levels=rev(metrics[1:4]),labels=rev(c('Fitted r (1-17 d)','Fitted K (1-17 d)','Volume (D7)','CPP (D7)')));wt(resp[c('condition','metric','mean','log2_ratio')],'figure4_overview.tsv')
lim<-max(2.4,ceiling(max(abs(resp$log2_ratio))*2)/2)
plot_overview<-function(){ggplot(resp,aes(condition,label))+geom_point(aes(size=abs(log2_ratio),fill=log2_ratio),shape=21,colour='#65717C',stroke=.18)+scale_fill_gradient2(low=rose,mid='#F3F4F3',high=green,midpoint=0,limits=c(-lim,lim),name='log2 ratio\nto HN')+guides(fill=guide_colourbar(display='rectangles'))+scale_size_continuous(range=c(.8,3.4),limits=c(0,lim),breaks=c(.5,1,2),name='Absolute\nchange')+labs(tag='f',title='Overall phenotype response',subtitle='Group means relative to HN',x=NULL,y=NULL)+theme_paper()+theme(legend.position='right',legend.box='vertical',legend.key.height=unit(2.3,'mm'),axis.text.y=element_text(size=5.8))}
