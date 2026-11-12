aa<-grep('^--file=',commandArgs(FALSE),value=TRUE)
source(file.path(dirname(sub('^--file=','',aa[1])),'helpers/setup.R'))
source(file.path(SCRIPT_DIR,'helpers/time_series.R'))
s<-rd('Morphology','CPP_summary.tsv');stopifnot(nrow(s)==54,all(s$n==3));s$condition<-factor(s$condition,levels=co)
s<-s[order(s$day,s$condition),] # Preserve the actual final point-drawing order.
pcpp<-ptime(s,'Cells per particle (CPP)')+geom_hline(yintercept=1,linetype='dotted',colour='#B7B7B7',linewidth=.3)+scale_y_continuous(breaks=1:5,limits=c(.78,5.22),expand=c(0,0))
export_fig('FigureS6',function()print(pcpp),108)
