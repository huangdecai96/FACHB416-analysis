aa<-grep('^--file=',commandArgs(FALSE),value=TRUE)
source(file.path(dirname(sub('^--file=','',aa[1])),'helpers/setup.R'))
source(file.path(SCRIPT_DIR,'helpers/time_series.R'))
photo<-rd('Photosynthesis','photosynthesis_all_days.tsv');ss<-do.call(rbind,lapply(c('fvfm','phipsii'),function(m)do.call(rbind,lapply(split(photo,list(photo$condition,photo$day),drop=TRUE),function(q)data.frame(condition=q$condition[1],day=q$day[1],metric=m,n=nrow(q),mean=mean(q[[m]]),se=sd(q[[m]])/sqrt(nrow(q)))))));stopifnot(nrow(ss)==120,all(ss$n==3));ss$condition<-factor(ss$condition,levels=co)
onephoto<-function(m,yl,tag){ptime(ss[ss$metric==m,],yl)+geom_vline(xintercept=7,linetype='dotted',linewidth=.35,colour='#9CA3AF')+scale_x_continuous(breaks=c(1,3,5,7,9,11,13,15,17),limits=c(.65,17.35),expand=c(0,0))+scale_y_continuous(breaks=seq(0,.8,.2),limits=c(0,.85),expand=c(0,0))+labs(tag=tag)}
pphoto<-(onephoto('fvfm',expression(F[v]/F[m]),'a')+onephoto('phipsii',expression(Phi[PSII]),'b')+plot_layout(guides='collect'))&theme(legend.position='bottom')
export_fig('FigureS7',function()print(pphoto),105)
