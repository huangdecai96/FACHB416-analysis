aa<-grep('^--file=',commandArgs(FALSE),value=TRUE)
source(file.path(dirname(sub('^--file=','',aa[1])),'helpers/setup.R'))
chr<-rd('Figure1','chromosomes.tsv');chr<-chr[match(sprintf('Chr%02d',1:17),chr$chr),]
genes1<-rd('Figure1','genes.tsv');rep1<-rd('Figure1','repeats_gc.tsv');te<-rd('Figure1','telomeres.tsv');cen<-rd('Figure1','centromeres.tsv');s01<-rd('Figure1','s0_positions.tsv')
stopifnot(nrow(chr)==17,sum(chr$length_bp)==101858869)
genes1<-genes1[genes1$seqid%in%chr$chr,];rep1<-rep1[rep1$seqid%in%chr$chr,];te<-te[match(chr$chr,te$New_ID),];cen<-cen[match(chr$chr,cen$seqid),];s01<-s01[s01$chr%in%chr$chr,]
# Retained final centromere-panel style (before later frame changes).
base_theme<-function(){theme_classic(base_size=8.5,base_family='Arial')+theme(text=element_text(colour=ink),axis.text=element_text(size=8,colour=ink),axis.title=element_text(size=8.5),axis.line=element_line(linewidth=.45,colour='black'),axis.ticks=element_line(linewidth=.4,colour='black'),axis.ticks.length=unit(1,'mm'),panel.border=element_rect(fill=NA,colour='black',linewidth=.4),plot.margin=margin(3,3,3,3),plot.title=element_text(size=9,face='bold'),plot.subtitle=element_text(size=7.5),plot.tag=element_text(size=12,face='bold'),legend.title=element_text(size=7.5),legend.text=element_text(size=7.5),legend.key.size=unit(3,'mm'))}
pcen<-ggplot(data.frame(chr=chr$chr,len=chr$length_mb,x=cen$midpoint_mb,y=17:1,score=cen$centromere_like_score),aes(y=y))+geom_segment(aes(x=0,xend=len,yend=y),colour='#C3CCD3',linewidth=1)+geom_point(aes(x=x,size=score),shape=15,colour='#86569C')+scale_size_continuous(range=c(1.8,3.2),name='Candidate score')+scale_y_continuous(breaks=17:1,labels=chr$chr)+labs(x='Chromosome position (Mb)',y=NULL)+base_theme()+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),panel.border=element_blank())
export_fig('FigureS2',function()print(pcen),140)
