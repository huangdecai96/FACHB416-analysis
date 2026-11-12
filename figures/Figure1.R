# R-generated scientific panels before the retained PowerPoint layout edits.
# See Non_R_figures.md; this is not a pixel-exact final-composition generator.
aa<-grep('^--file=',commandArgs(FALSE),value=TRUE)
source(file.path(dirname(sub('^--file=','',aa[1])),'helpers/setup.R'))
suppressPackageStartupMessages(library(circlize))
chr<-rd('Figure1','chromosomes.tsv');chr<-chr[match(sprintf('Chr%02d',1:17),chr$chr),]
genes1<-rd('Figure1','genes.tsv');rep1<-rd('Figure1','repeats_gc.tsv');te<-rd('Figure1','telomeres.tsv');cen<-rd('Figure1','centromeres.tsv');s01<-rd('Figure1','s0_positions.tsv')
stopifnot(nrow(chr)==17,sum(chr$length_bp)==101858869)
genes1<-genes1[genes1$seqid%in%chr$chr,];rep1<-rep1[rep1$seqid%in%chr$chr,];te<-te[match(chr$chr,te$New_ID),];cen<-cen[match(chr$chr,cen$seqid),];s01<-s01[s01$chr%in%chr$chr,]
draw1<-function(){
 old<-par(no.readonly=TRUE);on.exit(par(old));layout(matrix(c(1,2),1),widths=c(1.82,1));par(family='Arial',ps=7,mar=c(0,0,1.2,0),xpd=NA)
 circos.clear();circos.par(start.degree=90,gap.after=c(rep(1.4,16),7),track.margin=c(.003,.003),cell.padding=c(0,0,0,0),canvas.xlim=c(-1.17,1.17),canvas.ylim=c(-1.17,1.17),points.overflow.warning=FALSE)
 circos.initialize(chr$chr,xlim=cbind(0,chr$length_bp))
 cpal<-setNames(rep(c('#557D98','#91ADBF'),length.out=17),chr$chr)
 circos.trackPlotRegion(ylim=c(0,1),track.height=.065,bg.border=NA,panel.fun=function(x,y){c<-CELL_META$sector.index;len<-chr$length_bp[chr$chr==c];circos.rect(0,0,len,1,col=cpal[c],border='white',lwd=.4);circos.text(len/2,2.65,c,cex=.92,font=2,facing='clockwise',niceFacing=TRUE)})
 rmax<-ceiling(max(rep1$softmasked_percent)/10)*10;gmax<-ceiling(max(genes1$gene_density_per_mb)/100)*100;gcr<-range(rep1$gc_percent)
 circos.genomicTrack(rep1[,c('seqid','bin_start','bin_end','softmasked_percent')],ylim=c(0,rmax),track.height=.15,bg.border='#DEE3E7',panel.fun=function(region,value,...){circos.genomicRect(region,value,ybottom=0,ytop=value[[1]],col='#D58B63',border=NA)})
 circos.genomicTrack(genes1[,c('seqid','bin_start','bin_end','gene_density_per_mb')],ylim=c(0,gmax),track.height=.15,bg.border='#DEE3E7',panel.fun=function(region,value,...){circos.genomicRect(region,value,ybottom=0,ytop=value[[1]],col='#429B89',border=NA)})
 circos.genomicTrack(rep1[,c('seqid','bin_start','bin_end','gc_percent')],ylim=gcr,track.height=.09,bg.border='#DEE3E7',panel.fun=function(region,value,...){circos.genomicLines(region,value[[1]],col='#68747E',lwd=.55)})
 circos.genomicTrack(s01[,c('chr','start','end')],ylim=c(0,1),track.height=.055,bg.border='#DEE3E7',panel.fun=function(region,value,...){mid<-(region[[1]]+region[[2]])/2;circos.segments(mid,rep(.08,length(mid)),mid,rep(.92,length(mid)),col=rose,lwd=.9)})
 text(0,.29,'17 nuclear chromosomes',cex=1,font=2);text(0,.19,'101.86 Mb',cex=.95)
 legend(x=-.36,y=.08,legend=c(sprintf('Repeats: 0-%d%%',rmax),sprintf('Gene density: 0-%d / Mb',gmax),sprintf('GC: %.1f-%.1f%%',gcr[1],gcr[2]),'S0 candidate positions'),col=c('#D58B63','#429B89','#68747E',rose),lwd=c(3,3,1,3),bty='n',cex=.85,y.intersp=1.2,adj=0)
 mtext('a  Nuclear genome overview',side=3,line=-.2,adj=.06,cex=1,font=2);text(0,-1.26,'500-kb windows; tracks outside to inside',cex=.85)
 circos.clear();par(mar=c(3.9,2.8,2.8,.6),ps=7)
 plot(NA,xlim=c(0,15),ylim=c(.4,17.6),xlab='Chromosome position (Mb)',ylab='',yaxt='n',bty='n',cex.lab=.94,cex.axis=.88,xaxs='i')
 axis(2,at=17:1,labels=chr$chr,las=1,tick=FALSE,cex.axis=.86)
 for(i in 1:17){y<-18-i;co<-if(te$Telomere_class[i]=='Single telomere')'#B78A31'else'#56646E';segments(0,y,te$Length_Mb[i],y,col=co,lwd=1.0);points(c(0,te$Length_Mb[i]),c(y,y),pch=21,bg=ifelse(c(te$head_telomere[i],te$tail_telomere[i])=='yes',ink,'white'),cex=.55,col=ink)}
 mtext('b  Telomere-end evidence',side=3,line=1.4,adj=0,cex=1,font=2)
 mtext('15 double-end; 2 single-end',side=3,line=.2,adj=0,cex=.87)
 legend('bottomright',legend=c('Motif present','No terminal motif'),pch=21,pt.bg=c(ink,'white'),col=ink,bty='n',cex=.85)
}


export_fig('Figure1',draw1,125)
