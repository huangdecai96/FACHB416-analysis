pca <- rd('pca.tsv'); deg <- rd('deg.tsv')
theme5 <- function(scale=1) theme_paper(scale)+theme(
  plot.title=element_text(size=8.4*scale,face='bold',margin=margin(l=14*scale,b=3)),
  plot.subtitle=element_text(size=6.2*scale,margin=margin(b=5)),
  strip.background=element_blank(),strip.text=element_text(size=6.8*scale,face='bold'),
  panel.spacing=grid::unit(4,'mm'),plot.margin=margin(5,5,4,5))
shape_colours <- function() list(scale_fill_manual(values=paper_style$condition,guide='none'),
 scale_shape_manual(values=paper_style$concentration_shape,guide='none'))
p_fmt <- function(x) ifelse(x>=.001,formatC(x,format='f',digits=3),formatC(x,format='e',digits=1))

plot_pca <- function(scale=1) {
  dat <- pca;dat$conc <- substr(dat$condition,1,1)
  centre <- aggregate(cbind(PC1,PC2)~condition,dat,mean)
  centre$dx <- c(HA=-4,HN=-3,HU=0,LA=6,LN=-6,LU=13)[centre$condition]
  centre$dy <- c(HA=16,HN=-16,HU=16,LA=11,LN=-12,LU=-12)[centre$condition]
  ggplot(dat,aes(PC1,PC2))+geom_point(aes(fill=condition,shape=conc),
    size=2.35*scale,stroke=.35,colour=paper_style$text)+shape_colours()+
    geom_text(data=centre,aes(x=PC1+dx,y=PC2+dy,label=condition),size=2.5*scale,
      family=paper_style$font,fontface='bold',colour=paper_style$text)+
    scale_x_continuous(limits=c(-98,70),breaks=c(-75,-25,25),expand=c(0,0))+
    scale_y_continuous(limits=c(-100,72),breaks=c(-75,-25,25,65),expand=c(0,0))+
    labs(tag='a',title='Transcriptome PCA',subtitle='18 RNA samples; n = 3 per condition',
      x=sprintf('PC1 (%.1f%%)',pca$PC1_percent[1]),y=sprintf('PC2 (%.1f%%)',pca$PC2_percent[1]))+theme5(scale)
}
contrast_order <- c('HN_vs_LN','HU_vs_LU','HA_vs_LA','HU_vs_HN','HA_vs_HN','HA_vs_HU')
plot_deg <- function(scale=1) {
 dat<-deg;dat$contrast<-factor(dat$contrast,levels=rev(contrast_order))
 dat$direction<-factor(dat$direction,levels=c('Down','Up'))
 ggplot(dat,aes(n_plot,contrast,fill=direction))+
  annotate('rect',xmin=-4650,xmax=3900,ymin=2.55,ymax=3.45,fill='#EEF5F3',colour=NA)+
  geom_col(width=.62,colour='white',linewidth=.22)+geom_vline(xintercept=0,linewidth=.25,colour=paper_style$secondary)+
  geom_text(aes(x=n_plot+ifelse(n_plot<0,-130,130),label=n,hjust=ifelse(n_plot<0,1,0)),
    size=2.35*scale,family=paper_style$font,colour=paper_style$text)+
  scale_fill_manual(values=c(Down='#A84F65',Up='#2F73AE'),name=NULL)+
  scale_x_continuous(limits=c(-4650,3900),breaks=c(-4000,-2000,0,2000),labels=abs,expand=c(0,0))+
  scale_y_discrete(labels=gsub('_vs_',' vs ',rev(contrast_order)))+
  labs(tag='e',title='Support: genome-wide DEG context',subtitle='Editable support page; not a main-claim panel',x='DEGs (first condition relative to second)',y=NULL)+
  theme5(scale)+theme(legend.position='bottom',legend.key.height=grid::unit(2,'mm'),legend.margin=margin(0,0,0,0))
}
