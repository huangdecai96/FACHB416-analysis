aa<-grep('^--file=',commandArgs(FALSE),value=TRUE)
source(file.path(dirname(sub('^--file=','',aa[1])),'helpers/setup.R'))
source(file.path(SCRIPT_DIR,'helpers/phenotype_expression.R'))
base_theme<-function(){theme_classic(base_size=8.5,base_family='Arial')+theme(text=element_text(colour=ink),axis.text=element_text(size=12,colour=ink),axis.title=element_text(size=12.75,lineheight=.95),axis.title.y=element_text(margin=margin(r=2)),axis.title.x=element_text(margin=margin(t=3)),axis.line=element_blank(),axis.line.x=element_blank(),axis.line.y=element_blank(),axis.ticks=element_line(linewidth=.4,colour='black'),axis.ticks.length=unit(1,'mm'),panel.border=element_rect(fill=NA,colour='black',linewidth=.45),plot.margin=margin(3,3,3,3),plot.title=element_text(size=9,face='bold'),plot.subtitle=element_text(size=7.5),plot.tag=element_text(size=12,face='bold'),legend.title=element_text(size=7.5),legend.text=element_text(size=7.5),legend.key.size=unit(3,'mm'))}
polish<-function(p,headers=FALSE){
 p<-p+base_theme()+labs(caption=NULL)
 if(!headers)p<-p+labs(title=NULL,subtitle=NULL)
 for(i in seq_along(p$layers)){
  l<-p$layers[[i]]
  if(inherits(l$geom,'GeomText')||inherits(l$geom,'GeomTextRepel')){if(!is.null(l$aes_params$size))l$aes_params$size<-max(2.65,l$aes_params$size);l$aes_params$family<-'Arial'}
  if(inherits(l$geom,'GeomSegment')||inherits(l$geom,'GeomErrorbar')){if(!is.null(l$aes_params$linewidth))l$aes_params$linewidth<-max(.4,l$aes_params$linewidth)}
  p$layers[[i]]<-l
 };p
}
plot_metric<-function(m,tag){
 d<-raw[raw$metric==m,];s<-summ[summ$metric==m,];a<-cld[cld$metric==m,]
 d$x<-match(d$condition,co)+c(-.23,-.09,.05)[d$replicate];d$level<-substr(d$condition,1,1);s$x<-match(s$condition,co)+.22
 a<-merge(a,s,by=c('condition','metric'));tops<-tapply(d$value,d$condition,max);rng<-max(d$value,s$mean+s$se)
 a$y<-pmax(tops[a$condition],a$mean+a$se)+rng*.10;a$x<-match(a$condition,co)
 yl<-switch(m,growth_rate_r=expression(italic(r)~(d^{-1})),carrying_capacity_k=expression(atop('Carrying capacity',paste('(',10^4~cells~mL^{-1},')'))),day7_cell_volume=expression('Cell volume ('*mu*m^3*')'),day7_cpp='Cells per\nparticle (CPP)')
 p<-ggplot()+geom_point(data=d,aes(x,value,fill=condition,shape=level),size=2.1,stroke=.35,colour=ink)+geom_errorbar(data=s,aes(x=x,ymin=mean-se,ymax=mean+se),width=.13,linewidth=.4)+geom_point(data=s,aes(x,mean),shape=18,size=1.8)+geom_text(data=a,aes(x,y,label=letter),size=2.8,family='Arial')+scale_fill_manual(values=paper_style$condition,guide='none')+scale_shape_manual(values=paper_style$concentration_shape,guide='none')+scale_x_continuous(breaks=1:6,labels=co,limits=c(.5,6.6),expand=c(0,0))+scale_y_continuous(limits=c(0,max(a$y)*1.13),expand=c(0,0))+labs(tag=tag,x=NULL,y=yl,subtitle=sprintf('ANOVA P = %.3g',tests$anova_p[tests$metric==m]))+base_theme()
 b<-ggplot_build(p);stopifnot(all(abs(b$data[[1]]$x-d$x)<1e-8),all(abs(b$data[[1]]$y-d$value)<1e-8));p
}
plot_comp<-function(){d<-cs;d$condition<-factor(d$condition,levels=co);d$cell_class<-factor(d$cell_class,levels=c(1,2,4,8));d$label_y<-NA_real_;for(cc in co){ix<-which(d$condition==cc);ix<-ix[order(as.numeric(d$cell_class[ix]),decreasing=TRUE)];v<-100*d$mean[ix];d$label_y[ix]<-cumsum(v)-v/2};ggplot(d,aes(condition,100*mean,fill=cell_class))+geom_col(width=.88,colour='white',linewidth=.3)+geom_text(data=d[d$mean>.11,],aes(y=label_y,label=sprintf('%.1f%%',100*mean),colour=as.numeric(as.character(cell_class))>=4),size=2.65,fontface='bold',family='Arial',show.legend=FALSE)+scale_colour_manual(values=c('FALSE'='#000000','TRUE'='#FFFFFF'),guide='none')+scale_fill_manual(values=paper_style$colony,name='Cells per colony')+scale_y_continuous(breaks=c(0,50,100),expand=c(0,0))+coord_cartesian(ylim=c(0,100))+labs(tag='e',x='Nitrogen source and concentration',y='The proportion of\ncells in different\ncolonies (%)')+base_theme()+theme(legend.position='bottom',legend.key.width=unit(2.4,'mm'),legend.key.height=unit(2.4,'mm'),legend.spacing.x=unit(.5,'mm'),axis.title.x=element_text(size=11.55))}

comp_check<-ggplot_build(plot_comp());bc<-comp_check$data[[1]];lc<-comp_check$data[[2]]
for(i in seq_len(nrow(lc))){ix<-which(abs(bc$x-lc$x[i])<1e-8 & lc$y[i]>bc$ymin & lc$y[i]<bc$ymax);stopifnot(length(ix)==1,abs((bc$ymax[ix]-bc$ymin[ix])-as.numeric(sub('%','',lc$label[i])))<.051)}
write.table(lc[c('x','y','label')],file.path(ROOT,'verification/figure4e_label_centres.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
resp$label<-factor(resp$metric,levels=rev(metrics[1:4]),labels=rev(c('r','K','Cell volume','CPP')))
overview<-polish(plot_overview())+theme(legend.position='bottom',legend.box='horizontal',legend.key.height=unit(2,'mm'),legend.key.width=unit(2,'mm'),legend.title=element_text(size=7),legend.text=element_text(size=7),legend.spacing.x=unit(1,'mm'),axis.text.y=element_text(size=12),panel.border=element_blank(),axis.line=element_blank(),axis.ticks=element_blank())+guides(fill=guide_colourbar(display='rectangles',direction='horizontal',barwidth=unit(17,'mm'),barheight=unit(2,'mm'),title.position='top'),size=guide_legend(direction='horizontal',title.position='top',nrow=1))
cors$trait<-factor(cors$metric,levels=metrics,labels=c('r','K','Cell\nvolume\nDay 7','CPP\nDay 7','Single\ncells\nDay 7','Fv/Fm\nDay 7','\u03A6PSII\nDay 7'))
BODY_H<-38
make_grobs<-function(plots,pairs,fixed_width=FALSE){
 gs<-lapply(plots,ggplotGrob)
 for(g in gs)stopifnot(sum(g$layout$name=='panel')==1)
 for(ids in pairs){hh<-do.call(unit.pmax,lapply(gs[ids],function(g)g$heights));for(i in ids)gs[[i]]$heights<-hh}
 for(i in seq_along(gs)){pr<-gs[[i]]$layout[gs[[i]]$layout$name=='panel',];gs[[i]]$heights[pr$t]<-unit(BODY_H,'mm')}
 ids<-if(fixed_width)seq_along(gs)else 1:2
 # Align panel edges using an OUTER spacer; do not inflate the numeric tick-label column
 # with the long 'Cell volume' row label from the overview panel.
 pc<-gs[[1]]$layout$l[gs[[1]]$layout$name=='panel']
 left<-vapply(gs[ids],function(g)convertWidth(sum(g$widths[seq_len(pc-1)]),'mm',TRUE),numeric(1))
 right<-vapply(gs[ids],function(g)convertWidth(sum(g$widths[(pc+1):length(g$widths)]),'mm',TRUE),numeric(1))
 ml<-max(left);mr<-max(right);pw<-88-ml-mr;stopifnot(pw>=25)
 for(j in seq_along(ids)){i<-ids[j];g<-gs[[i]];g$widths[1]<-g$widths[1]+unit(ml-left[j],'mm');g$widths[length(g$widths)]<-g$widths[length(g$widths)]+unit(mr-right[j],'mm');g$widths[pc]<-unit(pw,'mm');gs[[i]]<-g}
 gs
}
cairo_pdf(file.path(ROOT,'work/layout_measure.pdf'),width=180/25.4,height=200/25.4,pointsize=8.5,family='Arial')
G4<-make_grobs(c(lapply(1:4,function(i)plot_metric(metrics[i],letters[i])),list(plot_comp(),overview)),list(1:2,3:4,5:6),TRUE)
pr_e<-G4[[5]]$layout[G4[[5]]$layout$name=='panel',];bar_inner_mm<-convertWidth(G4[[5]]$widths[pr_e$l],'mm',TRUE)*.88/6.2-.3
pushViewport(viewport(gp=gpar(fontfamily='Arial',fontface='bold',fontsize=2.65*.pt)));label_width_mm<-max(convertWidth(stringWidth(lc$label),'mm',TRUE));popViewport();stopifnot(label_width_mm<bar_inner_mm)
write.table(data.frame(max_label_width_mm=label_width_mm,bar_inner_width_mm=bar_inner_mm),file.path(ROOT,'verification/percentage_width_check.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
# Match the d-panel y-title baseline to panel b without moving either plot body.
ytitle_x<-function(g){idx<-which(g$layout$name=='ylab-l');col<-g$layout$l[idx];child<-g$grobs[[idx]]$children[[1]];pushViewport(viewport(width=g$widths[col],height=unit(BODY_H,'mm')));xx<-(convertX(grobX(child,0),'mm',TRUE)+convertX(grobX(child,180),'mm',TRUE))/2;popViewport();convertWidth(sum(g$widths[seq_len(col-1)]),'mm',TRUE)+xx}
title_b<-ytitle_x(G4[[2]]);title_d_old<-ytitle_x(G4[[4]])
idx_d<-which(G4[[4]]$layout$name=='ylab-l');G4[[4]]$grobs[[idx_d]]$children[[1]]$x<-G4[[4]]$grobs[[idx_d]]$children[[1]]$x+unit(title_b-title_d_old,'mm')
stopifnot(abs(ytitle_x(G4[[4]])-title_b)<1e-8)
write.table(data.frame(panel=c('b','d'),y_title_x_mm=c(title_b,ytitle_x(G4[[4]]))),file.path(ROOT,'verification/figure4_bd_title_alignment.tsv'),sep='\t',row.names=FALSE,quote=FALSE)

shift<-polish(plot_shift())+theme(legend.position='bottom',panel.border=element_blank(),axis.line.x=element_line(linewidth=.45,colour='black'),axis.line.y=element_blank(),axis.ticks.y=element_blank())
# Origin-centred, unframed PCA axes; all original sample coordinates are retained.
 pca<-polish(E5$plot_pca(.85))+labs(tag='a')+
 geom_hline(yintercept=0,linewidth=.32,linetype='dashed',colour='#7F8994')+geom_vline(xintercept=0,linewidth=.32,linetype='dashed',colour='#7F8994')+
 annotate('segment',x=c(-50,50),xend=c(-50,50),y=-2,yend=2,linewidth=.35)+
 annotate('segment',y=c(-50,50),yend=c(-50,50),x=-2,xend=2,linewidth=.35)+
 annotate('text',x=c(-50,50),y=-7,label=c('-50','50'),size=3.9,family='Arial')+
 annotate('text',x=-5,y=c(-50,50),label=c('-50','50'),hjust=1,size=3.9,family='Arial')+
 annotate('text',x=-4,y=-7,label='0',hjust=1,size=3.9,family='Arial')+
 scale_x_continuous(limits=c(-100,100),breaks=NULL,expand=c(0,0))+
 scale_y_continuous(limits=c(-100,100),breaks=NULL,expand=c(0,0))+
 theme(panel.border=element_blank(),axis.line=element_blank(),axis.ticks=element_blank(),axis.text=element_blank())
 candidates<-polish(plot_candidates())+theme(panel.border=element_blank(),axis.line.x=element_line(linewidth=.45,colour='black'),axis.line.y=element_blank(),axis.ticks.y=element_blank())
heat<-polish(plot_heat())+labs(tag='d')+scale_y_discrete(expand=expansion(add=0))+scale_x_discrete(expand=expansion(add=0))+theme(panel.border=element_blank(),axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_text(size=11.7),axis.text.y=element_text(size=12),legend.position='right',legend.key.height=unit(5,'mm'))
G5<-make_grobs(list(pca,shift,candidates,heat),list(1:2,3,4),FALSE)
# Align the overall left edge of a with c/d without changing either top panel's data width.
pad_a<-G5[[1]]$widths[1]-unit(3,'pt');G5[[1]]$widths[1]<-unit(3,'pt');G5[[1]]$widths[length(G5[[1]]$widths)]<-G5[[1]]$widths[length(G5[[1]]$widths)]+pad_a
# Panel letters are anchored to the overall figure grid, independently of tick/label widths.
for(i in 1:4){j<-which(G5[[i]]$layout$name=='tag');G5[[i]]$grobs[[j]]<-zeroGrob()}

dim_mm<-function(g)convertHeight(sum(g$heights),'mm',valueOnly=TRUE)
row4<-vapply(G4[c(1,3,5)],dim_mm,numeric(1));row5<-vapply(G5[c(1,3,4)],dim_mm,numeric(1))
Y4<-c(0,cumsum(row4+4)[1:2]);Y5<-c(0,cumsum(row5+4)[1:2]);H4<-sum(row4)+8+2;H5<-sum(row5)+8+2
cat('Measured heights',H4,H5,'mm\n')
trace<-list()
for(fig in c(4,5)){gs<-if(fig==4)G4 else G5;ys<-if(fig==4)rep(Y4,each=2)else c(Y5[1],Y5[1],Y5[2],Y5[3]);for(i in seq_along(gs)){g<-gs[[i]];p<-g$layout[g$layout$name=='panel',];trace[[paste(fig,i)]]<-data.frame(figure=fig,panel=letters[i],body_height_mm=convertHeight(g$heights[p$t],'mm',TRUE),body_width_mm=if(fig==4||i<=2)convertWidth(g$widths[p$l],'mm',TRUE)else NA_real_,body_top_mm=ys[i]+convertHeight(sum(g$heights[seq_len(p$t-1)]),'mm',TRUE))}}
gaptrace<-lapply(seq_along(G4),function(i){g<-G4[[i]];panel<-g$layout[g$layout$name=='panel',];yl<-g$layout[g$layout$name=='ylab-l',];data.frame(panel=letters[i],panel_left_mm=convertWidth(sum(g$widths[seq_len(panel$l-1)]),'mm',TRUE),ylabel_column_width_mm=convertWidth(g$widths[yl$l],'mm',TRUE))})
 write.table(do.call(rbind,gaptrace),file.path(ROOT,'verification/figure4_label_geometry.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
 dev.off();write.table(do.call(rbind,trace),file.path(ROOT,'verification/panel_geometry_mm.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
drawg<-function(g,x,y,w){h<-dim_mm(g);pushViewport(viewport(x=unit(x,'native'),y=unit(y,'native'),width=unit(w,'mm'),height=unit(h,'mm'),just=c('left','top')));grid.draw(g);popViewport()}
draw4<-function(){page(H4);for(i in 1:6)drawg(G4[[i]],if(i%%2)0 else 91,Y4[ceiling(i/2)],89);popViewport()}
draw5<-function(){page(H5);drawg(G5[[1]],0,Y5[1],89);drawg(G5[[2]],91,Y5[1],89);drawg(G5[[3]],0,Y5[2],180);drawg(G5[[4]],0,Y5[3],180);for(i in 1:4)grid.text(letters[i],unit(if(i==2)92 else 1,'native'),unit(c(Y5[1],Y5[1],Y5[2],Y5[3])[i]+2,'native'),just=c('left','top'),gp=gpar(fontfamily='Arial',fontsize=12,fontface='bold',col=ink));popViewport()}
export_fig('Figure4',draw4,H4)
