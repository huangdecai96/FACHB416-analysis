families<-rd('representative_families_summary.tsv');counts<-rd('representative_families_12species_long.tsv')
stable<-rd('stage_focal_stable_counts.tsv');runs<-rd('stage_focal_36_runs.tsv')
chain<-rd('S0B_candidate_member_chain.tsv')
stopifnot(nrow(families)==18,nrow(counts)==216,all(table(counts$family_id)==12),nrow(runs)==36,all(runs$contraction_dominant==1))
stopifnot(all(stable$stable_six_contracted==c(156,128)),all(stable$stable_six_expanded==c(49,46)))
stopifnot(sum(families$original_all_six_stable=='True')==17,families$original_models_supporting[families$family_id=='OG0000288']==4)
med<-aggregate(copies~family_id,counts[counts$species!='FCHB416',],median)
check<-merge(med,families,by='family_id');stopifnot(all(check$copies==check$other11_median))
stopifnot(all(families$FCHB416-families$other11_median==families$difference_vs_other11_median))
# Original author-approved order, followed by the newly added S0B candidate.
family_order<-c('OG0000007','OG0000212','OG0000019','OG0002200','OG0000195','OG0000001','OG0000011',
 'OG0000052','OG0000170','OG0006706','OG0000274','OG0000185','OG0000073','OG0000485','OG0008326','OG0007850','OG0000015','OG0000288')
short_names<-c('Protein kinase','Beta-helix','Peptidase C1','Peptidase C1','Peptidase S8','HA2-associated','MYND zinc finger',
 'Rieske-associated','Carbonic anhydrase','Sugar transporter','GT32-associated','RRM1-associated','Peptidase M43','Peptidase S28',
 'Peptidase C1 / I29','HA2 / helicase (mixed)','Histone-associated','S0B precursor transport')
species_order<-c('Araciborskii_CCAP276_35','Tdistendus_CCAP276_37','NREL46BD3','DOE0152z','FCHB416',
 'Tobliquus_SAG276_10','UTEX393','Ppectinatus_SAG2003','Shindakii_SAG47_86','Sobtusus_CCAP276_36',
 'Sproducto_capitatus_SAG21_81','Cmucronulatus_SAG2332')
species_labels<-c('T. sp. CCAP 276/35','T. distendus','T. sp. NREL46B-D3','DOE0152z','FACHB-416','SAG 276-10','UTEX 393',
 'T. sp. SAG 2003','S. hindakii','S. obtusus','S. producto-capitatus','C. mucronulatus')
f<-families[match(family_order,families$family_id),]; f$short_label<-short_names
f$y<-c(seq(19.3,13.3,-1),seq(11.8,2.8,-1),1.1)
counts$x<-match(counts$species,species_order);counts$y<-f$y[match(counts$family_id,f$family_id)]
counts$z_raw<-ave(counts$copies,counts$family_id,FUN=function(v)as.numeric(scale(v)))
stopifnot(all(is.finite(counts$z_raw)))
counts$z<-pmax(-2.5,pmin(2.5,counts$z_raw))
counts$text_colour<-ifelse(abs(counts$z)>1.3,'white',paper_style$text)
write.table(counts,file.path(output_dir,'heatmap_plotted_values.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(data.frame(check=c('family_rows','heatmap_cells','original_contracted','original_expanded','rescue_contracted','rescue_expanded','contraction_dominant_runs','S0B_six_model_support'),
 value=c(18,216,156,49,128,46,36,4)),file.path(output_dir,'numerical_checks.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
theme3<-function(scale=1)theme_paper(scale)+theme(plot.title=element_text(size=8.3*scale,face='bold',margin=margin(l=13*scale,b=3)),
 plot.subtitle=element_text(size=6.15*scale,margin=margin(b=5)),plot.margin=margin(5,5,5,5))

plot_heatmap<-function(scale=1){
 focal<-counts[counts$species=='FCHB416',]
 groups<-data.frame(ymin=c(12.8,2.3,.6),ymax=c(19.8,12.3,1.6),y=c(16.3,7.3,1.1),
  label=c('Stable\ncontractions','Stable\nexpansions','S0B\ncandidate'),colour=c(pal['contraction'],pal['expansion'],pal['candidate']))
 ggplot(counts,aes(x,y,fill=z))+geom_tile(width=1,height=1,colour='white',linewidth=.22)+
  geom_rect(data=data.frame(ymin=groups$ymin,ymax=groups$ymax),aes(xmin=4.5,xmax=5.5,ymin=ymin,ymax=ymax),
            inherit.aes=FALSE,fill=NA,colour=paper_style$text,linewidth=.33)+
  geom_text(data=focal,aes(label=copies,colour=I(text_colour)),size=2.2*scale,family=paper_style$font,fontface='bold')+
  geom_segment(data=groups,aes(x=12.66,xend=12.66,y=ymin,yend=ymax,colour=I(colour)),inherit.aes=FALSE,linewidth=.55)+
  geom_text(data=groups,aes(x=13.02,y=y,label=label,colour=I(colour)),inherit.aes=FALSE,size=2.15*scale,
            hjust=0,lineheight=1.1,family=paper_style$font,fontface='bold')+
  scale_fill_gradient2(low=pal['contraction'],mid='#F5F6F7',high=pal['expansion'],limits=c(-2.5,2.5),midpoint=0,
    name='Within-family z score',breaks=c(-2,0,2),guide=guide_colourbar(display='rectangles',nbin=48,barwidth=grid::unit(30,'mm'),barheight=grid::unit(2.2,'mm'),title.position='left'))+
  scale_x_continuous(breaks=1:12,labels=species_labels,limits=c(.5,15.25),expand=c(0,0))+
  scale_y_continuous(breaks=f$y,labels=paste(f$short_label,f$family_id,sep='  '),limits=c(.4,20),expand=c(0,0))+
  labs(tag='a',title='Family-copy profiles across Scenedesmaceae',
       subtitle='FACHB-416 copy numbers are labelled.',x=NULL,y=NULL)+
  theme3(scale)+theme(axis.line=element_blank(),axis.ticks=element_blank(),axis.text.x=element_text(angle=42,hjust=1,vjust=1,size=6.15*scale),
   axis.text.y=element_text(size=6.0*scale),legend.position='bottom',legend.text=element_text(size=6.0*scale),
   legend.title=element_text(size=6.0*scale),legend.margin=margin(0,0,0,0))
}

plot_overall<-function(scale=1){
 d<-data.frame(annotation=rep(c('Original','Annotation\nrescue'),each=2),direction=rep(c('Contracted','Expanded'),2),
    n=c(156,49,128,46))
 d$annotation<-factor(d$annotation,levels=c('Original','Annotation\nrescue'));d$direction<-factor(d$direction,levels=c('Contracted','Expanded'))
 ggplot(d,aes(annotation,n,fill=direction))+geom_col(width=.62,position=position_dodge(width=.68))+
  geom_text(aes(label=n),position=position_dodge(width=.68),vjust=-.55,size=2.7*scale,family=paper_style$font,fontface='bold')+
  scale_fill_manual(values=c(Contracted=unname(pal['contraction']),Expanded=unname(pal['expansion'])),name=NULL)+
  scale_y_continuous(limits=c(0,180),breaks=c(0,50,100,150),expand=c(0,0))+
  scale_x_discrete(expand=expansion(add=.55))+
  labs(tag='b',title='Contraction-dominated\nremodeling',subtitle=NULL,x=NULL,y='Number of families')+
  theme3(scale)+theme(legend.position='bottom',legend.text=element_text(size=6.0*scale),
    legend.key.size=grid::unit(2.4,'mm'),legend.spacing.x=grid::unit(.8,'mm'),
    axis.text.x=element_text(size=6.8*scale),axis.title.y=element_text(size=6.8*scale),plot.title=element_text(size=8.1*scale,face='bold',margin=margin(l=13*scale,b=3)))
}

selected_ids<-c('OG0000007','OG0000195','OG0000052','OG0000170','OG0000274','OG0006706','OG0000288')
cmp<-f[match(selected_ids,f$family_id),];cmp$yy<-7:1
cmp$colour<-ifelse(cmp$family_id=='OG0000288',pal['candidate'],ifelse(cmp$FCHB416<cmp$other11_median,pal['contraction'],pal['expansion']))
cmp$short_label[cmp$family_id=='OG0000007']<-'Protein kinase'
cmp$short_label[cmp$family_id=='OG0000288']<-'S0B precursor transport'
write.table(cmp,file.path(output_dir,'representative_copy_comparisons.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
plot_copies<-function(scale=1){
 ggplot(cmp,aes(y=yy))+
  annotate('rect',xmin=-Inf,xmax=Inf,ymin=.52,ymax=1.48,fill='#EAF4EF',colour=NA)+
  geom_segment(aes(x=other11_median,xend=FCHB416,yend=yy),colour='#A9B1BA',linewidth=.65)+
  geom_point(aes(x=other11_median),shape=21,fill='white',colour='#65717C',size=2.45*scale,stroke=.5)+
  geom_point(aes(x=FCHB416,fill=I(colour)),shape=21,colour='white',size=2.65*scale,stroke=.3)+
  geom_text(aes(x=18.0,label=paste0(FCHB416,' / ',other11_median),colour=I(colour)),hjust=.5,
    size=2.3*scale,fontface='bold',family=paper_style$font)+
  annotate('text',x=18,y=7.64,label='FACHB / ref.',size=1.95*scale,colour=paper_style$secondary,family=paper_style$font)+
  scale_x_continuous(limits=c(-.5,20.4),breaks=c(0,5,10,15),expand=c(0,0))+
  scale_y_continuous(breaks=cmp$yy,labels=paste(cmp$short_label,cmp$family_id,sep='\n'),limits=c(.5,8.0),expand=c(0,0))+
  labs(tag='c',title='Selected copy contrasts',subtitle=NULL,x='Gene copies',y=NULL)+
  theme3(scale)+theme(axis.line.y=element_blank(),axis.ticks.y=element_blank(),axis.text.y=element_text(size=6.25*scale,lineheight=1.03),
    axis.title.x=element_text(size=6.8*scale))
}
