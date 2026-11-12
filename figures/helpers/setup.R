suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(grid)})
aa<-grep('^--file=',commandArgs(FALSE),value=TRUE)
SCRIPT_DIR<-normalizePath(dirname(sub('^--file=','',aa[1])),winslash='/')
args<-commandArgs(TRUE)
if(length(args)<2)stop('Usage: Rscript figures/FigureN.R DATA_DIRECTORY NEW_OUTPUT_DIRECTORY')
DATA<-normalizePath(args[1],winslash='/',mustWork=TRUE)
ROOT<-normalizePath(args[2],winslash='/',mustWork=FALSE)
if(ROOT==DATA || startsWith(ROOT,paste0(DATA,'/')) || startsWith(DATA,paste0(ROOT,'/')))stop('Output must be separate from input data')
REPO<-normalizePath(file.path(SCRIPT_DIR,'..'),winslash='/')
if(ROOT==REPO || startsWith(REPO,paste0(ROOT,'/')))stop('Output must not contain the code repository')
if(dir.exists(ROOT)&&length(list.files(ROOT,all.files=TRUE,no..=TRUE)))stop('Use a new empty output directory')
for(n in c('figures','verification','work'))dir.create(file.path(ROOT,n),recursive=TRUE,showWarnings=FALSE)
OUT<-file.path(ROOT,'figures');output_dir<-file.path(ROOT,'verification')
source(file.path(SCRIPT_DIR,'helpers/paper_style.R'))
ink<-'#20242B';rose<-'#A84F65';green<-'#1B927F';blue<-'#2F73AE';muted<-'#5A626D'
rd<-function(fig,name)read.delim(file.path(DATA,fig,name),check.names=FALSE,stringsAsFactors=FALSE,fileEncoding='UTF-8',quote='',comment.char='')
export_fig<-function(nm,fun,h,default_dpi=500){
 dpi<-as.numeric(Sys.getenv('FIGURE_DPI',as.character(default_dpi)))
 cairo_pdf(file.path(OUT,paste0(nm,'.pdf')),width=180/25.4,height=h/25.4,family='Arial');fun();dev.off()
 ragg::agg_png(file.path(OUT,paste0(nm,'.png')),width=180,height=h,units='mm',res=dpi,background='white');fun();dev.off()
 capture.output(sessionInfo(),file=file.path(ROOT,'verification',paste0(nm,'_sessionInfo.txt')))
}
base_theme<-function(){theme_classic(base_size=8.5,base_family='Arial')+theme(text=element_text(colour=ink),axis.text=element_text(size=8,colour=ink),axis.title=element_text(size=8.5),axis.title.y=element_text(margin=margin(r=2)),axis.title.x=element_text(margin=margin(t=3)),axis.line=element_blank(),axis.line.x=element_blank(),axis.line.y=element_blank(),axis.ticks=element_line(linewidth=.4,colour='black'),axis.ticks.length=unit(1,'mm'),panel.border=element_rect(fill=NA,colour='black',linewidth=.45),plot.margin=margin(3,3,3,3),plot.title=element_text(size=9,face='bold'),plot.subtitle=element_text(size=7.5),plot.tag=element_text(size=12,face='bold'),legend.title=element_text(size=7.5),legend.text=element_text(size=7.5),legend.key.size=unit(3,'mm'))}
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
page<-function(h){grid.newpage();pushViewport(viewport(width=unit(180,'mm'),height=unit(h,'mm'),xscale=c(0,180),yscale=c(h,0)))}
gtext<-function(s,x,y,size=6)grid.text(s,unit(x,'native'),unit(y,'native'),just='left',gp=gpar(fontfamily='Arial',fontsize=size,col=ink))
plot_at<-function(p,x,y,w,h){pushViewport(viewport(x=unit(x,'native'),y=unit(y,'native'),width=unit(w,'mm'),height=unit(h,'mm'),just=c('left','top')));grid.draw(ggplotGrob(p));popViewport()}

