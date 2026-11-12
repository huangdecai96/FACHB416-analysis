# Final Figure2 export from the author-approved composition.
# The phylogeny/synteny composition was finalized manually; no alternative redraw is substituted.
aa<-grep('^--file=',commandArgs(FALSE),value=TRUE)
source(file.path(dirname(sub('^--file=','',aa[1])),'helpers/setup.R'))
im<-png::readPNG(file.path(DATA,'Figure2/approved_slide.png'))
keep<-im[,,1]<1|im[,,2]<1|im[,,3]<1
rr<-which(rowSums(keep)>0);cc<-which(colSums(keep)>0);pad<-18L
top<-max(1L,min(rr)-pad);bottom<-min(dim(im)[1],max(rr)+pad)
left<-max(1L,min(cc)-pad);right<-min(dim(im)[2],max(cc)+pad)
cropped<-im[top:bottom,left:right,,drop=FALSE]
stopifnot(all(rr>=top&rr<=bottom),all(cc>=left&cc<=right))
png::writePNG(cropped,file.path(OUT,'Figure2.png'))
check<-png::readPNG(file.path(OUT,'Figure2.png'))
stopifnot(identical(dim(cropped),dim(check)),all(cropped==check))
write.table(data.frame(left=left,top=top,right=right,bottom=bottom,all_nonwhite_pixels_preserved=TRUE),file.path(ROOT,'verification/Figure2_crop.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
capture.output(sessionInfo(),file=file.path(ROOT,'verification/Figure2_sessionInfo.txt'))
