#' generate sliding window ranges
#' @param region GRange object
#' @import GenomicRanges
#' @export
#' @examples
#' generate.sliding.window()
generate.sliding.window<-function(region,step.percent=0.01,step.base=10,up=1000,down=1000){
chr<-rep(seqnames(region),1/step.percent+up/step.base+down/step.base)
start<-start(region)
end<-end(region)
generate.sliding.window.start(start,end,step.percent,step.base,up,down)->win.start
generate.sliding.window.end(start,end,step.percent,step.base,up,down)->win.end
GRanges(seqnames=chr,range=IRanges(start=c(win.start),end=c(win.end)))->win
win
}

generate.sliding.window.start<-function(start,end,step.percent=0.01,step.base=10,up=1000,down=1000){
library(GenomicRanges)
if(up==0){
tbls.left<-NULL
}
else{
out.left<-NULL
for(i in 1:(up/step.base)){
round((start-up)+(step.base*(i-1)))->tmp
c(out.left,tmp)->out.left
}
matrix(out.left,nrow=length(start))->tbls.left
}

out.middle<-NULL

for(i in 1:(1/step.percent)){
round(start+((end-start+1)*step.percent*(i-1)))->tmp
c(out.middle,tmp)->out.middle
}
matrix(out.middle,nrow=length(start))->tbls.middle

if(down==0){
tbls.right<-NULL
}
else{
out.right<-NULL
for(i in 1:(down/step.base)){
round(end+(step.base*(i-1)))->tmp
c(out.right,tmp)->out.right
}
matrix(out.right,nrow=length(start))->tbls.right
}

cbind(tbls.left,tbls.middle,tbls.right)->tbls
tbls
}


generate.sliding.window.end<-function(start,end,step.percent=0.01,step.base=10,up=1000,down=1000){
library(GenomicRanges)

if(up==0){
tbls.left<-NULL
}
else{
out.left<-NULL
for(i in 1:(up/step.base)){
round((start-up-1)+(step.base*(i)))->tmp
c(out.left,tmp)->out.left
}
matrix(out.left,nrow=length(start))->tbls.left
}

out.middle<-NULL

for(i in 1:(1/step.percent)){
round(start-1+((end-start+1)*step.percent*(i)))->tmp
c(out.middle,tmp)->out.middle
}
matrix(out.middle,nrow=length(start))->tbls.middle

if(down==0){
tbls.right<-NULL
}
else{
out.right<-NULL
for(i in 1:(down/step.base)){
round(end-1+(step.base*(i)))->tmp
c(out.right,tmp)->out.right
}
matrix(out.right,nrow=length(start))->tbls.right
}

cbind(tbls.left,tbls.middle,tbls.right)->tbls
tbls
}


#' calculate overlap % between 2 bed files, need bedtools installed
#' @param bed1 first bed object
#' @param bed2 second bed object
#' @export
#' @examples
#' bedtools.coveragebed.range()
bedtools.coveragebed.range<-function(functionstring="coverageBed",bed1,bed2,opt.string="",mcol=TRUE)
{
  #create temp files
  a.file=tempfile()
  b.file=tempfile()
  out = tempfile()
  options(scipen =99) # not to use scientific notation when writing out
 data.frame(chr=as.character(seqnames(bed1)),start=start(bed1),end=end(bed1))->a.out
 if(mcol==TRUE){
 data.frame(chr=as.character(seqnames(bed2)),start=start(bed2),end=end(bed2),gene=mcols(bed2)$gene)->b.out
 }
 else{
 data.frame(chr=as.character(seqnames(bed2)),start=start(bed2),end=end(bed2))->b.out
 }

 write.table(a.out,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
 write.table(b.out,file=b.file,quote=F,sep="\t",col.names=F,row.names=F)

  # create the command string and call the command using system()
  command=paste(functionstring,"-a",a.file,"-b",b.file,opt.string,">",out,sep=" ")
  cat(command,"\n")
  try(system(command))

  res=read.table(out,header=F)
  return(res)
}



#' calculate overlap % between 2 bed files, need bedtools installed
#' @param bed1 first bed file
#' @param bed2 second bed file
#' @export
#' @examples
#' bedtools.coveragebed()

bedtools.coveragebed<-function(functionstring="coverageBed",bed1,bed2,opt.string="")
{
  #create temp files
  #a.file=tempfile()
  #b.file=tempfile()
  out   =tempfile()
  options(scipen =99) # not to use scientific notation when writing out
 
  #write.table(bed1,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
  #write.table(bed2,file=b.file,quote=F,sep="\t",col.names=F,row.names=F)
 
  # create the command string and call the command using system()
  command=paste(functionstring,"-a",bed1,"-b",bed2,opt.string,">",out,sep=" ")
  cat(command,"\n")
  try(system(command))
 
  res=read.table(out,header=F)
  #unlink(a.file);
  #unlink(b.file);
  unlink(out)
  return(res)
}

