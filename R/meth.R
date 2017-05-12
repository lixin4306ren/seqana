#' plot and calculate methylation along sliding windows from bsseq object
#' this function plot and calculate slidng windows methylation
#' @param BS.BM bsseq object
#' @param region GRange object
#' @import GenomicRanges data.table bsseq
#' @export
#' @examples
#' plot_meth_sliding_bsseq()
plot_meth_sliding_bsseq<-function(BS.BM,region,win.base=5000,step.base=5000,up.base=200000,down.base=200000,win.percent=0.01,step.percent=0.01,up.percent=1,down.percent=1,min.num=20,min.mc.depth=0,frame=TRUE,main="DNAm",col="black",append=FALSE,dist="base",x.axis=TRUE,ylim=c(0,1)){
small.num=0.00001
	f<-function(x,a,b,up,down){
		width(x)*a->win
		width(x)*b->step
		region.start<-start(x)-width(x)*up
		region.end<-end(x)+width(x)*down
		round(seq(region.start,region.end-win+1+small.num,step))->start
		round(seq(region.start+win-1,region.end+small.num,step))->end
		GRanges(seqnames=seqnames(x),range=IRanges(start=start,end=end))
	}

rowSums(getCoverage(BS.BM,type="Cov"))->Cov
rowSums(getCoverage(BS.BM,type="M"))->M

Cov[Cov > 0]->cg.x.cov
M[Cov>0]->cg.x.m
M[M > min.mc.depth]->mcg.x

cg1<- BS.BM@rowRanges[Cov>0]
mcg1<- BS.BM@rowRanges[M>min.mc.depth]

if(dist=="percent"){
win<-f(region[1],win.percent,step.percent,up.percent,down.percent)
for(i in 2:length(region)){
	tmp.win<-f(region[i],win.percent,step.percent,up.percent,down.percent)
	c(win,tmp.win)->win
}
}
else if(dist=="base"){
win<-generate.sliding.window(region,win.percent,step.base,up.base,down.base)
}
print("starting")

as.data.frame(findOverlaps(win,cg1))->tmp1
countOverlaps(win,mcg1)->mcg
countOverlaps(win,cg1)->c

which(countOverlaps(win,cg1)==0)->tmp2 #regions without Cs

data.table(id=as.character(tmp1$queryHits),cov=cg.x.cov[tmp1$subjectHits],m=cg.x.m[tmp1$subjectHits])->new1
setkey(new1,id)
as.integer(new1[,sum(cov),by=id]$id)->id1
data.table(id=new1[,sum(cov),by=id]$id,start=start(win)[id1],end=end(win)[id1],m=new1[,sum(m),by=id]$V1,cov=new1[,sum(cov),by=id]$V1)->tmp1
tmp1<-rbind(data.frame(tmp1),data.frame(id=as.character(tmp2),start=start(win)[tmp2],end=end(win)[tmp2],m=rep(0,length(tmp2)),cov=rep(0,length(tmp2))))
sort(as.integer(tmp1$id),index.return = T)->ind1
tmp1[ind1$ix,]->tmp1

tmp1$mcg=mcg
tmp1$c=c
tmp1$cov<min.num ->ind
tmp1$cov[ind]<-0
tmp1$m[ind]<-0
tmp1$meth=tmp1$m/tmp1$cov
matrix(tmp1$meth,nrow=length(region))->tbls

if(dist=="percent"){
	pos1<-1
	pos2<-up.percent/step.percent+1
	pos3<-ncol(tbls)-down.percent/step.percent
	pos4<-ncol(tbls)
axis=c(-up.percent,"0 %","100 %",down.percent)
}
else if(dist=="base"){
	pos1<-1
	pos2<-(ncol(tbls)-1/step.percent)/2+1
	pos3<-ncol(tbls)-(ncol(tbls)-1/step.percent)/2
	pos4<-ncol(tbls)

axis=c(-up.base,"0 %","100 %",down.base)
}

if (append==FALSE){
	if(frame==TRUE){
	plot(colMeans(tbls,na.rm=T),ylim=ylim,ylab="methylation level",xlab="",xaxt="n",type="l",main=main,col=col)
	}else{
	plot(colMeans(tbls,na.rm=T),ylim=ylim,ylab="methylation level",xlab="",xaxt="n",yaxt="n",type="l",main=main,col=col,axes=F)
	}
	axis(2)
}else{
	points(colMeans(tbls,na.rm=T),type="l",col=col)

}
if(x.axis==TRUE){axis(1,c(pos1,pos2,pos3,pos4),axis)}
list(tbls=tbls,x.axis.pos=c(pos1,pos2,pos3,pos4),x.axis=axis)->plot.out
return(plot.out)
}

#' calculate methylation level based on given Grange from bsseq object
#' @import GenomicRanges data.table bsseq
#' @export
#' @examples
#' meth_count_region_for_bsseq()


meth_count_region_for_bsseq<-function(BS,range,min.cov=5,min.cpg.num=1){
win<-range
cpg.range<-BS@rowRanges

f.tmp<-function(BS.tmp,win,min.cov,min.cpg.num){
if(is.null(names(win))){names(win)<-1:length(win)}
cov<-getCoverage(BS.tmp,type="Cov")
m<-getCoverage(BS.tmp,type="M")
as.data.frame(findOverlaps(win,cpg.range))->tmp1
which(countOverlaps(win,cpg.range)==0)->tmp2 #regions without sequenced Cs
data.table(id=names(win)[tmp1$queryHits],cov=cov[tmp1$subjectHits],m=m[tmp1$subjectHits])->new1
setkey(new1,id)
new1[,sum(cov),by=id]$id->id1
data.table(id=new1[,sum(cov),by=id]$id,chr=as.vector(seqnames(win[id1])),start=start(win[id1]),end=end(win[id1]),m=new1[,sum(m),by=id]$V1,cov=new1[,sum(cov),by=id]$V1)->tmp1
tmp1<-rbind(data.frame(tmp1),data.frame(id=names(win)[tmp2],chr=as.vector(seqnames(win[tmp2])),start=start(win)[tmp2],end=end(win)[tmp2],m=rep(0,length(tmp2)),cov=rep(0,length(tmp2))))
order(tmp1$chr,tmp1$start)->ind1
tmp1[ind1,]->tmp1
tmp1.range=GRanges(seqnames=tmp1$chr,range=IRanges(start=tmp1$start,end=tmp1$end))
countOverlaps(tmp1.range,cpg.range)->total.c
countOverlaps(tmp1.range,cpg.range)->covered.c
tmp1$total.c<-total.c
tmp1$covered.c<-covered.c

tmp1$cov<min.cov | tmp1$covered.c<min.cpg.num ->ind
tmp1$m_level[ind]<-NA
tmp1$m_level[!ind]<-tmp1$m[!ind]/tmp1$cov[!ind]
return(tmp1)
}

out<-f.tmp(BS[,1],win,min.cov,min.cpg.num)[,c(2,3,4,9)]
print (1)
if(length(sampleNames(BS))>1){    
for(i in 2:length(sampleNames(BS))){
print(i)
BS.tmp<-BS[,i]
f.tmp(BS.tmp,win,min.cov,min.cpg.num)->out.tmp
cbind(out,out.tmp[,9])->out
}
}
return(out)

}
