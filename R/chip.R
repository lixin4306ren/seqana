#' merge peaks based on replicates and find differentially peaks between 2 samples
#' 
#' @param peak1 sample1 rep1
#' @param peak2	sample1 rep2
#' @param peak3 sample2 rep1
#' @param peak4 sample2 rep2
#' @import GenomicRanges
#' @export
#' @examples
#' diff_peaks_replicate()

diff_peaks_replicate<-function(peak1=NULL,peak2=NULL,peak3=NULL,peak4=NULL,min.dis=10000,min.dis.for.diff=10000,remove.small0=0,remove.small0.for.diff=0,remove.small.for.diff=NULL,remove.small=NULL,remove.small.each.replicate=NULL,remove.small2=0,name1=NULL,name2=NULL,name3=NULL){
if(!is.null(remove.small)){merge_peaks(peak1,min.dis,remove.small=remove.small.each.replicate,return.range=T,remove.small0=remove.small0)->merge1}
else{merge_peaks(peak1,min.dis,return.range=T,remove.small0=remove.small0)->merge1}
if(!is.null(remove.small)){merge_peaks(peak2,min.dis,remove.small=remove.small.each.replicate,return.range=T,remove.small0=remove.small0)->merge2}
else{merge_peaks(peak2,min.dis,return.range=T,remove.small0=remove.small0)->merge2}
if(!is.null(remove.small)){merge_peaks(peak3,min.dis,remove.small=remove.small.each.replicate,return.range=T,remove.small0=remove.small0)->merge3}
else{merge_peaks(peak3,min.dis,return.range=T,remove.small0=remove.small0)->merge3}
if(!is.null(remove.small)){merge_peaks(peak4,min.dis,remove.small=remove.small.each.replicate,return.range=T,remove.small0=remove.small0)->merge4}
else{merge_peaks(peak4,min.dis,return.range=T,remove.small0=remove.small0)->merge4}

get_common_peak3(merge1,merge2,return.range=T)->merge.peak1
width(merge.peak1)>=remove.small ->ind1
merge.peak1[ind1]->merge.peak1
data.frame(chr=as.vector(seqnames(merge.peak1)),start=start(merge.peak1),end=end(merge.peak1))->tmp
file.name<-paste(name1,name3,"large.bed",sep=".")
write.table(tmp,file=file.name,quote=F,col.names=F,row.names=F,sep="\t")
get_common_peak3(merge3,merge4,return.range=T)->merge.peak2
width(merge.peak2)>=remove.small ->ind2
merge.peak2[ind2]->merge.peak2
data.frame(chr=as.vector(seqnames(merge.peak2)),start=start(merge.peak2),end=end(merge.peak2))->tmp
file.name<-paste(name2,name3,"large.bed",sep=".")
write.table(tmp,file=file.name,quote=F,col.names=F,row.names=F,sep="\t")


merge_peaks(peak1,min.dis=min.dis.for.diff,return.range=T,remove.small0=remove.small0.for.diff)->merge1.for.diff
merge_peaks(peak2,min.dis=min.dis.for.diff,return.range=T,remove.small0=remove.small0.for.diff)->merge2.for.diff
merge_peaks(peak3,min.dis=min.dis.for.diff,return.range=T,remove.small0=remove.small0.for.diff)->merge3.for.diff
merge_peaks(peak4,min.dis=min.dis.for.diff,return.range=T,remove.small0=remove.small0.for.diff)->merge4.for.diff

get_common_peak3(merge1.for.diff,merge2.for.diff,return.range=T)->merge.peak1.for.diff
width(merge.peak1.for.diff)>=remove.small.for.diff ->ind1
merge.peak1.for.diff[ind1]->merge.peak1.for.diff

get_common_peak3(merge3.for.diff,merge4.for.diff,return.range=T)->merge.peak2.for.diff
width(merge.peak2.for.diff)>=remove.small.for.diff ->ind2
merge.peak2.for.diff[ind2]->merge.peak2.for.diff


######output all domain
data.frame(chr=as.vector(seqnames(merge.peak1.for.diff)),start=start(merge.peak1.for.diff),end=end(merge.peak1.for.diff))->tmp
file.name<-paste(name1,name3,"all.bed",sep=".")
write.table(tmp,file=file.name,quote=F,col.names=F,row.names=F,sep="\t")
data.frame(chr=as.vector(seqnames(merge.peak2.for.diff)),start=start(merge.peak2.for.diff),end=end(merge.peak2.for.diff))->tmp
file.name<-paste(name2,name3,"all.bed",sep=".")
write.table(tmp,file=file.name,quote=F,col.names=F,row.names=F,sep="\t")


setdiff(merge.peak1,reduce(c(merge.peak2,merge.peak2.for.diff)))->peak.up
setdiff(merge.peak2,reduce(c(merge.peak1,merge.peak1.for.diff)))->peak.down
#setdiff(merge.peak1,merge.peak2)->peak.up
#setdiff(merge.peak2,merge.peak1)->peak.down

width(peak.up)>=remove.small2 ->ind1
width(peak.down)>=remove.small2 ->ind2
peak.up[ind1]->peak.up
peak.down[ind2]->peak.down
data.frame(chr=as.vector(seqnames(peak.up)),start=start(peak.up),end=end(peak.up),direction="Up")->out1
data.frame(chr=as.vector(seqnames(peak.down)),start=start(peak.down),end=end(peak.down),direction="Down")->out2
rbind(out1,out2)->out
file.name<-paste(name1,name2,name3,"diff.bed",sep=".")
write.table(out,file=file.name,quote=F,col.names=F,row.names=F,sep="\t")

}

#' merge neighbouring peaks to large one by given distance
#' @param peak peak bed file path
#' @param mis.dis distance cutoff, default 10000
#' @import GenomicRanges
#' @export
#' @examples
#' merge_peaks()
merge_peaks<-function(peak=peakfile,min.dis=10000,remove.small=NULL,remove.small0=0,return.range=FALSE){
######### first remove region less than remove.small0
read.table(peak)->x
GRanges(seqnames=x[,1],range=IRanges(start=x[,2],end=x[,3]))->y
width(y)>=remove.small0 ->ind
y[ind]->y
gaps(y)->gap
width(gap)< min.dis ->ind
reduce(c(y,gap[ind]))->z

######### after merge by mis.dis, remove merged region less than remove.small
if(!is.null(remove.small)){
width(z)>=remove.small -> ind
z[ind]->z
}
if(return.range==FALSE){
data.frame(chr=as.vector(seqnames(z)),start=start(z),end=end(z))->out
return(out)
}
else{
return(z)
}

}


#' extract common peaks based on replicates
#' @param read1 rep1
#' @param read1 rep2
#' @param min.overlap overlap % cutoff between two files, default 0.5
#' @import GenomicRanges
#' @export
#' @examples
#' get_common_peak3()
get_common_peak3<-function(peak1,peak2,return.range=FALSE,min.overlap=0.5){
bedtools.coveragebed.range(bed1=peak2,bed2=peak1,mcol=F)->cov1
bedtools.coveragebed.range(bed1=peak1,bed2=peak2,mcol=F)->cov2

cov1[,ncol(cov1)]>min.overlap ->ind1
cov2[,ncol(cov2)]>min.overlap ->ind2
cov1[ind1,c(1,2,3)]->peak1
cov2[ind2,c(1,2,3)]->peak2

GRanges(seqname=peak1$V1,range=IRanges(start=peak1$V2,end=peak1$V3))->peak1.range
GRanges(seqname=peak2$V1,range=IRanges(start=peak2$V2,end=peak2$V3))->peak2.range

peak1.range %over% peak2.range ->ind1
peak2.range %over% peak1.range ->ind2
peak1.range[ind1]->peak1.range
peak2.range[ind2]->peak2.range

c(peak1.range,peak2.range)-> both.peak.range
reduce(both.peak.range)-> both.peak.range
if(return.range==FALSE){
data.frame(chr=as.vector(seqnames(both.peak.range)),start=start(both.peak.range),end=end(both.peak.range))->both.peak
return(both.peak)
}
else{
return(both.peak.range)
}
}

#' this function calculate RPKM for chip-seq data, support bam and bed input format
#' @param read input file, could be bam or bed format
#' @param range GRange object
#' @import GenomicRanges Repitools
#' @export
#' @examples
#' cal.rpkm.chip()
cal.rpkm.chip<-function(read,range,input="bam",pair=F){
if(input=="bam"){
       	if(pair==F){BAM2GRanges(read)->chip}
       	else if(pair==T){BAMpair2GRanges(read)->chip}
}
else if(input=="bed"){
read.table(read)->tmp
GRanges(seqname=tmp$V1,range=IRanges(start=tmp$V2,end=tmp$V3))->chip
}

countOverlaps(range,chip)->tmp
tmp/(length(chip)/1000000)/(width(range)/1000)->rpkm
mcols(range)$rpkm<-rpkm
return(range)
}


#' peak annotate function
#' this function find nearest gene asscosiated with given peaks, deault database is human.gencodeV19
#' @param peaks peak file, could be bed file path or Grange object, default bed fiel path
#' @import bumphunter hash GenomicFeatures
#' @export
#' @examples
#' annotatePeaks()
annotatePeaks<-function(peaks,genedb="~/amber3/no_back_up/data/R_data/human.gencodeV19",peak.type="BED",by_gene=TRUE,min.distance=2000,pro.len=500,pro=TRUE,subchr=FALSE,head=FALSE){

if(peak.type=="BED"){
	if(head==TRUE){
		read.table(peaks,header=T)->peak
	}
	else{
		read.table(peaks)->peak
	}
	if(subchr==TRUE){
		GRanges(seqname=sub("chr","",peak[,1]),range=IRanges(start=peak[,2],end=peak[,3]))->peak.range
	}
	else{
		GRanges(seqname=peak[,1],range=IRanges(start=peak[,2],end=peak[,3]))->peak.range
	}
}
else if(peak.type=="RANGE"){
peaks->peak.range
}


names(peak.range)<-1:length(peak.range)

loadDb(genedb)->db
if(by_gene==TRUE){
genes(db)->gene
}
else{
transcriptsBy(db,by="gene")->gene
unlist(gene)->gene
}

annotateNearest(peak.range,gene)->anno
anno$peak_id<-as.numeric(names(peak.range))
abs(anno$dist)<=min.distance -> ind
anno[ind,]->anno
names(gene)[anno$subjectHits]->anno$gene

if(pro==TRUE){
if(by_gene==TRUE){
promoters(gene, upstream=pro.len, downstream=pro.len)->pro
}
else{
promoters(db, upstream=pro.len, downstream=pro.len)->pro
}
if(by_gene==TRUE){
hash(keys=mcols(gene)$gene_id,values=names(gene))->table
values(table,mcols(pro)$gene_id)->names(pro)
}
else{
hash(keys=mcols(gene)$tx_name,values=names(gene))->table
values(table,mcols(pro)$tx_name)->names(pro)
}

annotateNearest(peak.range,pro)->anno.pro
anno.pro$peak_id<-as.numeric(names(peak.range))
abs(anno.pro$dist)==0 & !is.na(anno.pro$dist)-> ind
anno.pro[ind,]->anno.pro
names(pro)[anno.pro$subjectHits]->anno.pro$gene
list(all=anno,promoter=anno.pro)->out
}
else{
list(all=anno)->out
}
return(out)
}


