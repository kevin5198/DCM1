#! /usr/bin/Rscript

library(getopt)
opt = getopt(matrix(c(
'dir','d',1,'character',
'group','u',1,'character',
'comp','c',1,'character',
'type','t',1,'character',
'clcol','l',2,'logical',
'clrow','r',2,'logical',
'colum','m',2,'logical',
'name','n',1,'integer',
'rowname','w',2,'logical',
'color','s',1,'character',
'size','z',1,'double',
'out','o',1,'character',
'help','h',0,'logical'
),byrow=TRUE, ncol=4));

usage<-function(){
	cat("This script is used to plot heatmap of differential expression.
Usage   Rscript",get_Rscript_filename(),"[options]
Options:
	-d, --dir	the file containing differential expression matrix or directory contain all compare of differential expression files(txt format).
	-u, --group	group information file.
	-c, --comp	compare information file.
	-t, --type	Gene/Transcript Type: protein_coding,lncRNA,sncRNA,...
	-l, --clcol	boolean values determining if columns should be clustered(default TRUE).
	-r, --clrow	boolean values determining if rows should be clustered(default TRUE).
	-m, --colum	scaled by colum(default FALSE).
	-n, --name	the column number of row names (default NULL, the names must be unique).
	-w, --rowname	show row names or not (default do not show).
	-s, --color	the colors of heatmap (default green,black,red).
	-z, --size	the gene names size in plot (default 5).
	-o, --out	uotput directory(default ./).
	-h, --help	display this help and exit.
	\n",sep=" ")
	q(status=1);
}

if (is.null(opt$dir) || is.null(opt$group) || !is.null(opt$help)){
	stop(usage())
}

if (is.null(opt$clcol)) {opt$clcol=TRUE}
if (is.null(opt$clrow)) {opt$clrow=TRUE}
if (is.null(opt$rowname)) {opt$rowname=FALSE}
if (is.null(opt$name)) {opt$name=NULL}
if (is.null(opt$colum)) {opt$colum=FALSE}
if (is.null(opt$out)) {opt$out="./"}
if (is.null(opt$color)) {opt$color="green,black,red"}
if (is.null(opt$size)) {opt$size=5}

library(pheatmap)
library(gplots)
colors=strsplit(opt$color,",")[[1]]
group<-read.table(opt$group,header=T,sep="\t",check.names=F,colClasses = "character")
gr<-unlist(strsplit(as.character(group$Group[1]),","))
ga<-data.frame(rep(group$Sample[1],length(gr)),gr)
names(ga)<-names(group)
if(nrow(group)>=2){
	for(i in 2:nrow(group)){
		gr<-unlist(strsplit(as.character(group$Group[i]),","))
		g<-data.frame(rep(group$Sample[i],length(gr)),gr)
		names(g)<-names(group)
		ga<-rbind(ga,g)
	}
}
	
if (dir.exists(opt$dir)){
	if(is.null(opt$comp)) {
		stop("Can not find compare information file.")
	}

##########################
#main programe
	comp<-read.table(opt$comp,header=T,sep="\t",check.names=F,colClasses = c(rep("character", 2), rep("numeric", 3)))
	library(pheatmap)
	library(gplots)

	for(i in 1:nrow(comp)){
		values<-c("p_value","q_value")
		value<-min(c(comp$p_value[i],comp$q_value[i]))
		pq<-values[which(c(comp$p_value[i],comp$q_value[i])==value)]
		file<-paste(opt$dir,"/",comp$Test[i],"_vs_",comp$Control[i],".txt",sep="")
	
		if(!file.exists(file)){
			next
		}
		dife<-read.table(file,header=T,sep="\t",row.names=opt$name,check.names=F)

		if(!is.null(opt$type)){
			lncRNA<-c("processed_transcript","lincRNA","3prime_overlapping_ncrna","antisense","non_coding","sense_intronic","sense_overlapping","TEC","known_ncrna","macro_lncRNA","bidirectional_promoter_lncrna")
			sncRNA<-c("snRNA","snoRNA","rRNA","Mt_tRNA","Mt_rRNA","misc_RNA","miRNA","ribozyme","sRNA","scaRNA","vaultRNA")
		#	dife<-read.table(file,header=T,sep="\t",row.names=1,check.names=F)
			otp<-c("Gene_Type","Trans_Type")
			tm<-otp[which(otp %in% names(dife))]
			if(length(tm)==0){
				tm<-"Gene_Type"
			}
			if(opt$type == "all"){
				type<-type<-unique(dife[,eval(tm)])
			}else{
				type<-unlist(strsplit(opt$type,","))
				if("lncRNA" %in% type){
					type<-c(type,lncRNA)
				}
				if("sncRNA" %in% type){
					type<-c(type,sncRNA)
				}
			}
	
			if(any(type %in% dife[,eval(tm)])){
				dif<-dife[which(dife[,eval(tm)] %in% type),]
			}else{
				stop("The type provided doesn't exist.\n")
			}
		}else{
			dif<-dife
		}

		vs<-as.character(t(comp[i,1:2]))
		g<-ga[which(ga$Group %in% vs),]
		anno_col<-data.frame(Group=g[,2])

		if(all(c(dif$p_value,dif$q_value)==1) || all(c(comp$p_value[i],comp$q_value[i])==1)){
			up<-dif[which(dif$Fold_Change >=comp$fc[i]),]
			down<-dif[which(dif$Fold_Change <=1/comp$fc[i]),]
			updown<-rbind(up,down)
			heatdf<-updown[,which(colnames(updown) %in% g[,1])]
			#heatdf<-updown[,11:ncol(updown)]
		}else{
			sig<-dif[which(dif[,which(names(dif)==pq)]<=value),]
			up<-sig[which(sig$Fold_Change >=comp$fc[i]),]
			down<-sig[which(sig$Fold_Change <=1/comp$fc[i]),]
			updown<-rbind(up,down)
			heatdf<-updown[,which(colnames(updown) %in% g[,1])]
			#heatdf<-updown[,11:ncol(updown)]
		}

		#vs<-as.character(t(comp[i,1:2]))
		#g<-ga[which(ga$Group %in% vs),]
		#anno_col<-data.frame(Group=g[,2])
		#rownames(anno_col)<-g[,1]

		if(length(which(anno_col[,1]==vs[1]))>1 && length(which(anno_col[,1]==vs[2]))>1){
			sc<-"row"
		}else if(opt$colum){
			sc<-"column"
		}else{
			sc<-"none"
			heatdf<-log2(heatdf+1)
		}

		if(any(duplicated(g[,1]))){
			anno_col<-NA
		}else{
			rownames(anno_col)<-g[,1]
		}
		if(nrow(heatdf)>1){
			pheatmap(heatdf,show_rownames=opt$rowname,annotation_col=anno_col,scale=sc,annotation_legend=F,color = colorpanel(128, colors[1], colors[2], colors[3]),border_color = NA,cluster_cols=opt$clcol,cluster_rows=opt$clrow,filename=paste(opt$out,"/heatmap_",vs[1],"_vs_",vs[2],".pdf",sep=""),width=9,height=12,fontsize_row=opt$size)
			pheatmap(heatdf,show_rownames=opt$rowname,annotation_col=anno_col,scale=sc,annotation_legend=F,color = colorpanel(128, colors[1], colors[2], colors[3]),border_color = NA,cluster_cols=opt$clcol,cluster_rows=opt$clrow,filename=paste(opt$out,"/heatmap_",vs[1],"_vs_",vs[2],".png",sep=""),width=9,height=12,fontsize_row=opt$size)
		}else{
			cat('no enough data',file=paste(opt$out,"/heatmap_",vs[1],"_vs_",vs[2],".log",sep=""))
		}
	}
}else if(file.exists(opt$dir)){
	f<-read.table(opt$dir,header=T,sep="\t",check.names=F,row.names=opt$name)
	mfs<-f[,which(names(f) %in% ga$Sample)]
	mm<-apply(mfs,1,sum)
	mfs<-mfs[mm>0,]
	anno_col<-data.frame(Group=ga$Group)
	rownames(anno_col)<-ga$Sample

	if(all(table(anno_col)>1)){
		sc<-"row"
	}else if(opt$colum){
		sc<-"column"
	}else{
		sc<-"none"
		mfs<-log2(mfs+1)
	}

	pheatmap(mfs,show_rownames=opt$rowname,annotation_col=anno_col,scale=sc,annotation_legend=F,color = colorpanel(128, colors[1], colors[2], colors[3]),border_color = NA,cluster_cols=opt$clcol,cluster_rows=opt$clrow,filename=paste(opt$out,"/heatmap.pdf",sep=""),width=9,height=12,fontsize_row=opt$size)
	pheatmap(mfs,show_rownames=opt$rowname,annotation_col=anno_col,scale=sc,annotation_legend=F,color = colorpanel(128, colors[1], colors[2], colors[3]),border_color = NA,cluster_cols=opt$clcol,cluster_rows=opt$clrow,filename=paste(opt$out,"/heatmap.png",sep=""),width=9,height=12,fontsize_row=opt$size)
}

file.remove("Rplots.pdf")

