require(ggplot2, quiet=TRUE)
require(getopt, quiet=TRUE)

opt = getopt(matrix(c(
'path','p',1,'character',
'comp','c',1,'character',
'out','o',1,'character',
'help','h',0,'logical'
),byrow=TRUE, ncol=4));

usage<-function(){
	cat("
Usage: Rscript Volcano_plot.r [-p ./] [-c ./compare.txt] [-o ./]
Options:
	-p, --path	data file path.
	-c, --comp	compare.txt file.
	-o, --out	prefix of out files.
	-h, --help	display this help and exit.
	\n",sep=" ")
	q(status=1);
}

if (!is.null(opt$help)){
	stop(usage())
}

if (is.null(opt$path)) {opt$path="./"}
if (is.null(opt$comp)) {opt$comp="./compare.txt"}
if (is.null(opt$out)) {opt$out="./"}

comp<-read.table(opt$comp,header= T,sep="\t",colClasses=c("character","character","numeric","numeric","numeric")) # 指定每列数据类型，某列全为T或F会被读为TRUE或FALSE

for(i in 1:nrow(comp)){
	d<-read.table(paste(opt$path,comp[i,1],"_vs_",comp[i,2],".txt",sep=""),header=T,sep="\t")
	sig<-d
	du<-sig$fold_change
	du[which(sig$P_value > comp[i,4])]<-1
	du[which(du >= comp[i,3])]<-1000
	du[which(du <= 1/comp[i,3])]<-"down"
	du[which(du == "1000")]<-"up"
	vol2<-data.frame(log2FC=log2(sig$fold_change),log10pval=-log10(sig$P_value),du=du)
	upl<-paste("up_regulated peaks (",length(grep("up",du)),")",sep="")
	downl<-paste("down_regulated peaks (",length(grep("down",du)),")",sep="")
	if(nrow(vol2)==0){
		break
	}
	p<-ggplot(vol2,aes(x=log2FC,y=log10pval,colour=du))
	p<-p+geom_point(shape=20,size=1)
	p<-p+scale_colour_manual(values=c("up"="brown3","down"="seagreen3"),breaks = c("up","down"))
	p<-p+geom_hline(yintercept = -log10(comp[i,4]),colour="grey35",linetype=2,size=0.3)
	p<-p+geom_vline(xintercept = c(-log2(comp[i,3]),log2(comp[i,3])),colour="grey35",linetype=2,size=0.3)
	p<-p+theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5),legend.position="none")
	p<-p+labs(x="log2(Fold Change)",y=paste("-log10(","Pvalue",")",sep=""),title=paste(comp[i,1],"_vs_",comp[i,2],sep=""))
	p<-p+theme(plot.margin=margin(15, 40, 10, 10))
	p<-p+scale_x_continuous(limits=c(-6,6),breaks=c(-6,-4,-2,0,2,4,6))
	p<-p+scale_y_continuous(limits=c(0,20),breaks=c(0,5,10,15,20))
	p<-p+annotate("text",x=c(-4,4),y=20,label=c("190 down regulated m6A peaks","106 up regulated m6A peaks"),color=c("seagreen3","brown3"))
	ggsave(paste(opt$out,comp[i,1],"_vs_",comp[i,2],"_volcano.pdf",sep=""),plot=p,device="pdf",width=7,height=7,units="in")
	ggsave(paste(opt$out,comp[i,1],"_vs_",comp[i,2],"_volcano.png",sep=""),plot=p,device="png",width=7,height=7,units="in",dpi=600)
}
