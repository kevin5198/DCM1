library("ggplot2",quietly=TRUE,warn.conflicts=FALSE)

FoldChange <- 1.5
pvalue <- 0.05
Test_name <- "Exp"
Control_name <- "Control"
type <- "genes"
color <- c("brown3","grey35","seagreen3")

volcano<-function(x,log2FC=1,value=2,level=0.05,times=1.5,title="Test vs Control",labe=NULL){
	du2<-x[,log2FC]
	du2[which(x[,value]>level)]<-0
	du2[which(du2>=log2(times))]<-10000
	du2[which(du2<=(-log2(times)))]<-"down"
	du2[which(du2=="10000")]<-"up"
	du2[which(!du2 %in% c("up","down"))]<-"nde"
	vol2<-data.frame(x[,log2FC],-log10(x[,value]),du2,id=x[,3])
	upl<-paste("up_regulated ",type," (",length(grep("up",du2)),")",sep="")
	downl<-paste("down_regulated ",type," (",length(grep("down",du2)),")",sep="")
	ndel<-paste("not differential expressed (",length(grep("nde",du2)),")",sep="")
	p<-ggplot(vol2,aes(x=vol2[,1],y=vol2[,2],colour=du2))
	p<-p+geom_point(shape=20,size=1)
	p<-p+scale_colour_manual(values=c("up"=color[1],"nde"=color[2],"down"=color[3]),breaks = c("up","nde","down"),labels=c(upl,ndel,downl))
	p<-p+geom_hline(yintercept = -log10(level),colour="green4",linetype=2,size=0.3)
	p<-p+geom_vline(xintercept = c(-log2(times),log2(times)),colour="green4",linetype=2,size=0.3)
	p<-p+theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5))
	p<-p+labs(x="log2(Fold Change)",y=paste("-log10(",names(x)[value],")",sep=""),title=title)
	p<-p+theme(legend.position="top",legend.justification=0.9,legend.text = element_text(size = 10))
	p<-p+guides(col = guide_legend(ncol = 1,title=NULL))
	p<-p+theme(plot.margin=margin(15, 40, 10, 10))
	p<-p+xlim(-(max(abs(vol2[,1]))),max(abs(vol2[,1])))
	if(!is.null(labe)){
		lab<-merge(vol2,labe,by.x=4,by.y=1)
		lab$du2<-gsub("up",color[1],lab$du2)
		lab$du2<-gsub("down",color[3],lab$du2)
		lab$du2<-gsub("nde",color[2],lab$du2)
		p<-p+geom_text(data=lab,aes(x=lab[,2],y=lab[,3],label=lab[,5]),colour="black",inherit.aes = FALSE,show.legend=F,size=2,vjust =-0.5)
	}
	return(p)
}

dif <- read.table("data.txt", header=T, sep="\t", check.names=FALSE, colClasses = "character")
dif$log2FC <- as.numeric(dif$log2FC)
dif$p_value <- as.numeric(dif$p_value)
dif[,which(names(dif)==paste(Test_name,"_value",sep=""))] <- as.numeric(dif[,which(names(dif)==paste(Test_name,"_value",sep=""))])
dif[,which(names(dif)==paste(Control_name,"_value",sep=""))] <- as.numeric(dif[,which(names(dif)==paste(Control_name,"_value",sep=""))])

p <- scatter(dif,test=which(names(dif)==paste(Test_name,"_value",sep="")),control=which(names(dif)==paste(Control_name,"_value",sep="")),times=FoldChange,xlab=Control_name,ylab=Test_name,labe=NULL,log2FC=which(names(dif)=="log2FC"))
ggsave(paste("scatter_", Test_name, "_vs_", Control_name,".pdf", sep=""), plot=p, device="pdf", width=7, height=7, units="in")
ggsave(paste("scatter_", Test_name, "_vs_", Control_name,".png", sep=""), plot=p, device="png", width=7, height=7, units="in", dpi=600)

p <- volcano(dif,log2FC=which(names(dif)=="log2FC"),value=which(names(dif)=="p_value"),level=pvalue,times=FoldChange,title=paste(Test_name,"vs",Control_name),labe=NULL)
ggsave(paste("volcano_", Test_name, "_vs_", Control_name, ".pdf",sep=""), plot=p, device="pdf", width=7, height=8, units="in")
ggsave(paste("volcano_", Test_name, "_vs_", Control_name, ".png",sep=""), plot=p, device="png", width=7, height=8, units="in", dpi=600)
