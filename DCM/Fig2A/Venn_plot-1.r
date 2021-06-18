require(VennDiagram, quiet=TRUE)
require(getopt, quiet=TRUE)

opt = getopt(matrix(c(
'file','f',1,'character',
'group','g',1,'character',
'title','t',1,'character',
'out','o',1,'character',
'help','h',0,'logical'
),byrow=TRUE, ncol=4));

usage<-function(){
	cat("
Usage: Rscript Venn_plot.r  [-g group1,group2] [-t title] [-o ./]
Options:
	-g, --group	group name in compare.txt, comma split.
	-t, --title	plot title.
	-o, --out	output path.
	-h, --help	display this help and exit.
	\n",sep=" ")
	q(status=1);
}

if (is.null(opt$group) || is.null(opt$title) || is.null(opt$out) || !is.null(opt$help)){
	stop(usage())
}

groups <- unlist(strsplit(as.character(opt$group),split=","))

name1 <- seq(1,5297)
name2 <- seq(1303,6270)

venn_listdata = list(group1 = name1, group2 = name2)
names(venn_listdata) = c(groups[1], groups[2])

#venn.diagram(venn_listdata,fill=c("orchid3","dodgerblue"),main=opt$title,width=20,height=20,units="cm",main.cex=2,margin=0.05,cex=1.2,cat.cex=1.2,resolution=300,filename=paste(opt$out,"Venn_plot.png",sep=""))
venn.diagram(venn_listdata,fill=c("#7B68EE","#F08080"),main=opt$title,width=20,height=20,units="cm",main.cex=2,margin=0.05,cex=1.2,cat.cex=1.2,resolution=300,filename=paste(opt$out,"Venn_plot.png",sep=""))

venn.plot <- venn.diagram(venn_listdata,fill=c("#7B68EE","#F08080"),main=opt$title,width=20,height=20,units="cm",main.cex=2,margin=0.05,filename = NULL)
pdf(file = paste(opt$out,"Venn_plot.pdf",sep=""))
grid.draw(venn.plot)
dev.off()
