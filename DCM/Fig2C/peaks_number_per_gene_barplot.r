require(getopt, quiet=TRUE)

opt = getopt(matrix(c(
'file','f',1,'character',
'name','n',1,'character',
'type','t',1,'character',
'out','o',1,'character',
'help','h',0,'logical'
),byrow=TRUE, ncol=4));

usage<-function(){
	cat("
Usage: Rscript utr5_cds_utr3_barplot.r [-f ./data.txt] [-n sample] [-t mRNA] [-o ./]
Options:
	-f, --file	data file name.
	-n, --name	sample or compare name.
	-t, --type	Transcript type.
	-o, --out	output path.
	-h, --help	display this help and exit.
	\n",sep=" ")
	q(status=1);
}

if (is.null(opt$file) || is.null(opt$name) || is.null(opt$type) || is.null(opt$out) || !is.null(opt$help)){
	stop(usage())
}

MyPlot <- function(){
	df <- read.table(opt$file, header=T, sep="\t")
	data <- df$count
	par(mar = c(5, 5, 5, 8))
	ymax <- 3000
	ytick <- 1000
	bp <- barplot(
		data,
		col="#000000",
		xlab="Number of m6A peaks",
		ylab="Number of genes",
		ylim=c(0,ymax),
		yaxt="n",
		names.arg = c("1","2","3","4","5","6","7","8","9","10",">10"),
		cex.main=2,
		cex.names=1.2,
		cex.axis=1.2,
		cex.lab=1.5
	)
	axis(2,seq(0, ymax, ytick))
	text(bp, data+max(data)*1.15*0.025, paste(round(data/sum(data)*100,2),"%",sep=""), adj=0.5, cex=1.2, srt=0, xpd=TRUE) # 以原图的x轴为横坐标，而不是数值，数值对应不上
	usr <- par("usr")
	x <- usr[2]*0.98
	y <- usr[4]*1.05
	dev.off()
}

png(file = paste(opt$out,opt$name,"_peaks_per_",opt$type,".png",sep=""), width=3800, heigh=2400, res=301)
MyPlot()

pdf(file = paste(opt$out,opt$name,"_peaks_per_",opt$type,".pdf",sep=""), width=12, height=8)
MyPlot()
