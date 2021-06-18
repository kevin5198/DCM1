require(getopt, quiet=TRUE)

opt = getopt(matrix(c(
'file','f',1,'character',
'out','o',1,'character',
'help','h',0,'logical'
),byrow=TRUE, ncol=4));

usage<-function(){
	cat("
Usage: Rscript barplot.r [-f ./data.txt] [-o ./]
Options:
	-f, --file	data file name.
	-o, --out	output path.
	-h, --help	display this help and exit.
	\n",sep=" ")
	q(status=1);
}

if (is.null(opt$file) || is.null(opt$out) || !is.null(opt$help)){
	stop(usage())
}

MyPlot <- function(){
	df <- read.table(opt$file, header=T, sep="\t")
	data <- rbind(df$Non_m6A, df$Common, df$Unique)
	par(mar = c(3, 5, 2, 2))
	bp <- barplot(
		data,
		col=c("#4682B4", "#8B0000", "#DAA520"),
		ylab="Gene number",
		ylim=c(0,20000),
		yaxt="n",
		cex.main=2,
		cex.names=1,
		cex.axis=1,
		cex.lab=1.5,
		names.arg=c('NC', 'DCM'),
	)
	axis(2,seq(0, 20000, 5000))
	text(0.7, 7204, '14408', col='white', adj=0.5, cex=1.4, srt=0)
	text(0.7, 16023, '3230', col='white', adj=0.5, cex=1.4, srt=0)
	text(0.7, 17954, '633', col='white', adj=0.5, cex=1.4, srt=0)
	text(1.9, 7331, '14663', col='white', adj=0.5, cex=1.4, srt=0)
	text(1.9, 16278, '3230', col='white', adj=0.5, cex=1.4, srt=0)
	text(1.9, 18130, '474', col='white', adj=0.5, cex=1.4, srt=0)
	usr <- par("usr")
	x <- usr[2]*0.12
	y <- usr[4]*1.035
	legend(x, y, legend = c('Unique m6A genes', 'Common m6A genes', 'Non-m6A genes'), fill=c("#DAA520", "#8B0000", "#4682B4"), cex=1, box.col="white", xpd=TRUE)
	dev.off()
}

png(file = paste(opt$out,"barplot.png",sep=""), width=1200, heigh=2400, res=301)
MyPlot()

pdf(file = paste(opt$out,"barplot.pdf",sep=""), width=6, height=12)
MyPlot()
