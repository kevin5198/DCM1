require(getopt, quiet=TRUE)

opt = getopt(matrix(c(
'file','f',1,'character',
'datatype','d',1,'character',
'num','n',1,'numeric',
'cex','c',1,'numeric',
'main','m',1,'numeric',
'barwidth','b',1,'numeric',
'figwidth','w',1,'numeric',
'figheight','i',1,'numeric',
'out','o',1,'character',
'help','h',0,'logical'
),byrow=TRUE, ncol=4));

usage<-function(){
	cat("
Usage: Rscript GO_seq_BP_CC_MF_barplot.r [-f ./data.txt] [-d Enrichment Score] [-n 10] [-c 0.7] [-m 2] [-b 0.8] [-w 10] [-i 9] [-o ./]
Options:
	-f, --file	data file name, default: ./data.txt
	-d, --datatype	data type, default: Enrichment Score
	-n, --num	the number of term to plot, default: 10
	-c, --cex	font size, default: 0.7
	-m, --main	title size, default: 2
	-b, --barwidth	bar width, default: 0.8
	-w, --figwidth	figure width, default: 10
	-i, --figheight	figure height, default: 9
	-o, --out	output path, default: ./
	-h, --help	display this help and exit.
	\n",sep=" ")
	q(status=1);
}

if (!is.null(opt$help)){stop(usage())}
if (is.null(opt$file)){opt$file <- "./data.txt"}
if (is.null(opt$datatype)){opt$datatype <- "Enrichment Score"}
if (is.null(opt$num)){opt$num <- 10}
if (is.null(opt$cex)){opt$cex <- 0.7}
if (is.null(opt$main)){opt$main <- 2}
if (is.null(opt$barwidth)){opt$barwidth <- 0.8}
if (is.null(opt$figwidth)){opt$figwidth <- 10}
if (is.null(opt$figheight)){opt$figheight <- 9}
if (is.null(opt$out)){opt$out <- "./"}

mydata <- read.table(file = opt$file, header = TRUE, check.names = FALSE, sep = "\t", na.strings="N/A",quote = "")
mydata <- data.frame(mydata)

GOSumPlots <- function(mydata){
	barCol <- c(rep("darkred", opt$num), rep("darkgreen", opt$num), rep("darkblue", opt$num))
	barMES <- mydata[, c("Term", sub(' ','_',opt$datatype))]

	barDES <- barMES[, 2]
	names(barDES) <- c(1:length(barDES))

	legText1 <- NULL
	legText2 <- NULL
	legText3 <- NULL

	if ("darkred" %in% unique(barCol)) legText1 <- "BP"
	if ("darkgreen" %in% unique(barCol)) legText2 <- "CC"
	if ("darkblue" %in% unique(barCol)) legText3 <- "MF"

	legText = c(legText1, legText2, legText3)
	rm(legText1, legText2, legText3)

	if (length(legText) == 3){
		pFix = 3
	} else if (length(legText) == 2){
		pFix = 2
	} else if (length(legText) == 1){
		pFix = 2
	}

	pdf(paste(opt$out,"GO_",sub(' ','',opt$datatype),".pdf",sep=""), width = opt$figwidth, height = opt$figheight)
	par(mai = c(3, 2, 1, 0.7))
	barplot(barDES, width = opt$barwidth, space = 0.25, col = barCol,
			legend.text = legText,
			args.legend = list(x = -pFix, fill = unique(barCol)),
			ylab = opt$datatype, names.arg = "")
	title("Sig GO terms of DE gene", cex.main = 2, font.main= opt$main)
	par(srt=45)
	text(x = seq(0.8, length(barDES) - 0.2, 1), y = -0.2,
		labels = mydata$Term, adj = c(1, 0), xpd = TRUE, cex = opt$cex)
	dev.off()

	png(paste(opt$out,"GO_",sub(' ','',opt$datatype),".png",sep=""), width = opt$figwidth, height = opt$figheight, units = "in", res = 400)
	par(mai = c(3, 2, 1, 0.7))
	barplot(barDES, width = opt$barwidth, space = 0.25, col = barCol,
			legend.text = legText,
			args.legend = list(x = -pFix, fill = unique(barCol)),
			ylab = opt$datatype, names.arg = "")
	title("Sig GO terms of DE gene", cex.main = 2, font.main= opt$main)
	par(srt=45)
	text(x = seq(0.8, length(barDES) - 0.2, 1), y = -0.2,
		labels = mydata$Term, adj = c(1, 0), xpd = TRUE, cex = opt$cex)
	dev.off()

	rm(barMES, barDES, mydata)
}

GOSumPlots(mydata)
