require(ggplot2, quiet=TRUE)
require(getopt, quiet=TRUE)

opt = getopt(matrix(c(
'file','f',1,'character',
'datatype','d',1,'character',
'num','n',1,'numeric',
'size','s',1,'numeric',
'barwidth','b',1,'numeric',
'figwidth','w',1,'numeric',
'percent','p',1,'numeric',
'out','o',1,'character',
'help','h',0,'logical'
),byrow=TRUE, ncol=4));

usage<-function(){
	cat("
Usage: Rscript Pathway_Enrichment_Score_Barplot.r [-f ./data.txt] [-d Enrichment Score] [-n 10] [-s 12] [-b 0.9] [-w 8.5] [-p 0.6] [-o ./]
Options:
	-f, --file	data file name, default: ./data.txt
	-d, --datatype	data type, default: Enrichment Score
	-n, --num	the number of term to plot, default: 10
	-s, --size	font size, default: 12
	-b, --barwidth	bar width, default: 0.9
	-w, --figwidth	figure width, default: 8.5
	-p, --percent	figure height(percent of num), default: 0.6
	-o, --out	output path, default: ./
	-h, --help	display this help and exit.
	\n",sep=" ")
	q(status=1);
}

if (!is.null(opt$help)){stop(usage())}
if (is.null(opt$file)){opt$file <- "./data.txt"}
if (is.null(opt$datatype)){opt$datatype <- "Enrichment Score"}
if (is.null(opt$num)){opt$num <- 10}
if (is.null(opt$size)){opt$size <- 12}
if (is.null(opt$barwidth)){opt$barwidth <- 0.9}
if (is.null(opt$figwidth)){opt$figwidth <- 8.5}
if (is.null(opt$percent)){opt$percent <- 0.6}
if (is.null(opt$out)){opt$out <- "./"}

mydata <- read.table(file = opt$file, header = TRUE, check.names = FALSE, sep = "\t", na.strings="N/A",quote = "")
mydata <- data.frame(mydata)

if ( min(mydata$Pvalue) != max(mydata$Pvalue)){
		dotplot <- ggplot(data = mydata, aes(x = mydata$Enrichment.Score, y = reorder(paste(mydata$Term," [",mydata$ID,"]",sep=""), mydata$Enrichment.Score))) + 
				geom_point(aes(colour = Pvalue, size = Count), alpha = 1, position = position_jitter(w = 0.0, h = 0.0)) + 
				labs(x = paste('\n',opt$datatype,' (-log10(p_value))\n'), y = "\n \n", title = '\nSig pathway of DE gene\n') +
				theme(plot.title = element_text(hjust = 0)) +
				theme(axis.text.x = element_text(color = "black", size = opt$size, angle = 0, hjust = 1, vjust = 1), axis.text.y = element_text(color = "black", size = opt$size)) +
				scale_colour_gradient(limits = c(min(mydata$Pvalue), max(mydata$Pvalue)), low = "red", high = "#0066FF");
}else{
		dotplot <- ggplot(data = mydata, aes(x = mydata$Enrichment.Score, y = reorder(paste(mydata$Term," [",mydata$ID,"]",sep=""), mydata$Enrichment.Score))) + 
				geom_point(aes(colour = Pvalue, size = Count), alpha = 1, position = position_jitter(w = 0.0, h = 0.0)) + 
				labs(x = paste('\n',opt$datatype,' (-log10(p_value))\n'), y = "\n \n", title = '\nSig pathway of DE gene\n') +
				theme(plot.title = element_text(hjust = 0)) +
				theme(axis.text.x = element_text(color = "black", size = opt$size, angle = 0, hjust = 1, vjust = 1), axis.text.y = element_text(color = "black", size = opt$size)) +
				scale_colour_gradient(limits = c(min(mydata$Pvalue), max(mydata$Pvalue)*1.01), low = "red", high = "#0066FF");
	}
	ggsave(paste(opt$out,sub(' ','_',opt$datatype), '.png', sep = ''), dotplot, width = opt$figwidth, height = opt$percent * (opt$num + 1), dpi = 300)
	ggsave(paste(opt$out,sub(' ','_',opt$datatype), '.pdf', sep = ''), dotplot, width = opt$figwidth, height = opt$percent * (opt$num + 1), dpi = 300)
