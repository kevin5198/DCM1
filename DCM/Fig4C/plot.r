plot_four_quadrant <- function() {
	par(mar = c(5, 5, 5, 5))
	data <- read.table('./data.txt', header=T, comment.char="", sep='\t')
	colors <- c("#808080", "#B22222", "#00BFFF", "#FF7F50", "#2E8B57")
	plot(x=data$RNAlog2fc,y=data$meriplog2fc,col=alpha(colors[data$color], data$transparency), type='p', pch=16, cex=data$size, cex.lab=1.2, xlab = 'Log2FC(Gene Expression)', ylab='Log2FC(Methylation)')
	abline(v=-0.263034405833794, lty=2, col = '#808080')
	abline(v=0.263034405833794, lty=2, col = '#808080')
	abline(h=-0.263034405833794, lty=2, col = '#808080')
	abline(h=0.263034405833794, lty=2, col = '#808080')
	text(x=4.5*0.8, y=10*0.8, labels = 'Hyper-up: 50', col='#B22222', cex=1.5)
	text(x=-2.5*0.8, y=10*0.8, labels = 'Hyper-down: 62', col='#00BFFF', cex=1.5)
	text(x=4.5*0.8, y=-4*0.8, labels = 'Hypo-up: 224', col='#FF7F50', cex=1.5)
	text(x=-2.5*0.8, y=-4*0.8, labels = 'Hypo-down: 137', col='#2E8B57', cex=1.5)
	dev.off()
}

png(file = 'FC_1.2_Pvalue_0.05.png', width=2400, heigh=2400, res=301)
plot_four_quadrant()

pdf(file = 'FC_1.2_Pvalue_0.05.pdf', width=8, height=8)
plot_four_quadrant()
