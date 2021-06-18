#! /usr/bin/Rscript

myPie <- function(file, title){
	par(mai = c(0, 0, 0.5, 0))
	data<-read.table(file,header=T,sep="\t")
	data2<-data$Count
	name<-c()
	for(i in 1:nrow(data)){
		name<-c(name,paste(round(data$Count[i]*100/sum(data$Count),2),"%",sep=""))
	}
	names(data2)<-name
	pie(data2, main=title, clockwise=T, init.angle=0, col=colors, cex=1.5, radius=0.8)
}

colors = c('#808080','#CC3333','#336699','#CD853F','#3CB371')

pdf('All_region_stat.pdf', width=15, height=5)
layout(matrix(c(1,2,3,4,4,4), 2, 3, byrow=TRUE), widths=c(1,1,1), heights=c(5,1))
myPie('./Total_region_stat.txt', 'Total m6A distribution')
myPie('./NC_region_stat.txt', 'NC m6A distribution')
myPie('./DCM_region_stat.txt', 'DCM m6A distribution')
plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
legend(x=-0.6, y=2, legend=c("5'UTR","TSS","CDS","Stop_codon","3'UTR"), text.width=0.2, cex=1.5, fill=colors, ncol=5, box.col="white", xpd=TRUE)

png('All_region_stat.png', width=15, height=5, units='in',res=300)
layout(matrix(c(1,2,3,4,4,4), 2, 3, byrow=TRUE), widths=c(1,1,1), heights=c(5,1))
myPie('./Total_region_stat.txt', 'Total m6A distribution')
myPie('./NC_region_stat.txt', 'NC m6A distribution')
myPie('./DCM_region_stat.txt', 'DCM m6A distribution')
plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
legend(x=-0.6, y=2, legend=c("5'UTR","TSS","CDS","Stop_codon","3'UTR"), text.width=0.2, cex=1.5, fill=colors, ncol=5, box.col="white", xpd=TRUE)
