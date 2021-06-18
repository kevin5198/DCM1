#!/home/novelbio/wangzl/software/R/bin/R
# /home/novelbio/wangzl/software/R/bin/Rscript cumulative.R m6a mrna group (all gene) m6a/ac4c
# /home/novelbio/wangzl/software/R/bin/Rscript cumulative.R diff/diff_peak_ann.xls T1VSNC.log2FC1.FDR0.05.fpkm.anno.xls T1NC m6a/ac4c
library(ggplot2)
args=commandArgs(T)
#m6a ctrl treat all
dp=read.table(args[1],header=TRUE,sep="\t",quote="",check.names=FALSE)
dpm=dp[dp$gene_biotype=="protein_coding",]
#mrna
de=read.table(args[2],header=TRUE,sep="\t",quote="",check.names=FALSE,row.names=1)
dem=de[de$gene_biotype=="protein_coding",]
modi<-args[4] #modification
dem$group=paste0(modi,"(-)")
dem[na.omit(match(as.character(unique(dpm$ensembl_gene_id)),rownames(dem))),]$group=paste0(modi,"(+)")
color<-c("black","red")
p<-ggplot(dem, aes(log2FC, colour = group)) + stat_ecdf(geom="line",size=2,)+theme_bw()+ylab("Cumulative fraction")+xlab('log2FC')+
  theme(legend.position = c(0.05,0.95),legend.justification =c(0.05,1),
        legend.title=element_blank(),
        legend.text=element_text(size=15,colour = c("black","red")),
        plot.title =  element_text(size=20,lineheight=2, face="bold",hjust = 0.5,vjust = 1),
        axis.text=element_text(size=15,face = "bold.italic"),
        panel.grid=element_blank(),
        plot.margin=unit(rep(3,4),'lines'),
        axis.title= element_text(size=15, color="black", face= "bold", vjust=0.5, hjust=0.5))+
  scale_color_manual(values = color)+
  ggtitle("Cumulative Differential \nmRNA Abundance")+
  geom_hline(yintercept = c(0,0),linetype=2,size=1)+
  geom_hline(yintercept = c(1,0),linetype=2,size=1)+
  geom_vline(xintercept = 0,linetype=2,size=1)
ks<-ks.test(dem[dem$group==paste0(modi,"(+)"),]$log2FC,dem[dem$group==paste0(modi,"(-)"),]$log2FC)
print(ks)

ggsave(paste0(args[3],"_Cumulative_Distribution_Fraction.pdf"),p)