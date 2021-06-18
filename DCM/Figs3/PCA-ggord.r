library(ggord)
library(ggplot2)
library(egg)
pca_data=read.table('data.txt',header=T,sep='\t',row=1)
pca_data=t(as.matrix(pca_data))
sample.name <- gsub("^X", "", rownames(pca_data))
sample.name <- sub("\\.", "-", sample.name)
pca_group=factor(c(rep('Wt',5),rep('db_db',5)))
pca=prcomp(pca_data,center=TRUE,retx=T)
p <- ggord(pca, grp_in=pca_group,size=10,txt=NULL,vectyp=0,obslab=FALSE,parse=TRUE)
ggsave(file="PCA_plot.png", plot=set_panel_size(p=p+geom_text(label=sample.name,size=3)+theme(text=element_text(size=15)), width=unit(24,'cm'), height=unit(18,'cm')), width=12, height=9, dpi=300)
ggsave(file="PCA_plot.pdf", plot=set_panel_size(p=p+geom_text(label=sample.name,size=3)+theme(text=element_text(size=15)), width=unit(24,'cm'), height=unit(18,'cm')), width=12, height=9, dpi=300)
