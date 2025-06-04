# nichetags install guide
you can install it from github:

```
install.packages(c("Polychrome","stringdist","stringr","ggplot2","cowplot","igraph"))
if(!require(devtools)){
	install.packages("devtools")
}
if(!require(nichetags)){
	library(devtools)
	install_github("woshijiaomu/nichetags")
}
```

after installation，attach the package in R：

```
library(nichetags)
```

Quick Start Guide：

```
#attach R package,load sample data
library(nichetags)
data(exam)

#calculate connectome
nnt=Dnichenetwork(tag_expression,cell_clusters)
summary(nnt)
#clone and the cells it contains
clone2cell=nnt$clone
clone2cell.vec=unlist(lapply(names(clone2cell),function(x){
  y=clone2cell[[x]]
  res=rep(x,nrow(y))
  names(res)=rownames(y)
  res
}))
clone2cell.df=data.frame(clone_code=clone2cell.vec,cell=names(clone2cell.vec))
clone2cell.df$clone_ID=nnt$cloneID[clone2cell.df$clone_code]
write.csv(clone2cell.df,file = "clone2cell.csv")

niche2clone=nnt$niche
niche2clone.vec=unlist(lapply(names(niche2clone),function(x){
  res=paste(x,names(niche2clone[[x]]),sep=":")
  res
}))
niche2clone.df=as.data.frame(str_split_fixed(niche2clone.vec,":",n=2))
colnames(niche2clone.df) = c("niche","clone_code")
niche2clone.df$clone_ID=nnt$cloneID[niche2clone.df$clone_code]
write.csv(niche2clone.df,file = "niche2clone.csv")

#nnt=nichenetwork(tag_expression,cell_clusters,share_method="mean")
#dnnt=Dnichenetwork(tag_expression,cell_clusters,direction = T)

#draw quanlity control figures
print_nichetag(nnt, file = "nichetag_qc.pdf")
print_clustertag(nnt, file = "cluster_qc.pdf")
tag_cancer_noncancer(nnt, file = "celltype_qc.pdf")
tag_cci(nnt, file = "cci_qc.pdf")

#draw niche-niche network
print_nichenet(nnt,file="nichenichenetwork.pdf",vsize = 10)

```
