#' print niche network
#'
#' @param nnt    a list contains connectome information of all niches, the result of Dnichenetwork
#' @param file    file names to be saved
#' @param vsize    vertex size
#' @param esize    edge size
#' @param vertex.label.cex    vertex label size
#' @param axes   TorF, if or not draw axes
#' @param width    file size
#' @param height    file size
#' @param weighted    layout caculation parameter
#' @param direction    layout caculation parameter
#'
#' @return    pdf file
#' @import igraph
#' @importFrom scales hue_pal
#' @export
#'
#' @examples
#' print_nichenet(nnt,file="nichenet-w2.pdf",vsize = 10)
print_nichenet<-function(nnt,file="nichenet.pdf",vsize=1,esize=1,
                         vertex.label.cex =1, axes=FALSE,
                         width=7,height=7,weighted=FALSE,direction=FALSE){

  niche.size=nnt[["niche_size"]]
  niche_share=nnt[["niche_interaction"]]
  g <- graph_from_data_frame(niche_share[,1:2], directed = direction)
  # 设置顶点大小
  V(g)$size <- log10(1+niche.size[V(g)$name])  # 通过节点名称匹配大小
  E(g)$size <- log10(1+niche_share[,3])

  nicheID=1:length(V(g)$name)
  names(nicheID)=V(g)$name

  #library(matlab)
  #color_pallete=jet.colors(length(V(g)$name))
  #library(scales)
  color_pallete=hue_pal()(length(V(g)$name))

  pdf(file,width = width,height = height)
  #E(g)$curved <- seq(-0.5, 0.5, length.out=ecount(g))
  #E(g)$curved <- rep( 1, ecount(g))
  # 绘制网络图
  set.seed(666)
  if(weighted){
    lo=layout_with_fr(g,  grid = "nogrid", niter = 100000, weights = E(g)$size)
    #print(lo)
    #print(norm_coords(lo))
  }else{
    lo=layout_with_fr(g, grid = "nogrid", niter = 100000)
  }
  plot(g,
       axes = axes,
       layout = lo,
       vertex.frame.width = 0.5,
       vertex.label = V(g)$name,
       vertex.label.cex = vertex.label.cex,
       edge.arrow.size = 0.2,
       edge.width = esize*E(g)$size,  # 使用 size 作为边的宽度
       vertex.size = vsize*V(g)$size, # 使用 vector 作为顶点大小
       edge.color = "grey",
       vertex.color = color_pallete)
  dev.off()
}

#' print clone network
#'
#' @param nnt a list contains connectome information of all niches, the result of Dnichenetwork
#' @param file pdf file name to be saved
#' @param mark_groups igrap parameter
#' @param margin igraph parameter
#' @param vertex.label.cex vertex label size
#' @param vsize vertex size
#' @param esize edge size
#' @param ecolor edge color
#' @param width file size
#' @param height file size
#' @param weighted layout parameter
#' @param direction layout parameter
#' @param niter layout parameter
#' @param axes TorF, if or not draw axes
#'
#' @return df.niche2cloneID
#' @import igraph
#' @importFrom Polychrome createPalette
#' @importFrom scales hue_pal
#' @importFrom stringr str_split_fixed str_split
#' @export
#'
#' @examples
#' print_nichenetwork(nnt,file="nichenetwork-w3-26.pdf",vertex.label.cex=0.1,vsize=1.3,esize=1,ecolor=adjustcolor("grey80", alpha.f = 0.5),mark_groups = T)
#' print_nichenetwork(nnt,file="nichenetwork-w3-27.pdf",vertex.label.cex=0.1,vsize=1.3,esize=1,ecolor=NA,mark_groups = T)
#' print_nichenetwork(nnt,file="nichenetwork-w3-28.pdf",vertex.label.cex=0.1,vsize=1.3,esize=1,ecolor=NA)
#' print_nichenetwork(nnt,file="nichenetwork-w3-29.pdf",vertex.label.cex=0.1,vsize=1.3,esize=1)

print_nichenetwork<-function(nnt,file="nichenetwork.pdf",
                             mark_groups = FALSE,margin = c(0, 0, 0, 0),
                             vertex.label.cex=1,vsize=1,esize=1,ecolor="grey80",
                             width=7,height=7,
                             weighted=FALSE,direction=FALSE,niter = 100000,
                             axes=FALSE){
  shareData2=nnt[["clone_interaction"]]
  clone.size=nnt[["clone_size"]]
  celltypes=nnt[["code2celltype"]]
  niche=nnt[["niche"]]
  #将互作clone分开
  shareCloneAB=as.data.frame(str_split_fixed(shareData2$clonePair,"-",n=2))
  #shareData3=data.frame(c1=shareCloneAB[,1],c2=shareCloneAB[,2],edge)
  #取出所有包含交叉的clone,无互作的克隆不在igraph图中显示
  shareCloneID=unique(c(shareCloneAB[,1],shareCloneAB[,2]))
  #包含交叉的clone的size用其中的细胞数表示
  shareCloneSize=clone.size[shareCloneID]
  #将shareCloneID，shareCloneSize等作为输出结果输出到list
  #使用igraph画clone互作network，input是shareCloneAB，shareData2$share，shareCloneSize
  # 创建图对象
  g <- graph_from_data_frame(shareCloneAB, directed = direction)
  cloneID=1:length(V(g)$name)
  names(cloneID)=V(g)$name
  niche2cloneID=lapply(niche,function(x){cloneID[intersect(names(x),names(cloneID))]})
  newniche=lapply(niche,function(x){x[intersect(names(x),names(cloneID))]})
  # 设置顶点大小
  V(g)$size <- vsize*log10(1+shareCloneSize[V(g)$name])  # 通过节点名称匹配大小
  #E(g)$size <- log2(1+shareData2$share)
  E(g)$size <- shareData2$share

  pdf(file,width = width,height = height)
  last_code<- str_split_fixed(V(g)$name, "_",n=2)[,2]
  color_pallete <- createPalette(nchar(last_code[1]),
                                 seedcolors = c("#E69F00" ,"#56B4E9" ,"#009E73", "#F0E442", "#0072B2", "#D55E00","#CC79A7"))
  char_vectors <- strsplit(last_code, split = "")
  #print(char_vectors)
  clone_df=sapply(char_vectors,as.integer)
  colors=apply(clone_df,2,function(x){color_pallete[as.logical(x)]})
  clone_cluster=apply(clone_df,2,function(x){nnt[["code2celltype"]][as.logical(x)]})
  # 绘制网络图
  set.seed(666)
  if(weighted){
    lo=layout_with_fr(g,  grid = "nogrid", niter = niter, weights = E(g)$size)
    #print(lo)
    #print(norm_coords(lo))
  }else{
    lo=layout_with_fr(g, grid = "nogrid", niter = niter)
  }

  if(mark_groups){
    #library(Polychrome)
    first_code <-str_split_fixed(V(g)$name, "_",n=2)[,1]
    char_vectors <- strsplit(first_code, split = "")
    clone_df=as.data.frame(sapply(char_vectors,as.integer))
    colnames(clone_df)=1:ncol(clone_df)
    clone_1tag=clone_df[,which(apply(clone_df,2,sum)==1)]
    clone_tag=apply(clone_1tag,2,function(x){nnt[["code2tag"]][as.logical(x)]})
    clone1tag_list=list()
    for(tag in nnt[["code2tag"]]){
      clone1tag_list[[tag]]=as.integer(names(clone_tag[clone_tag==tag]))
    }
    mark.groups=clone1tag_list[sapply(clone1tag_list,length)>0]
    #library(scales)
    mark.col=hue_pal()(length(clone1tag_list))
  }else{
    mark.groups = list()
    mark.col = rainbow(length(mark.groups), alpha = 0.3)
  }

  plot(g,
       axes = axes,
       layout = lo,
       mark.groups = mark.groups,
       mark.col=mark.col,
       #mark.border=mark.colors,
       mark.shape =1/2,
       vertex.frame.width = 0.5,
       vertex.label = 1:length(V(g)$name),
       vertex.label.cex = vertex.label.cex,
       edge.arrow.size = 0.2,
       edge.width = esize*log10(1+E(g)$size),  # 使用 size 作为边的宽度
       vertex.size = V(g)$size, # 使用 vector 作为顶点大小
       edge.color = ecolor,
       vertex.color = colors,
       margin = margin
  )

  legend("topleft", legend = celltypes,col = color_pallete,pch = 21, pt.bg = color_pallete, pt.cex = 1,cex = 0.5,bty = "n")

  if(mark_groups){legend("topright", legend = names(clone1tag_list),col = mark.col,
                         pch = 21, pt.bg = mark.col, pt.cex = 1,cex = 0.5,bty = "n")}

  dev.off()

  df.niche2cloneID=data.frame()
  for(name in names(niche2cloneID)){
    if(length(newniche[[name]])>0){
      df=data.frame(tag=name,tagnum=newniche[[name]],
                    cloneID=niche2cloneID[[name]],
                    code=names(niche2cloneID[[name]]))
      df.niche2cloneID=rbind(df.niche2cloneID,df)
    }
  }
  rownames(df.niche2cloneID)=NULL
  return(df.niche2cloneID)
}



#' print tag vs cell number and types
#'
#' @param nnt a list contains connectome information of all niches, the result of Dnichenetwork
#' @param file pdf
#'
#' @return pdf
#' @import ggplot2
#' @import cowplot
#' @export
#'
#' @examples
#' print_nichetag(nnt)
print_nichetag<-function(nnt,file="nichetag.pdf"){
  #library(cowplot)
  niche2tagnum=nnt[["niche_tags"]]
  niche2cellnum=nnt[["niche_size"]]
  niche2typenum=nnt[["niche_celltypes"]]
  data0=data.frame(celltag=names(niche2tagnum),tagnumber=log10(niche2tagnum))
  data0$celltag=factor(data0$celltag,levels = data0$celltag)
  data1=data.frame(celltag=names(niche2cellnum),cellnumber=niche2cellnum)
  data1$celltag=factor(data1$celltag,levels = data1$celltag)
  data2=data.frame(celltag=names(niche2typenum),celltypes=niche2typenum)
  data2$celltag=factor(data2$celltag,levels = data2$celltag)
  p0 <- ggplot(data0, aes(x = celltag, y =tagnumber)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    theme_minimal() + ylab("log10(tagnumber)")+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  p1 <- ggplot(data1, aes(x = celltag, y =cellnumber)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  p2 <- ggplot(data2, aes(x = celltag, y =celltypes)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  pg=plot_grid(p0,p1,p2,ncol=1,nrow = 3,hjust = "hv")
  ggsave(plot = pg,filename=file)
}


#' print cluster vs tag number and types
#'
#' @param nnt a list contains connectome information of all niches, the result of Dnichenetwork
#' @param file pdf
#'
#' @return pdf
#' @import ggplot2
#' @import cowplot
#' @export
#'
#' @examples
#' print_clustertag(nnt)
print_clustertag<-function(nnt,file="clustertag.pdf"){
  tag_expression=nnt[["tag_matrix"]]
  cell_clusters=droplevels(nnt[["cell_type"]])
  plot.list=list()
  for (cluster_id in unique(cell_clusters)){
    print(cluster_id)
    cluster_cells <- names(cell_clusters)[cell_clusters == cluster_id]
    tag_expression_cluster <- tag_expression[cluster_cells, ]
    tag_counts_per_cell <- rowSums(tag_expression_cluster)
    unique_tag_counts  <- rowSums(tag_expression_cluster > 0)
    tag_counts_df <- data.frame(
      tag_count = tag_counts_per_cell,
      unique_tag_count = unique_tag_counts)
    tag_count_plot <- ggplot(tag_counts_df, aes(x = tag_count)) +
      geom_density(fill = "blue", alpha = 0.5) +
      labs(title = cluster_id,
           x = "Number of Tags per Cell",
           y = "Density") +
      theme_minimal()
    plot.list[[paste0(cluster_id,"_1")]]=tag_count_plot
    unique_tag_count_plot <- ggplot(tag_counts_df, aes(x = unique_tag_count)) +
      geom_density(fill = "red", alpha = 0.5) +
      labs(title = cluster_id,
           x = "Kinds of Tags per Cell",
           y = "Density") +
      theme_minimal()
    plot.list[[paste0(cluster_id,"_2")]]=unique_tag_count_plot
  }
  #library(cowplot)
  pg=plot_grid(plotlist = plot.list,ncol=2,align = "hv")
  ggsave(plot = pg,filename = file,height = 1.5*length(unique(cell_clusters)))
}
