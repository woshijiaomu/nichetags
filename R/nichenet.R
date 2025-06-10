#' Calculate niche-niche connections with connectome data from nichetag method
#'
#' This fuction use tag expression matrix to infer the connection among niches labeled with different nichetags
#'
#' @param scObject    Seurat, input seurat object in which tag expression matrix is merged in gene expression matrix
#' @param groupby    character, a colname in meta data of seurat object, used for clone defination
#' @param share_method    min,max,mean, mothed to calculate connection strength between different clones
#' @param direction    TorF, if or not to distinguish tag-sender or tag-receiver
#'
#' @return    a list contains connectome information of all niches:
#' \describe{
#'   \item{tag_matrix}{data.frame, tag expression matrix of all cells}
#'   \item{cell_type}{factor, cell type or state used for clone defination}
#'   \item{cell_tags}{numeric, total expression of each tag}
#'   \item{niche_size}{integer, covered cell number of each tag}
#'   \item{niche_celltypes}{integer, covered cell types or cell states of each tag}
#'   \item{niche_interaction}{data.frame, overlapped cell numbers of tag-tag pairs}
#'   \item{code2celltype}{character,corresponding celltype of each celltypes code positons}
#'   \item{code2tag}{character, corresponding tag of each tag code positons}
#'   \item{cell_barcode}{character, code of each cell, tagcode+celltypecode}
#'   \item{clone}{list, eche element represent a tag expression matrix of each clone}
#'   \item{clone_size}{integer, covered cell number of each clone}
#'   \item{niche}{list, mean tag expression of covered cells of all covered clone of each tag}
#'   \item{clone_interaction}{data.frame, overlapped tag expression of clone-clone pairs}
#'   \item{cloneID}{integer, order of clones when make graph with graph_from_data_frame in igraph }
#' }
#' @importFrom SeuratObject FetchData
#' @importFrom stats complete.cases sd ecdf
#' @importFrom igraph graph_from_data_frame
#' @export
#'
#' @examples
#' nnt=Dnichenetwork(scObject,groupby="cell_clusters")
#' summary(nnt)
#' nnt=nichenetwork(scObject,groupby="cell_clusters",share_method="mean")
#' dnnt=Dnichenetwork(scObject,groupby="cell_clusters",direction = T)
Dnichenetwork<-function(scObject,groupby="cell_clusters",share_method="min",direction=FALSE){
  res=list()
  tags <- grep("^[ATCG]{8}$", rownames(scObject), value = TRUE)
  tag_expression <- FetchData(scObject, vars = tags,layer= "counts")
  tag_expression=quality_control(tag_expression)
  tag_expression=tag_expression[,order(apply(tag_expression,2,sum),decreasing = T)]
  cell_clusters=scObject@meta.data[,groupby]
  names(cell_clusters)=colnames(scObject)
  cell_clusters=droplevels(cell_clusters)
  res[["tag_matrix"]]=tag_expression
  res[["cell_type"]]=cell_clusters
  niche=apply(tag_expression,2,function(x){x[x>0]})
  niche.size=sapply(niche,function(x){length(x)})
  res[["niche_tags"]]=apply(tag_expression,2,sum)
  res[["niche_size"]]=niche.size
  niche.celltypes=sapply(niche,function(x){length(unique(cell_clusters[names(x)]))})
  res[["niche_celltypes"]]=niche.celltypes
  niche_pairs=combn(names(niche),2)
  niche_share=data.frame()
  for(n in 1:ncol(niche_pairs)){
    n1=niche_pairs[1,n]
    n2=niche_pairs[2,n]
    novlp=intersect(names(niche[[n1]]),names(niche[[n2]]))
    dfovlp=cbind(niche[[n1]][novlp],niche[[n2]][novlp])
    niche_share=rbind(niche_share,df2ovlp=data.frame(niche1=n1,niche2=n2,value=sum(apply(dfovlp,1,min))))
  }
  rownames(niche_share)=NULL
  niche_share=niche_share[niche_share$value>0,]
  res[["niche_interaction"]]=niche_share
  cell.typeBarcode<-function(cell_clusters){
    cell.cluster=droplevels(cell_clusters)
    cell.cluster.code=as.data.frame(matrix(rep(0,length(cell.cluster)*length(unique(cell.cluster))),
                                           ncol=length(unique(cell.cluster))))
    colnames(cell.cluster.code)=sort(unique(cell.cluster))
    rownames(cell.cluster.code)=names(cell_clusters)
    for(cc in colnames(cell.cluster.code)){
      index=which(cell.cluster==cc)
      cell.cluster.code[index,cc]=1
    }
    cell.typeBarcode=apply(cell.cluster.code,1,function(x){paste0(x,collapse="")})
    res=list()
    res[["typeBarcode"]]=cell.typeBarcode
    res[["celltype"]]=colnames(cell.cluster.code)
    res
  }
  cell.typeBarcode.list=cell.typeBarcode(cell_clusters)
  cell.typeBarcode=cell.typeBarcode.list[["typeBarcode"]]
  res[["code2celltype"]]=cell.typeBarcode.list[["celltype"]]
  cell.cloneBarcode=apply(tag_expression,1,function(x){paste0(as.integer(x>0),collapse="")})
  res[["code2tag"]]=colnames(tag_expression)
  cell.2barcode=paste(cell.cloneBarcode,cell.typeBarcode,sep="_")
  res[["cell_barcode"]]=cell.2barcode
  clone=list()
  for(cloneBarcode in unique(cell.2barcode)){
    clone[[cloneBarcode]]=tag_expression[cell.2barcode==cloneBarcode,]
  }
  clone.size=sapply(clone,function(x){nrow(x)})
  res[["clone"]]=clone
  res[["clone_size"]]=clone.size
  clone.mean=data.frame(t(sapply(clone,function(x){apply(x,2,mean)})))
  clone.mean.barcode.cell=apply(clone.mean,2,function(x){sort(x[x>0],decreasing = T)})
  res[["niche"]]=clone.mean.barcode.cell
  if(direction==T){
    overlap=lapply(clone.mean.barcode.cell,findoverlap)
    sr2value=data.frame()
    for(ol in names(overlap)){
      if(!all(is.na(overlap[[ol]]))){
        sr2value=rbind(sr2value,data.frame(sr=names(overlap[[ol]]),value=overlap[[ol]],tag=ol))
      }
    }
    rownames(sr2value)=NULL
    shareData=data.frame(clonePair=sr2value$sr,share=sr2value$value)
  }

  if(direction==F){
    print(direction)
    print(share_method)
    if(share_method=="min"){
      shareNumber=unlist(lapply(clone.mean.barcode.cell,function(a){
        if(length(a)>1){apply(combn(a,2),2,min)}else{return(NA)}
      }))
    }else if(share_method=="mean"){
      shareNumber=unlist(lapply(clone.mean.barcode.cell,function(a){
        if(length(a)>1){apply(combn(a,2),2,mean)}else{return(NA)}
      }))
    }else if(share_method=="max"){
      shareNumber=unlist(lapply(clone.mean.barcode.cell,function(a){
        if(length(a)>1){apply(combn(a,2),2,max)}else{return(NA)}
      }))
    }else{
      stop("share_method wrong")
    }
    shareClone=unlist(lapply(clone.mean.barcode.cell,function(a){if(length(a)>1){apply(combn(names(a),2),2,function(x){paste(sort(x),collapse="-")})}else{return(NA)}}))
    shareNumber=shareNumber[!is.na(shareNumber)]
    shareClone=shareClone[!is.na(shareClone)]
    shareData=data.frame(clonePair=shareClone,share=shareNumber)
  }

  shareData2=aggregate(share ~ clonePair, data = shareData, sum)
  res[["clone_interaction"]]=shareData2
  shareCloneAB=as.data.frame(str_split_fixed(shareData2$clonePair,"-",n=2))
  g <- graph_from_data_frame(shareCloneAB, directed = direction)
  cloneID=1:length(V(g)$name)
  names(cloneID)=V(g)$name
  res[["cloneID"]]=cloneID
  return(res)
}


#' find overlap of clones
#'
#' @param clonetag     clones
#' @param share_method    min,max,mean, mothed to calculate connection strength between different clones
#'
#' @return    overlap
#' @export
#'
#' @examples
#' findoverlap()
findoverlap<-function(clonetag,share_method="min"){
  clonetag=sort(clonetag,decreasing = T)
  sender=names(clonetag)[clonetag>=10]
  receivor=setdiff(names(clonetag),sender)
  if(length(sender)>0 & length(receivor)>0){
    pairs=combn(names(clonetag),2)
    setaside=which(pairs[1,] %in% sender)
    pairs2=pairs[,setaside,drop=F]
    print(dim(pairs))
    print(dim(pairs2))
    pairstag=rbind(clonetag[pairs2[1,]],clonetag[pairs2[2,]])
    if(share_method=="min"){shareNumber=apply(pairstag,2,min)}
    if(share_method=="mean"){shareNumber=apply(pairstag,2,mean)}
    if(share_method=="max"){shareNumber=apply(pairstag,2,max)}
    names(shareNumber)=apply(pairs2,2,function(x){paste(x,collapse="-")})
    return(shareNumber)
  }else{
    return(NA)
  }
}

#' Extract clone connections from nnt for Gephi analysis
#'
#' @param nnt a list contains connectome information of all niches
#' @param file output file name
#'
#' @returns a file for Gephi to get more circular layout
#' @importFrom stringr str_split_fixed
#' @export
#'
#' @examples
#' Gephi_prepare(nnt)
Gephi_prepare<-function(nnt,file="Gephi_input.csv"){
  id.df=str_split_fixed(nnt$clone_interaction$clonePair,"-",n=2)
  id.df2=apply(id.df,2,function(x){nnt$cloneID[x]})
  res=cbind(id.df2,nnt$clone_interaction$share)
  colnames(res)=c("source","target","weight")
  write.csv(res,file,row.names = F)
}



#' Get clone expression matrix, for each gene use mean value of all cells of the same clone
#'
#' @param scObject Seurat, input seurat object in which tag expression matrix is merged in gene expression matrix
#' @param nnt a list contains connectome information of all niches
#' @param seurat_layer layer in scObject used to calculate expression level of genes
#'
#' @returns clone expression matrix of all nontag genes
#' @importFrom SeuratObject FetchData
#' @export
#'
#' @examples
#' clonematrix=Clone_expr(scObject,nnt,seurat_layer="counts")
Clone_expr<-function(scObject,nnt,seurat_layer="counts"){
  nontags <- grep("^[ATCG]{8}$", rownames(scObject), value = TRUE,invert = TRUE)
  expr=FetchData(scObject, vars = nontags,layer= seurat_layer)
  res=sapply(nnt[["clone"]],function(x){
    cells=rownames(x)
    apply(expr[cells,],2,mean)
  })
  res=as.data.frame(res)
  res=res[,colnames(res) %in% names(nnt$cloneID)]
  colnames(res)=paste0("Clone",nnt[["cloneID"]][colnames(res)])
  res
}

