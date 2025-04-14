#' Calculate niche-niche connections with connectome data from nichetag method
#'
#' This fuction use tag expression matrix to infer the connection among niches labeled with different nichetags
#'
#' @param tag_expression    dataframe, input tag expression matrix in which cells are rownames and tag expression are colnames
#' @param cell_clusters    factor, input clusters of all cells in tag_expression matrix
#' @param share_method    min,max,mean, mothed to calculate connection strength between different clones
#' @param direction    TorF, if or not to distinguish tag-sender or tag-receiver
#'
#' @return    a list contains connectome information of all niches
#' @importFrom stats complete.cases sd ecdf
#' @export
#'
#' @examples
#' nnt=Dnichenetwork(tag_expression,cell_clusters)
#' summary(nnt)
#' nnt=nichenetwork(tag_expression,cell_clusters,share_method="mean")
#' dnnt=Dnichenetwork(tag_expression,cell_clusters,direction = T)
Dnichenetwork<-function(tag_expression,cell_clusters,share_method="min",direction=FALSE){
  res=list()
  tag_expression=tag_expression[,order(apply(tag_expression,2,sum),decreasing = T)]
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
