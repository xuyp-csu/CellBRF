#' evalcluster
#'
#' Three criteria (NMI, RI, ARI) to evaluate clustering performance.
#' @param truelabel A numeric vector of true labels of each sample.
#' @param predlabel A numeric vector of predicted labels of each sample.
#' @return NMI Value of normalized mutual information.
#' @return RI Value of rand index.
#' @return ARI Value of adjusted rand index.
#' @keywords clustering, validation
#' @export
#' @author Validation methods implemented by Yunpei Xu, xu_yunpei@csu.edu.cn, Central South University.
#' @examples
#' truelabel = sample(1:3, size=10, replace=TRUE)
#' predlabel = sample(1:3, size=10, replace=TRUE)
#' evalcluster(truelabel,predlabel)



evalcluster<-function(truelabel,predlabel){
  if(length(truelabel)!=length(predlabel))
    stop("truelabel and predlabel must have the same length")
  total = length(truelabel)
  x_ids = unique(truelabel)
  y_ids = unique(predlabel)
  #Mutual information
  MI = 0.0
  for(idx in x_ids){
    for(idy in y_ids){
      idxOccur = which(truelabel==idx)
      idyOccur = which(predlabel==idy)
      idxyOccur = intersect(idxOccur,idyOccur)
      if(length(idxyOccur)>0){
        MI = MI + (length(idxyOccur)/total)*log2((length(idxyOccur)*total)/(length(idxOccur)*length(idyOccur)));
      }
    }
  }
  
  #Normalized Mutual information
  Hx = 0; #Entropies
  for(idx in x_ids){
    idxOccurCount = length(which(truelabel==idx));
    Hx = Hx - (idxOccurCount/total) * log2(idxOccurCount/total);
  }
  Hy = 0;#Entropies
  for(idy in y_ids){
    idyOccurCount = length(which(predlabel==idy));
    Hy = Hy - (idyOccurCount/total) * log2(idyOccurCount/total);
  }
  nmi = 2 * MI / (Hx+Hy)
  
  #(adjusted) Rand Index
  tab = table(truelabel,predlabel)
  conv_df = as.data.frame.matrix(tab)
  n <- sum(tab)
  ni <- apply(tab, 1, sum)
  nj <- apply(tab, 2, sum)
  n2 <- choose(n, 2)
  nis2 <- sum(choose(ni[ni > 1], 2))
  njs2 <- sum(choose(nj[nj > 1], 2))
  ri = 1 + (sum(tab^2) - (sum(ni^2) + sum(nj^2))/2)/n2
  ari=c(sum(choose(tab[tab > 1], 2)) - (nis2 * njs2)/n2)/((nis2 + njs2)/2 - (nis2 * njs2)/n2)
  
  out = c(nmi,ri,ari)
  names(out)=c("NMI","RI","ARI")
  return(out)
}

