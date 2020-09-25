
ghdist<-function(n,g=0,h=0){
  #
  # generate n observations from a g-and-h dist.
  #
  #   this function is based on "Rallfun-v34.txt"
  #   which can be download from https://dornsife.usc.edu/labs/rwilcox/software/
  x<-rnorm(n)
  if (g>0){
    ghdist<-(exp(g*x)-1)*exp(h*x^2/2)/g
  }
  if(g==0)ghdist<-x*exp(h*x^2/2)
  ghdist
}

.pbos<-function(x,beta=.2,omhatx){
#   Compute the percentage bend correlation between x and y.
#
#   beta is the bending constant for omega sub N.
#
#   this function is based on "Rallfun-v34.txt"
#   which can be download from https://dornsife.usc.edu/labs/rwilcox/software/
  psi<-(x-median(x))/omhatx
  i1<-length(psi[psi<(-1)])
  i2<-length(psi[psi>1])
  sx<-ifelse(psi<(-1),0,x)
  sx<-ifelse(psi>1,0,sx)
  .pbos<-(sum(sx)+omhatx*(i2-i1))/(length(x)-i1-i2)
  .pbos
}

fisher_transfer_test<-function(x1,z1,x0,z0,biasCorrection=TRUE){
#  compute p-value with Fisher Z-transformation test
  n1<-length(x1)
  n0<-length(x0)
  p1=cor(x1,z1,method = 'pearson')
  p0=cor(x0,z0,method = 'pearson')
  if(biasCorrection==TRUE){
    p1=p1-p1/(2*(n1-1))
    p0=p0-p0/(2*(n0-1))
  }

  w1=0.5*log((1+p1)/(1-p1))
  w0=0.5*log((1+p0)/(1-p0))
  delta=sqrt(1/(n1-3)+1/(n0-3))
  T=abs(w1-w0)/delta
  return(2-2*pnorm(T))
}

st1<-function(x1,z1,x0,z0){
  #  Compute p-value for the  equal correlation test
  #  with Pearson correlation based on a logistic regression model
  #  corresponding to two independent groups
  n1<-length(x1)
  n0<-length(x0)
  p1<-n1/(n1+n0)
  p2<-1-p1
  x1_mean<-mean(x1)
  x0_mean<-mean(x0)
  z1_mean<-mean(z1)
  z0_mean<-mean(z0)
  t1=(n1-1)*var(x1)
  t2=(n1-1)*var(z1)
  c1<-sqrt(t1*t2)
  t3=(n0-1)*var(x0)
  t4=(n0-1)*var(z0)
  c0<-sqrt(t3*t4)
  w1<-n1*(x1-x1_mean)*(z1-z1_mean)/c1
  w0<-n0*(x0-x0_mean)*(z0-z0_mean)/c0
  u=sum(w1)*p2-sum(w0)*p1
  w=c(w1,w0)
  u_var=p1*p2*var(w)*(n1+n0-1)
  T=u^2/u_var
  mysign=sign(u)
  signedStat = mysign*T
  pval=1-pchisq(T,1)
  res = list(stat=T, pval=pval, signedStat = signedStat)
  return(res)
}

st2<-function(x1,z1,x0,z0){
  #  Compute p-value for the equal correlation test
  #  with mad-replacing-Pearson correlation
  #  based on a logistic regression model
  #  corresponding to two independent groups
  n1<-length(x1)
  n0<-length(x0)
  p1<-n1/(n1+n0)
  p2<-1-p1
  x1_median<-median(x1)
  x0_median<-median(x0)
  z1_median<-median(z1)
  z0_median<-median(z0)
  c1<-mad(x1)*mad(z1)
  c0<-mad(x0)*mad(z0)
  w1<-(x1-x1_median)*(z1-z1_median)/c1
  w0<-(x0-x0_median)*(z0-z0_median)/c0
  u=sum(w1)*p2-sum(w0)*p1
  w=c(w1,w0)
  u_var=p1*p2*var(w)*(n1+n0-1)
  T=u^2/u_var
  mysign=sign(u)
  signedStat = mysign*T
  pval=1-pchisq(T,1)
  res = list(stat=T, pval=pval, signedStat = signedStat)
  return(res)
}

st3<-function(x1,z1,x0,z0){
  #  Compute p-value for the equal correlation test
  #  with percentage bend correlation
  #  based on a logistic regression model
  #  corresponding to two independent groups
  n1<-length(x1)
  n0<-length(x0)
  p1<-n1/(n1+n0)
  p2<-1-p1
  b=0.2
  x1_median<-median(x1)
  x0_median<-median(x0)
  z1_median<-median(z1)
  z0_median<-median(z0)
  temp1<-sort(abs(x1- x1_median))
  omhatx<-temp1[floor((1-b)*n1)]
  if(omhatx==0){
    omhatx=temp1[n1]
  }
  temp2<-sort(abs(z1- z1_median))
  omhatz<-temp2[floor((1-b)*n1)]
  if(omhatz==0){
    omhatz=temp2[n1]
  }
  a1<-(x1-.pbos(x1,b,omhatx))/omhatx
  b1<-(z1-.pbos(z1,b,omhatz))/omhatz
  a1<-ifelse(a1<=-1,-1,a1)
  a1<-ifelse(a1>=1,1,a1)
  b1<-ifelse(b1<=-1,-1,b1)
  b1<-ifelse(b1>=1,1,b1)

  temp3<-sort(abs(x0- x0_median))
  omhatx<-temp3[floor((1-b)*n0)]
  if(omhatx==0){
    omhatx=temp3[n0]
  }
  temp4<-sort(abs(z0- z0_median))
  omhatz<-temp4[floor((1-b)*n0)]
  if(omhatz==0){
    omhatz=temp4[n0]
  }
  a0<-(x0-.pbos(x0,b,omhatx))/omhatx
  b0<-(z0-.pbos(z0,b,omhatz))/omhatz
  a0<-ifelse(a0<=-1,-1,a0)
  a0<-ifelse(a0>=1,1,a0)
  b0<-ifelse(b0<=-1,-1,b0)
  b0<-ifelse(b0>=1,1,b0)
  c1<-sqrt(sum(a1^2)*sum(b1^2))
  c0<-sqrt(sum(a0^2)*sum(b0^2))
  w1<-n1*a1*b1/c1
  w0<-n0*a0*b0/c0
  u=sum(w1)*p2-sum(w0)*p1
  w=c(w1,w0)
  u_var=p1*p2*var(w)*(n1+n0-1)
  T=u^2/u_var
  mysign=sign(u)
  signedStat = mysign*T
  pval=1-pchisq(T,1)
  res = list(stat=T, pval=pval, signedStat = signedStat)
  return(res)
}

st4<-function(x1,z1,x0,z0){
  #  Compute p-value for the  equal correlation test
  #  with Spearman corrletion
  #  based on a logistic regression model
  #  corresponding to two independent groups
  n1<-length(x1)
  n0<-length(x0)
  p1<-n1/(n1+n0)
  p2<-1-p1
  c1<-n1^2-1
  c0<-n0^2-1
  rank_x1<-rank(x1,ties.method = "average")
  rank_z1<-rank(z1,ties.method = "average")
  rank_x0<-rank(x0,ties.method = "average")
  rank_z0<-rank(z0,ties.method = "average")
  w1=(n1^2-1-6*(rank_x1-rank_z1)^2)/c1
  w0=(n0^2-1-6*(rank_x0-rank_z0)^2)/c0
  u=sum(w1)*p2-sum(w0)*p1
  w=c(w1,w0)
  u_var=p2*p1*var(w)*(n1+n0-1)
  T=u^2/u_var
  mysign=sign(u)
  signedStat = mysign*T
  pval=1-pchisq(T,1)
  res = list(stat=T, pval=pval, signedStat = signedStat)
  return(res)
}



st5=function (x1, z1, x0, z0) {
  # created on May 21, 2020
  #  (1) a revised st5 function with output of test statistic,
  #      p-value, and signed test statistic.
  #      The sign of the test statistic is the same as the sign
  #      of U statistic. Since st5 is the average of st3 and st4,
  #      the sign of U statistic with the larger absolute value is returned
  #
  # test if the corr(x1, z1) is equal to corr(x0, z0) via a robust method
  n1 <- length(x1)
  n0 <- length(x0)
  p1 <- n1/(n1 + n0)
  p2 <- 1 - p1
  b = 0.2
  x1_median <- median(x1)
  x0_median <- median(x0)
  z1_median <- median(z1)
  z0_median <- median(z0)
  temp1 <- sort(abs(x1 - x1_median))
  omhatx <- temp1[floor((1 - b) * n1)]
  if (omhatx == 0) {
    omhatx = temp1[n1]
  }
  temp2 <- sort(abs(z1 - z1_median))
  omhatz <- temp2[floor((1 - b) * n1)]
  if (omhatz == 0) {
    omhatz = temp2[n1]
  }
  a1 <- (x1 - .pbos(x1, b, omhatx))/omhatx
  b1 <- (z1 - .pbos(z1, b, omhatz))/omhatz
  a1 <- ifelse(a1 <= -1, -1, a1)
  a1 <- ifelse(a1 >= 1, 1, a1)
  b1 <- ifelse(b1 <= -1, -1, b1)
  b1 <- ifelse(b1 >= 1, 1, b1)
  temp3 <- sort(abs(x0 - x0_median))
  omhatx <- temp3[floor((1 - b) * n0)]
  if (omhatx == 0) {
    omhatx = temp3[n0]
  }
  temp4 <- sort(abs(z0 - z0_median))
  omhatz <- temp4[floor((1 - b) * n0)]
  if (omhatz == 0) {
    omhatz = temp4[n0]
  }
  a0 <- (x0 - .pbos(x0, b, omhatx))/omhatx
  b0 <- (z0 - .pbos(z0, b, omhatz))/omhatz
  a0 <- ifelse(a0 <= -1, -1, a0)
  a0 <- ifelse(a0 >= 1, 1, a0)
  b0 <- ifelse(b0 <= -1, -1, b0)
  b0 <- ifelse(b0 >= 1, 1, b0)
  c1 <- sqrt(sum(a1^2) * sum(b1^2))
  c0 <- sqrt(sum(a0^2) * sum(b0^2))
  w11 <- n1 * a1 * b1/c1
  w01 <- n0 * a0 * b0/c0
  u1 = sum(w11) * p2 - sum(w01) * p1
  w1 = c(w11, w01)
  u_var1 = p1 * p2 * var(w1) * (n1 + n0 - 1)
  T1 = u1^2/u_var1
  c1 <- n1^2 - 1
  c0 <- n0^2 - 1
  rank_x1 <- rank(x1, ties.method = "average")
  rank_z1 <- rank(z1, ties.method = "average")
  rank_x0 <- rank(x0, ties.method = "average")
  rank_z0 <- rank(z0, ties.method = "average")
  w12 = (n1^2 - 1 - 6 * (rank_x1 - rank_z1)^2)/c1
  w02 = (n0^2 - 1 - 6 * (rank_x0 - rank_z0)^2)/c0
  u2 = sum(w12) * p2 - sum(w02) * p1
  w2 = c(w12, w02)
  u_var2 = p1 * p2 * var(w2) * (n1 + n0 - 1)
  T2 = u2^2/u_var2
  myT = 0.5 * (T1 + T2)
  pval = 1 - stats::pchisq(myT, 1)


  if(T1>T2)
  {
    mysign=sign(u1)
  } else {
    mysign = sign(u2)
  }
  signedStat = mysign*myT

  res = list(stat=myT, pval=pval, signedStat = signedStat)

  return(res)
}

st6<-function(x1,z1,x0,z0){
  #  Compute p-value for the  equal correlation test
  #  with another way to combine Spearman corrletion and percentage bend correlation
  #  based on a multiple logistic regression model
  #  corresponding to two independent groups
  n1<-length(x1)
  n0<-length(x0)
  p1<-n1/(n1+n0)
  p2<-1-p1
  b=0.2
  x1_median<-median(x1)
  x0_median<-median(x0)
  z1_median<-median(z1)
  z0_median<-median(z0)
  temp1<-sort(abs(x1- x1_median))
  omhatx<-temp1[floor((1-b)*n1)]
  if(omhatx==0){
    omhatx=temp1[n1]
  }
  temp2<-sort(abs(z1- z1_median))
  omhatz<-temp2[floor((1-b)*n1)]
  if(omhatz==0){
    omhatz=temp2[n1]
  }
  a1<-(x1-.pbos(x1,b,omhatx))/omhatx
  b1<-(z1-.pbos(z1,b,omhatz))/omhatz
  a1<-ifelse(a1<=-1,-1,a1)
  a1<-ifelse(a1>=1,1,a1)
  b1<-ifelse(b1<=-1,-1,b1)
  b1<-ifelse(b1>=1,1,b1)

  temp3<-sort(abs(x0- x0_median))
  omhatx<-temp3[floor((1-b)*n0)]
  if(omhatx==0){
    omhatx=temp3[n0]
  }
  temp4<-sort(abs(z0- z0_median))
  omhatz<-temp4[floor((1-b)*n0)]
  if(omhatz==0){
    omhatz=temp4[n0]
  }
  a0<-(x0-.pbos(x0,b,omhatx))/omhatx
  b0<-(z0-.pbos(z0,b,omhatz))/omhatz
  a0<-ifelse(a0<=-1,-1,a0)
  a0<-ifelse(a0>=1,1,a0)
  b0<-ifelse(b0<=-1,-1,b0)
  b0<-ifelse(b0>=1,1,b0)
  c1<-sqrt(sum(a1^2)*sum(b1^2))
  c0<-sqrt(sum(a0^2)*sum(b0^2))
  w11<-n1*a1*b1/c1
  w01<-n0*a0*b0/c0
  u1=sum(w11)*p2-sum(w01)*p1
  w1=c(w11,w01)

  c1<-n1^2-1
  c0<-n0^2-1
  rank_x1<-rank(x1,ties.method = "average")
  rank_z1<-rank(z1,ties.method = "average")
  rank_x0<-rank(x0,ties.method = "average")
  rank_z0<-rank(z0,ties.method = "average")
  w12=(n1^2-1-6*(rank_x1-rank_z1)^2)/c1
  w02=(n0^2-1-6*(rank_x0-rank_z0)^2)/c0
  u2=sum(w12)*p2-sum(w02)*p1
  w2=c(w12,w02)

  w1_mean=mean(w1)
  w2_mean=mean(w2)
  a11=sum((w1-w1_mean)^2)
  a12=sum((w1-w1_mean)*(w2-w2_mean))
  a21=a12
  a22=sum((w2-w2_mean)^2)
  u=matrix(c(u1,u2),nrow = 1,ncol = 2,byrow = TRUE)
  cov_u=p1*p2*matrix(c(a11,a12,a21,a22),nrow = 2,ncol = 2,byrow = TRUE)
  st6<-tryCatch( {1-pchisq(u%*%solve(cov_u)%*%t(u),2)},error=function(e){
    1-pchisq(u%*%ginv(cov_u)%*%t(u),1)
  })

  res = list(u=u, pval=st6, cov_u=cov_u)
  return(res)

}


construct_network<-function(es,
                            cor_method = 'st5',
                            var.grp,
                            pseudo_adjust_cutoff = FALSE,
                            pAdjMethod = 'fdr',
                            cutoff = 0.05,
                            nPseudo = 25){
  # input:
  # es: an ExpressionSet object of microRNA dataset
  # cor_method: The method for equal correlation, 'ST5' is recommand.
  # var.grp: phenotype variable name indicating case-control status,0 as control, 1 as case.
  # pAdjMethod:if pAdjMethod=FALSE, the function will not do mutiple testing adjustment. If pAdjMethod="fdr"/"BH"/"BY"/"holm"/"hochberg"/"hommel"/"bonferroni"/"BH"/"BY", the specific method will be used for adjusting p-value.
  # cutoff: if p value is smaller than the cutoff, there will be an edge between the two nodes.
  # pseudo_adjust_cutoff: if the value is true, pseudo probes will be used for adjusting the cutoff
  # output:
  # my_graph: obtained network in igraph object
  # my_dat: obtained netork as data frame with 3 columns: edge id, node_id1,node_id2
  mat1=Biobase::exprs(es)
  pDat1=Biobase::pData(es)
  fDat1=Biobase::fData(es)
  geneid=featureNames(es)

  status = pDat1[, c(var.grp)]
  pos0=which(status==0) # controls
  pos1=which(status==1) # cases
  nProbes=dim(mat1)[1] # number of genes/probes


  probes=1:nProbes
  r1_pos=c()
  r2_pos=c()
  r1=c()
  # for each probe, calculate test statistic between this probe and other probes
  for(i in 1:(nProbes-1)){
    temp1=rep(probes[i],time=nProbes-i)
    temp2=probes[(i+1):nProbes]
    r1_pos=c(r1_pos,temp1)
    r2_pos=c(r2_pos,temp2)
    x11=mat1[temp1,pos1]
    x10=mat1[temp1,pos0]
    z11=mat1[temp2,pos1]
    z10=mat1[temp2,pos0]
    if(i!=(nProbes-1)){
      x11<-tapply(x11,rep(1:base::nrow(x11),base::ncol(x11)),function(i)i)
      x10<-tapply(x10,rep(1:base::nrow(x10),base::ncol(x10)),function(i)i)
      z11<-tapply(z11,rep(1:base::nrow(z11),base::ncol(z11)),function(i)i)
      z10<-tapply(z10,rep(1:base::nrow(z10),base::ncol(z10)),function(i)i)
      switch (cor_method,
              st1 = {temp1=as.numeric(mapply(st1,x11,z11,x10,z10)[2,])},
              st2 = {temp1=as.numeric(mapply(st2,x11,z11,x10,z10)[2,])},
              st3 = {temp1=as.numeric(mapply(st3,x11,z11,x10,z10)[2,])},
              st4 = {temp1=as.numeric(mapply(st4,x11,z11,x10,z10)[2,])},
              st5 = {temp1=as.numeric(mapply(st5,x11,z11,x10,z10)[2,])},
              st6 = {temp1=as.numeric(mapply(st6,x11,z11,x10,z10)[2,])},
              stop("wrong cor_method!")
      )
    }else{
      switch(cor_method,
             st1={temp1=st1(x11,z11,x10,z10)$pval},
             st2={temp1=st2(x11,z11,x10,z10)$pval},
             st3={temp1=st3(x11,z11,x10,z10)$pval},
             st4={temp1=st4(x11,z11,x10,z10)$pval},
             st5={temp1=st5(x11,z11,x10,z10)$pval},
             st6={temp1=st6(x11,z11,x10,z10)$pval},
             stop('wrong cor_method')
      )

    }
    r1=c(r1,temp1)
  }

  # p-values for pseudo genes
  pvalPseudo = NULL
  # decide if a pair of probes is differentially correlated
  if(pseudo_adjust_cutoff==TRUE){
    #nr=nrow(es) # number of probes
    if(nPseudo >= nProbes/2)
    {
      stop("nPseudo must be smaller than half of number of genes!")
    }
    #nPseudo=50
    # randomly pick 'nPseudo' genes as pseudo genes
    pos.pick=sample(1:nProbes, size=nPseudo, replace=FALSE)
    mat1.pick=mat1[pos.pick,]
    fDat1.pick=fDat1[pos.pick,]
    nSubj=ncol(es) # number of subjects
    # for each pseudo gene gene, permute subject ids
    pos=sample(1:nSubj, size=nSubj, replace=FALSE)
    mat1.pick.s=mat1.pick[1,pos, drop=FALSE]
    for(i in 2:nPseudo){
      pos=sample(1:nSubj, size=nSubj, replace=FALSE)
      mat1.pick.s=rbind(mat1.pick.s,mat1.pick[i,pos])
    }
    # make sure to change the row names
    rownames(mat1.pick.s)=paste("perm", rownames(mat1.pick), sep=".")
    rownames(fDat1.pick)=rownames(mat1.pick.s)
    # also need to change the content of fDat.pick
    # suppose column 'ID' is the probe id
    fDat1.pick$ID=paste("perm", fDat1.pick$ID, sep=".")

    # add dat.pick.s to dat
    # mat2 include both original p1+p2 genes and nPseudo pseudo genes
    mat2=rbind(mat1,mat1.pick.s)

    fDat1$ID = NA
    fDat2=rbind(fDat1,fDat1.pick)
    #nrow=dim(mat1)[1] # number of genes
    #nPseudo=50
    #cutoff1=c()
    # for each pseudo gene, we obtain its test statistic for differential correlation test
    for(i in (nProbes-nPseudo+1):nProbes){
      # non-pseudo genes in cases
      x11=mat1[c(1:(nProbes-nPseudo)),pos1]
      # non-pseudo genes in controls
      x10=mat1[c(1:(nProbes-nPseudo)),pos0]
      # i-the pseudo genes in cases
      z11=mat1[rep(i,time=(nProbes-nPseudo)),pos1]
      # i-th pseudo genes in controls
      z10=mat1[rep(i,time=(nProbes-nPseudo)),pos0]
      x11<-tapply(x11,rep(1:base::nrow(x11),base::ncol(x11)),function(i)i)
      x10<-tapply(x10,rep(1:base::nrow(x10),base::ncol(x10)),function(i)i)
      z11<-tapply(z11,rep(1:base::nrow(z11),base::ncol(z11)),function(i)i)
      z10<-tapply(z10,rep(1:base::nrow(z10),base::ncol(z10)),function(i)i)
      # obtain p-value for testing differential correlation
      switch (cor_method,
              st1 = {temp1=as.numeric(mapply(st1,x11,z11,x10,z10)[2,])},
              st2 = {temp1=as.numeric(mapply(st2,x11,z11,x10,z10)[2,])},
              st3 = {temp1=as.numeric(mapply(st3,x11,z11,x10,z10)[2,])},
              st4 = {temp1=as.numeric(mapply(st4,x11,z11,x10,z10)[2,])},
              st5 = {temp1=as.numeric(mapply(st5,x11,z11,x10,z10)[2,])},
              st6 = {temp1=as.numeric(mapply(st6,x11,z11,x10,z10)[2,])},
              stop("wrong cor_method!")
      )

      # we expect the p-values obtained based on pseudo genes are large
      # because pseudo genes are non-differentially correlated.
      # 'cutoff' percentile of p-value for testing differential correlation
      #pvalPseudo=c(pvalPseudo, stats::quantile(temp1, prob=cutoff))
      pvalPseudo=c(pvalPseudo, temp1)
    }
    #alpha1=stats::quantile(pvalPseudo, prob=cutoff)
    alpha1 = min(pvalPseudo, na.rm=TRUE)
    r1Adj = r1
  }else{
    alpha1=cutoff
    r1Adj=stats::p.adjust(r1, pAdjMethod)
  }

  adjacency_matrix1=matrix(0,nrow = nProbes, ncol = nProbes)
  rownames(adjacency_matrix1) = geneid
  colnames(adjacency_matrix1) = geneid

  pvalMat=matrix(0,nrow = nProbes, ncol = nProbes)
  rownames(pvalMat) = geneid
  colnames(pvalMat) = geneid

  pAdjMat=matrix(0,nrow = nProbes, ncol = nProbes)
  rownames(pAdjMat) = geneid
  colnames(pAdjMat) = geneid

  j=1
  my_dat=data.frame()
  # if p-value < alpha1, then the 2 genes are differentially correlated
  for(i in 1:length(r1)){
    if(r1Adj[i]<alpha1){
      adjacency_matrix1[r1_pos[i],r2_pos[i]]=1
      adjacency_matrix1[r2_pos[i],r1_pos[i]]=1
      my_dat=rbind(my_dat,data.frame(edge_id=j,node_id1=r1_pos[i],node_id2=r2_pos[i]))
      j=j+1
    }
    pvalMat[r1_pos[i],r2_pos[i]] = r1[i]
    pvalMat[r2_pos[i],r1_pos[i]] = r1[i]

    pAdjMat[r1_pos[i],r2_pos[i]] = r1Adj[i]
    pAdjMat[r2_pos[i],r1_pos[i]] = r1Adj[i]

  }
  my_graph=igraph::graph_from_adjacency_matrix(as.matrix(adjacency_matrix1),
                                               mode = "undirected",
                                               weighted = NULL)
  res = list(graph = my_graph,
             network_dat = my_dat,
             pvalMat = pvalMat,
             pAdjMat = pAdjMat,
             pvalPseudo = pvalPseudo,
             alpha1 = alpha1)
  invisible(res)
}

generate_data<-function(n1=50,n2=60,p1=5,p2=100){
  ### The function is to generate expression level matrixes of control subjects and case subjects.
  ### X matrix is for diseased subjects with the default sample size n1=50. Z matrix is for control subjects with the default sample size n2=60.
  ### X is generated from multivariate normal distribution N (0, SigmaX), where SigmaX is a block matrix ((SigmaP1, 0), (0, SigmaP2)), sigmaP1 is the p1*p1 matrix and SigmaP2 is the p2*p2 matrix. Z is generated from multivariate normal distribution N (0, SigmaZ), where SigmaZ is a block matrix ((E_P1, 0), (0, SigmaP2)) and E_P1 is p1*p1 identity matrix.
  #n1=50
  #n2=60
  #p1=5
  #p2=100
  Sigma1=clusterGeneration::rcorrmatrix(d=p1)
  Sigma0=clusterGeneration::rcorrmatrix(d=p2)
  Ip1=diag(p1)
  SigmaX=as.matrix(Matrix::bdiag(Sigma1,Sigma0))
  SigmaZ=as.matrix(Matrix::bdiag(Ip1,Sigma0))

  # X - n1 x (p1+p2) matrix (cases)
  X=MASS::mvrnorm(n = n1, mu=numeric(p1+p2), Sigma=SigmaX, tol = 1e-6, empirical = FALSE)
  # Z - n2 x (p1+p2) matrix (controls)
  Z=MASS::mvrnorm(n = n2, mu=numeric(p1+p2), Sigma=SigmaZ, tol = 1e-6, empirical = FALSE)

  # rows are genes
  # columns are subjects
  dat = t(rbind(X, Z))

  geneid=paste("g", 1:nrow(dat), sep="")
  sid=paste("s", 1:ncol(dat), sep="")

  rownames(dat)=geneid
  colnames(dat)=sid

  # grp=0 indicates control; grp=1 indicates case
  grp=c(rep(1, n1), rep(0, n2))
  pDat = data.frame(sid=sid, grp=grp)
  rownames(pDat) = sid

  # memGenes=1 indicates the gene is differentially correlated with at least one another gene between cases and controls
  #  that is, correlation of this gene with at least one other gene in cases is different from that in controls.
  # memGenes=0 indicates the gene is non-differentially correlated with any other genes between cases and controls
  #  that is, correlation of this gene with other genes in cases is the same as that in controls.
  memGenes=c(rep(1, p1), rep(0, p2))
  fDat = data.frame(geneid=geneid, memGenes=memGenes)
  rownames(fDat)=geneid

  es = genEset(ex=dat, pDat=pDat, fDat = fDat, annotation = "")

  rownames(SigmaX)=geneid
  colnames(SigmaX)=geneid

  rownames(SigmaZ)=geneid
  colnames(SigmaZ)=geneid

  res = list(es=es, covCase = SigmaX, covCtrl = SigmaZ)

  invisible(res)
}
