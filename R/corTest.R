
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
  return(1-pchisq(T,1))
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
  w1<-n1*(x1-x1_median)*(z1-z1_median)/c1
  w0<-n0*(x0-x0_median)*(z0-z0_median)/c0
  u=sum(w1)*p2-sum(w0)*p1
  w=c(w1,w0)
  u_var=p1*p2*var(w)*(n1+n0-1)
  T=u^2/u_var
  return(1-pchisq(T,1))
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
  return(1-pchisq(T,1))
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
  c1<-n1*(n1^2-1)
  c0<-n0*(n0^2-1)
  rank_x1<-rank(x1,ties.method = "average")
  rank_z1<-rank(z1,ties.method = "average")
  rank_x0<-rank(x0,ties.method = "average")
  rank_z0<-rank(z0,ties.method = "average")
  w1=n1^2-1-6*(rank_x1-rank_z1)^2/c1
  w0=n0^2-1-6*(rank_x0-rank_z0)^2/c0
  u=sum(w1)*p2-sum(w0)*p1
  w=c(w1,w0)
  u_var=p2*p1*var(w)*(n1+n0-1)
  T=u^2/u_var
  return(1-pchisq(T,1))
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
  c1 <- n1 * (n1^2 - 1)
  c0 <- n0 * (n0^2 - 1)
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

  c1<-n1*(n1^2-1)
  c0<-n0*(n0^2-1)
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
}
