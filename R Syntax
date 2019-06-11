data=read.csv("sell.csv",header=T,sep=",")
attach(data)

#Definisi tambahan
n = nrow(X)
p = ncol(X)
H=X%*%solve(t(X)%*%X)%*%t(X)

#Definisi variabel
Y=log(data$sell)
x1=log(data$lot)
x2=data$bdm
x3=data$fb
x4=data$sty
x5=data$drv
x6=data$rec
x7=data$ffin
x8=data$gwh
x9=data$ca
x10=data$gar
x11=data$reg
x0=matrix(c(1),nrow(data),ncol=1)
y_h=H%*%Y

#persamaan
X0=matrix(c(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11),nrow(data),11)
X=matrix(c(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11),nrow(data),12)
X1=matrix(c(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,y_h^2,y_h^3),nrow(data),14)


#Regresi
reg <- function(Y,X){
  betop <- (solve(t(X)%*%X))%*%(t(X)%*%Y)
  print(betop)
}
reg(Y,X)


#mencari model terbaik dengan metode Backward
BackStep <-function(X1,Y,p,n){
  betop1 <- (solve(t(X1)%*%X1))%*%(t(X1)%*%Y) 
  betop1
  H <- X1%*%(solve(t(X1)%*%X1))%*%(t(X1))
  e <- Y - (H%*%Y)
  dfe1<- n-p
  sse1 <- t(e)%*%e
  mse1 <- sse1/dfe1
  
  se1 <- sqrt(c(mse1)*diag(solve(t(X1)%*%X1)))
  t1 <- betop1/se1
  t1
  p_v <- 2*(1-pt(abs(t1),n-p))
  p_v
  hasil <- data.frame(t1,p_v)
  print(hasil)
}
BackStep(X,Y,2,511)
#dengan maksimum di x6 namun p-value x6 < alpha, x6 tidak perlu dikeluarkan dalam model

#multikolineritas
multi <- function(X){
  rka <- diag(solve(cor(X)))
  kesimpulana <- ifelse(rka>5,"ada multicolinear","tidak ada multicolinear")
  data <- data.frame(VIF=rka,Kesimpulan=kesimpulana)
  return(data)
}
multi(X0)

#mencari outlier
out <- function(X,Y){
  n <- nrow(X)
  p <- ncol(X)
  H <- X%*%(solve(t(X)%*%X))%*%(t(X))
  e <- Y - (H%*%Y)
  sse <- t(e)%*%e
  dfe <- n-p
  mse <- sse/dfe
  hii <- diag(H)
  kesimpulan <- ifelse(hii>2*p/n,"Outlier","bukan outlier")
  kesim <- c(kesimpulan)
  ti <- e*sqrt((n-p-1)/(c(sse)*(1-hii)-e^2))
  ttab <- qt(1-0.05/(2*n),n-p-1)
  kesi <- matrix(c(ifelse(abs(ti)>ttab,"OUTLIER","BUKAN OUTLIER")),n,1)
  dffits <- ti*sqrt(hii/(1-hii))
  ke <-c(ifelse(dffits>1,"INFLUENCE","TIDAK INFLUENCE"))
  DI <- (e^2/(p*c(mse)))*(hii/(1-hii)^2)
  ftab <- qf(0.5,p,n-p)
  k <- c(ifelse(DI>ftab,"INFLUENSIAL","TIDAK INFLUENSIAL"))
  nomer=matrix(c(1:546),546,1)
  data <- data.frame(No=nomer,Hii=hii,Outlier_di_X=kesim,ti,Outlier_di_Y=kesi,DFFITS=dffits,INFLUENSIAL=ke,COOKS_DISTANCE=DI,INFLUENSIALL=k)
  print(data)
}
outlier=out(X,Y)
#outlier di X
outX=outlier[outlier$Hii>2*p/n,]
#outlier di Y
outY=outlier[outlier$ti>qt(1-0.05/(2*n),n-p-1),]

nomer
#Evaluasi Reset uji ramsay
Fcon <- function(X0,X1,Y,q){
  n  <- nrow(X0)
  I <- diag(1,n,n)
  H <- X0%*%solve(t(X0)%*%X0)%*%t(X0)
  e <- (I-H)%*%Y
  sse0 <- 1-(t(e)%*%e)/(n-1)/var(Y)
  e1 <- (I-X1%*%solve(t(X1)%*%X1)%*%t(X1))%*%Y
  sse1 <- 1-(t(e1)%*%e1)/(n-1)/var(Y)
  Fhitung<- (sse1-sse0)/q/(1-sse1)/(n-p-q)
  p_val<- 1-pf(Fhitung,q,n-p-q)
  kesim <- ifelse(p_val<0.05,"Terjadi Kesalahan Spesifikasi","Tidak Terjadi Kesalahan Spesifikasi")
  d <- data.frame(Fhitung,Pvalue=p_val,Kesimpulan=kesim)
  print(d)
}
Fcon(X,X1,Y,2)

#Heteroskedastisitas
hetero <- function(X,Y){
  n <- nrow(X)
  p <- ncol(X)
  H <- X%*%(solve(t(X)%*%X))%*%(t(X))
  e <- Y - (H%*%Y)
  Yh <-e^2
  e_h <- Yh-(H%*%Yh)
  r_2 <- 1-(t(e_h)%*%e_h)/(n-1)/var(Yh)
  chisq <- n*r_2
  pv <- 1-pchisq(chisq,p)
  k <- ifelse(pv<0.05,"Tidak Heteroskedastisitas","Heteroskedastisitas")
  data <- data.frame(chisq,Pvalue=pv,Kesimpulan=k)
  print(data)
}
hetero(X,Y)


#Uji t
bk <- function(Y,X){
  n <- nrow(X)
  p <- ncol(X)
  betop <- (solve(t(X)%*%X))%*%(t(X)%*%Y)
  H <- X%*%(solve(t(X)%*%X))%*%(t(X))
  e <- Y - (H%*%Y)
  sse <- t(e)%*%e
  dfe <- n-p
  mse <- sse/dfe
  se <- sqrt(c(mse)*diag(solve(t(X)%*%X)))
  t <- betop/se
  p_v <- 2*(1-pt(abs(t),dfe))
  p_vkanan <- 1 - pt(abs(t),dfe)
  p_vkiri <- pt(abs(t),dfe)
  ke <- ifelse(p_v<0.05,"Ho Ditolak","Ho Diterima")
  pn <- ifelse(p_v<0.05,"Memberikan pengaruh yang berbeda","Tidak memberikan pengaruh yang berbeda")
  tabell <- data.frame(betop,t,se,p_v,p_vkanan,p_vkiri,Kesimpulan=ke,Penjelasan=pn)
  print(tabell)
}
bk(Y,X)

#anava
an <- function(Y,X){
  n <- nrow(X)
  p <- ncol(X)
  j <- matrix(c(1),n,n)
  H <- X%*%(solve(t(X)%*%X))%*%(t(X))
  e <- Y - (H%*%Y)
  ssto <- (t(Y)%*%Y)-((t(Y)%*%j%*%Y)/n)
  sse <- t(e)%*%e
  ssr <- ssto - sse
  dft <- n-1
  dfe <- n-p
  dfr <- p-1 
  mse <- sse/dfe
  msr <- ssr/dfr
  F <- msr/mse
  p_value <- 1-pf(F,dfr,dfe)
  SV <- c("Regresi","Erorr","Total")
  Df <- c(dfr,dfe,dft)
  Ss <- c(ssr,sse,ssto)
  Ms <- c(msr,mse,"-")
  f <- c(F,"-","-")
  pv <- c(p_value,"-","-")
  ke <- ifelse(p_value<0.05,"Ho Ditolak","Ho Diterima")
  kes <- c(ke,"-","-")
  pe <- ifelse(p_value<0.05,"Terdapat pengaruh yang diberikan variabel","Tidak terdapat pengaruh yang diberikan variabel")
  pen <- c(pe,"-","-")
  anava <- data.frame(SV,Df,Ss,Ms, F=f, "P value"=pv, Kesimpulan = kes, Penjelasan= pen)
  print(anava)
} 
an(Y,X)
