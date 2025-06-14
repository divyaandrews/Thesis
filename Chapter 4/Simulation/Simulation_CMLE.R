rm(list=ls(all=TRUE))

library(GA)
x<-c()
M<-c()

mse1<-function(x){
  return(mean((x-1)^2)) 
}
mse2<-function(x){
  return(mean((x-0.4)^2)) 
}
mse3<-function(x){
  return(mean((x-0.2)^2)) 
}
mse4<-function(x){
  return(mean((x-0.55)^2)) 
}



M[1] <- 2
x[1]<- 1

rnoge <-function(alpha,beta){
  
  u = runif(1)
  if(u < alpha) return(0)
  else return(1+rgeom(1,1/beta))
  
}

#set.seed(1234)
ingdata<-function(n,alpha0, alpha1, beta1,phi){
  for(i in 2:n){
    M[i]<- alpha0 + alpha1*x[i-1] + beta1*M[i-1]
    x[i]<- rnoge(phi,M[i])
  }
  return(list(x=x, M = M))
}

indf<- function(y){
  ifelse(y==0,1,0) 
}


ll<- function(par){
  alpha0 = par[1]; alpha1 = par[2]; beta1 = par[3]; phi = par[4];
  x = x_t
  M = lam
  f = numeric(n)
  for(i in 2:n){
    f[i]<- indf(x[i])*log(phi)+ 
      (1-indf(x[i]))*(log(1-phi)+ (x[i]-1)*log(1-(1/M[i])) - log(M[i]))
  }
  return(sum(f, na.rm=TRUE))
}


nlist<-c(50,200,500)
iterr = 100

#SE1 = matrix(NA, length(nlist), 4)
m1=matrix(NA,length(nlist),4)
mse_1=matrix(NA,length(nlist),4)

for(k in 1:length(nlist)){
  n = nlist[k]
  stval1= matrix(NA, iterr, 4)
  for(j in 1:iterr){
    x_t <- ingdata(n,1, 0.4, 0.2,0.55)$x
    lam <- ingdata(n,1, 0.4, 0.2,0.55)$M
    GA1 <- ga(type = "real-valued", fitness = ll, lower = c(0.01,0.01,0.01,0.01), upper = c(2,1,1,1),
              maxiter = 100, popSize = 500, monitor = FALSE)
    estml1 <- suppressWarnings(optim(summary(GA1)$solution[1,], ll, method="L-BFGS-B", lower=rep(0.01,4), upper=c(2,0.9999,0.9999,0.99)))
    stval1[j,] <- estml1$par
  }
  
  
  #SE1[k,] = apply(stval1, 2, sd)
  
  m1[k,] = apply(stval1, 2, mean)
  
  mse_1[k,1] = mse1(stval1[,1])
  mse_1[k,2] = mse2(stval1[,2])
  mse_1[k,3] = mse3(stval1[,3])
  mse_1[k,4] = mse4(stval1[,4])
  
}

#SE1

df1<-data.frame(m1[1,],mse_1[1,],m1[2,],mse_1[2,],m1[3,],mse_1[3,],m1[4,],mse_1[4,])
para<-c("alpha_0=1", "alpha_1= 0.4", "beta_1=0.2", "phi=0.55")
samp<-c("para",rep("n=50",2),rep("n=100",2),rep("n=200",2),rep("n=500",2))
df2<-cbind(para,df1)
df2<-rbind(samp,df2)
#colnames(df2) <- c("Sample Size", "Parameter", "CMLE", "MSE")
dataset_names <- list('Sheet1'= df2) 

#export each data frame to separate sheets in same Excel file
openxlsx::write.xlsx(dataset_names, file = '~/Documents/NEW_1_0.4_0.2_0.55.xlsx')
