rm(list=ls())

source("https://raw.githubusercontent.com/karlrohe/fastRG/master/fastRDPG.R")
source('~/Dropbox/GettingCash/nsf2018/code/bootstrapFunctions.R')


set.seed(1)
k = 100
n = 10000
d = 20000


Xhard = matrix(rgamma(n = n*k, shape = .15, rate = .1), ncol = k)*cbind(rep(2,n),matrix(rbinom(n = n*(k-1), prob = 3/k,size = 1), ncol = k-1))
Y = matrix(rgamma(n = d*k, shape = .15, rate = .1), ncol = k)*cbind(rep(2,d),matrix(rbinom(n = d*(k-1), prob = 10/k,size = 1), ncol = k-1))
Ahard = fastRG(Xhard,Y = Y,S = diag(rep(1,k)), avgDeg = 30)


Xmid = matrix(rgamma(n = n*k, shape = .15, rate = .1), ncol = k)*cbind(rep(1.4,n),matrix(rbinom(n = n*(k-1), prob = 3/k,size = 1), ncol = k-1))
Y = matrix(rgamma(n = d*k, shape = .15, rate = .1), ncol = k)*cbind(rep(1.4,d),matrix(rbinom(n = d*(k-1), prob = 10/k,size = 1), ncol = k-1))
Amid = fastRG(Xmid,Y = Y,S = diag(rep(1,k)), avgDeg = 30)

Xeasy = matrix(rgamma(n = n*k, shape = .15, rate = .1), ncol = k)*cbind(rep(1,n),matrix(rbinom(n = n*(k-1), prob = 3/k,size = 1), ncol = k-1))
Y = matrix(rgamma(n = d*k, shape = .15, rate = .1), ncol = k)*cbind(rep(1,d),matrix(rbinom(n = d*(k-1), prob = 10/k,size = 1), ncol = k-1))
Aeasy = fastRG(Xeasy,Y = Y,S = diag(rep(1,k)), avgDeg = 30)





makeScree = function(A, title){
  
  
  A@x[A@x>1] = 1
  # isSymmetric(A)
  
  s = makeL_rect(A) %>% svds(150)
  
  plot(s$d[-1], main = title, ylab = "eigenvalue", xlab = "2:150")
}

pdf("scree.pdf", height = 2, width = 5)
par(mfrow = c(1,3))
makeScree(Aeasy, title = "Graph 1")
makeScree(Amid, title = "Graph 2")
makeScree(Ahard, title = "Graph 3")
dev.off()


makeSummary = function(A, X,  title){
  
  A@x[A@x>1] = 1
  # isSymmetric(A)
  
  s = makeL_rect(A) %>% svds(150)
  
  # u = svd(X)$u
  trueCor = (t(s$u)%*%X)^2 %>% rowSums 
  
  # pboot = poissonBoot_rect(A, s)
  eboot = edgeBoot_rect(A)
  
  
  # highLight = c(.2,0)
  # redd = trueCor<highLight[1] & trueCor>highLight[2]
  
  # redd = (eboot$d.test[-1] < eboot$d.perm)
  # redd = (pboot$coru.boot < pboot$coru.perm) & (eboot$d.test < eboot$d.perm)
  
  
  plot((eboot$d.test[-1]), main = paste(title, "\nBootstrapped eigenvalues"), xlab = "Ordered eigenvalues, 2:150", ylab = "")
  # plot((pboot$coru.boot)[-1],(eboot$d.test[-1]), main = "BootDiagnostics")
  # points((pboot$coru.boot)[redd],(eboot$d.test[redd]),pch = "|")
  lines(c(-10,10^6), eboot$d.perm*c(1,1), col = "red")
  # lines(pboot$coru.perm*c(1,1), c(eboot$d.perm, -10), col = "red")
  
  plot(2:150, trueCor[-1], main = "Unobserved:\nis eigenvector just noise?", ylab = "v' E(A) v")
  # points(which(redd), trueCor[redd], pch = "|")
  # points(y = rep(0, 149)[!redd],x = which(!redd), pch  = "|")
  # points(y = rep(-.1, 149)[redd],x = which(redd), pch  = "|")
  # abline(highLight[1], 0)
  # abline(highLight[2],0)
}


pdf("summaryPlots.pdf", height= 3.5, width = 9)
par(mfcol = c(2,3))
makeSummary(Aeasy,Xeasy, title = "Graph 1")
makeSummary(Amid,Xmid, title = "Graph 2")
makeSummary(Ahard,Xhard, title = "Graph 3")
dev.off()








Xhard = matrix(rgamma(n = n*k, shape = .1, rate = .1), ncol = k)*cbind(rep(2,n),matrix(rbinom(n = n*(k-1), prob = 4/k,size = 1), ncol = k-1))
Ahard = fastRG(Xhard,S = diag(rep(1,k)), avgDeg = 50)
makeSummary(Ahard,Xhard, title = "Graph 3")

L3 = makeL3(Ahard) %>% makeL
s3 = svds(L3, 150)
s = makeL(Ahard) %>% svds(150)
plot(s$d[-1], s3$d[-1])
makeSummary(L3,Xhard, title = "Graph 3")
