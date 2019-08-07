# rm(list=ls())

library(rARPACK)
library(magrittr)
library(Matrix)

hash = function(el, nc){
  el[,1] + el[,2]/nc
}

makeL = function(A, tau= NULL){
  rs = rowSums(A)
  if(is.null(tau)) tau = mean(rs)
  D = Diagonal(length(rs),1/sqrt(rs + tau))
  D%*%A%*%D
}

makeL_rect = function(A, tau= NULL){
  rs = rowSums(A)
  cs = colSums(A)
  if(is.null(tau)) tau = c(mean(rs), mean(cs))
  
  Dr = Diagonal(nrow(A),1/sqrt(rs + tau[1]))
  Dc = Diagonal(ncol(A),1/sqrt(cs + tau[2]))
  Dr%*%A%*%Dc
}


poissonBoot_rect = function(A, s = NULL){
  if(is.null(s)) s = makeL_rect(A) %>% svds(150)
  Asamp = A
  Asamp@x = as.numeric(rpois(length(A@x), 1))
  Lsamp = makeL_rect(Asamp)
  
  bootVals = apply(rbind(s$u, s$v), 2, function(x){
    uu = x[1:nrow(s$u)];
    vv = x[-(1:nrow(s$u))];
    return(as.numeric(t(uu)%*%Lsamp%*%vv))
  }) 
  s.samp = svds(Lsamp, 150)
  bootCorU = (t(s.samp$u)%*%s$u)^2 %>% rowSums 
  
  
  Asamp2 = Asamp
  Asamp2@x = as.numeric(rpois(length(Asamp@x), Asamp@x))
  Lsamp2 = makeL_rect(Asamp2)
  
  s.samp2 = svds(Lsamp2, 150)
  bootCorU2 = (t(s.samp2$u)%*%s.samp$u)^2 %>% rowSums 
  
  
  
  # superm = apply(s.samp$u[,-1], 1, function(x){
  #   sample(x)*sample(c(-1,1), length(x), replace = T)
  # }) %>% t

  # kk = ncol(s.samp$u)-1
  # superm = apply(s.samp$u[,-1],1, function(x){
  #   rnorm(kk, rep(0,kk),sd = sd(x))
  # }) %>% t
  # # svperm = t(superm)%*%Lsamp %>% t
  # # superm = Lsamp%*%svperm
  # # vp = svd(svperm)$u
  # 
  # superm = matrix(rnorm(n*150, 0, ))
  # 
  # # superm = s.samp$u[sample(nrow(s.samp$u))]
  # bootCorUperm = (t(superm)%*%s$u)^2 %>% rowSums
  # quantile(bootCorUperm, .95)
  # 
  return(list(val.boot = bootVals, coru.boot = bootCorU, coru.perm = quantile(bootCorU2, .95)))
  
  # bootCorV = (t(s.samp$v)%*%s$v)^2 %>% rowSums 
  # return(list(val.boot = bootVals, coru.boot = bootCorU, corv.boot = bootCorV))
}

poissonBoot_symmetric = function(A, ei = NULL){
  if(is.null(ei)) ei = makeL(A) %>% eigs(150)
  Asamp = A
  Asamp@x = as.numeric(rpois(length(A@x), 1))
  Asamp = Asamp+t(Asamp)
  Lsamp = makeL(Asamp)
  
  bootVals = apply(ei$vec, 2, function(x){
    return(as.numeric(t(x)%*%Lsamp%*%x))
  }) 
  ei.samp = eigs(Lsamp, 150)
  bootCor = (t(ei.samp$vectors)%*%ei$vectors)^2 %>% rowSums 
  return(list(val.boot = bootVals, cor.boot = bootCor))
}



edgeBoot_rect = function(A){
  
  TT = as(A, "dgTMatrix")
  el = cbind(i= TT@i, j= TT@j, x = TT@x)
  
  m = nrow(el)
  samp = sample(m, .05*m)
  
  el.samp = el[samp,]
  el.fit = el[-samp, ]

  
  A.test = spMatrix(nrow(A), ncol(A), i = el.samp[,1] +1, j= el.samp[,2] +1, x= rep(1, nrow(el.samp)))
  A.test@x[A.test@x>1] = 1
  
  A.fit= spMatrix(nrow(A), ncol(A), i = el.fit[,1]+1, j= el.fit[,2]+1 , x= rep(1, nrow(el.fit)))
  A.fit@x[A.fit@x>1] = 1
  
  L.fit = makeL_rect(A.fit)
  L.test = makeL_rect(A.test)
  
  s = svds(L.fit, 150)
  uv = rbind(s$u, s$v)
  bootVals = apply(uv, 2, function(x){
    uu = x[1:nrow(s$u)];
    vv = x[-(1:nrow(s$u))];
    return(as.numeric(t(uu)%*%L.test%*%vv))
  }) 
  
  # uvperm = apply(uv,1, function(x){
  #   sample(x)*sample(c(-1,1), length(x), replace = T)
  # }) %>% t
  # 
  kk = ncol(uv)
  uvperm = apply(uv,1, function(x){
    rnorm(kk, rep(0,kk),sd = abs(x))
  }) %>% t

  bootValsPerm = apply(uvperm, 2, function(x){
    uu = x[1:nrow(s$u)];
    vv = x[-(1:nrow(s$u))];
    return(as.numeric(t(uu)%*%L.test%*%vv))
  }) 
  return(list(d.fit = s$d, d.test = bootVals, d.perm = 2*sd(bootValsPerm)))
}


edgeBoot_symmetric = function(A){
  
  TT = as(A, "dgTMatrix")
  el = cbind(i= TT@i, j= TT@j, x = TT@x)
  
  m = nrow(el)
  samp = sample(m, .05*m)
  
  el.samp = el[samp,]
  h.samp = c(hash(el.samp[,2:1], ncol(TT)), hash(el.samp[,1:2], ncol(TT)))
  h = hash(el, ncol(TT))
  removeThese = match(h.samp, h)
  el.fit = el[-removeThese, ]
  
  A.test = spMatrix(nrow(A), ncol(A), i = c(el.samp[,1], el.samp[,2])+1, j= c(el.samp[,2], el.samp[,1]) +1, x= rep(1, 2*nrow(el.samp)))
  A.test@x[A.test@x>1] = 1
  
  A.fit= spMatrix(nrow(A), ncol(A), i = c(el.fit[,1], el.fit[,2])+1, j= c(el.fit[,2], el.fit[,1])+1 , x= rep(1, 2*nrow(el.fit)))
  A.fit@x[A.fit@x>1] = 1
  
  L.fit = makeL(A.fit)
  L.test = makeL(A.test)
  
  ei.fit = eigs(L.fit, 150)
  
  bootVal = apply(ei.fit$vec, 2, function(x){
    return(as.numeric(t(x)%*%L.test%*%x))
  })
  
  
  return(list(val.fit = ei.fit$val, val.test = bootVal))
}







makeL3 = function(A){
  # For L = the regularized graph Laplacian, make [L*L].*L  <- matlab notation
  rs = rowSums(A)
  D = Diagonal(nrow(A), 1/sqrt(rs + mean(rs)))
  L = D%*%A; L = L%*%D
  # L  <- as(L, "symmetricMatrix")
  LL = L%*%L
  
  # the most difficult part is the elementwise multiplication in [L*L].*L.
  # I talked to Doug.  He suggested this using the hash function defined above.
  
  # now, convert them both to "triple" format.  So, each nonzero entry of LL makes one additional element in (tLL@i, tLL@j, tLL@x). the first element is the row, then the column, then the actual value.
  tLL = as(LL, Class = "TsparseMatrix")
  tL =  as(L, Class = "TsparseMatrix")
  
  # first we need to find where the elements of LL and L are both nonzero...
  hl = hash(cbind(tL@i, tL@j), ncol(L))
  h2 = hash(cbind(tLL@i, tLL@j), ncol(L))
  keep = match(h2, hl)
  
  
  include = complete.cases(keep)
  keepi = keep[include]
  # this does the multiplication for all of the places where LL and L both have nonzero elements.
  trip = cbind(tLL@i[include],tLL@j[include], tLL@x[include]*tL@x[keepi])
  L3 = spMatrix(nrow=nrow(L), ncol = ncol(L), trip[,1]+1, trip[,2]+1, trip[,3]*1000)
  diag(L3) = 0
  return(L3)
}
