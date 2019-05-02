
#library(fossil)
#library(cvxbiclustr)
#library(cvxclustr)
# library(MASS)
# library(Matrix)

k_row = k_col = 5

sylvester= function(A, B, C, eps = 0.0001){
  
  library(MASS)
  library(Matrix)
  
  A1 = Schur(A)
  Q1 = A1$Q; R1 = A1$T
  
  A2 = Schur(B)
  Q2 = A2$Q; R2 = A2$T 
  C = t(Q1) %*% C %*% Q2
  
  Rsq = R1 * R1
  I = diag(dim(A)[1])
  
  k = 1
  n = dim(R2)[1]
  
  while(k < n + 1){
    if(k < n){
      if(abs(R2[k+1, k]) < eps){
        left = R1 + R2[k,k] * I
        right = C[,k]
        temp = matrix(0, dim(X)[1],1)
        if(k == 1){
          temp = temp
        }else{
          temp = (X[,1:(k-1)]) %*% matrix(R2[1:(k-1),k],k-1,1)
        }
        temp = matrix(temp, dim(C)[1],1)
        X[,k] = ginv(left) %*% (right - temp)
        # mytry = myTryCatch(solve(left))
        #if(is.null(mytry$error) == 0){er = c(er,tt)}
        k = k+1
      }else{
        r11 = R2[k,k]; r12 = R2[k, k+1]; r21 = R2[k+1, k]; r22 = R2[k+1, k+1]
        temp2 = matrix(0, dim(X)[1],1);temp3=matrix(0, dim(X)[1],1)
        if(k == 1){
          temp2 = temp2
          temp3 = temp3
        }else{
          temps = X[,1:(k-1)] %*% matrix(R2[1:(k-1),k:(k+1)],k-1,2)
          temp2 = temps[,1]
          temp3 = temps[,2]
        }
        b1 = C[,k] - temp2 
        b2 = C[,k+1] - temp3
        b1_prime = R1 %*% b1 + r22 * b1 - r21 * b2
        b2_prime = R1 %*% b2 + r11 * b2 - r12 * b1
        b_prime = matrix(0, dim(X)[1],2)
        b_prime[,1] = b1_prime; b_prime[,2] = b2_prime
        X[,k:(k+1)] = ginv(R1 %*% R1 + (r11 + r22) * R1 +
                             (r11*r22 - r12*r21) * I) %*% b_prime 
        k = k+2
      }
    }else{
      if(abs(R2[1, k]) > eps){
        left = R1 + R2[k,k] * I
        right = C[,k]
        temp = matrix(0, dim(X)[1],1)
        if(k == 1){
          temp = temp
        }else{
          temp = (X[,1:(k-1)]) %*% matrix(R2[1:(k-1),k],k-1,1)
        }
        temp = matrix(temp, dim(C)[1],1)
        X[,k] = ginv(left) %*% (right - temp)
        k = k+1
      }else{
        R22 = R2
        R22 = cbind(R2, rep(0,dim(R2)[1]))
        R22 = rbind(R22,rep(0,dim(R2)[1]+1))
        r11 = R22[k,k]; r12 = R22[k, k+1]; r21 = R22[k+1, k]; r22 = R22[k+1, k+1]
        temp2 = matrix(0, dim(X)[1],1);temp3=matrix(0, dim(X)[1],1)
        
        if(k == 1){
          temp2 = temp2
          temp3 = temp3
        }else{
          temps = X[,1:(k-1)] %*% matrix(R22[1:(k-1),k:(k+1)],k-1,2)
          temp2 = temps[,1]
          temp3 = temps[,2]
        }
        
        b1 = C[,k] - temp2 
        b2 = - temp3
        b1_prime = R1 %*% b1 + r22 * b1 - r21 * b2
        b2_prime = R1 %*% b2 + r11 * b2 - r12 * b1
        b_prime = matrix(0, dim(X)[1],2)
        b_prime[,1] = b1_prime; b_prime[,2] = b2_prime
        GOD = ginv(R1 %*% R1 + (r11 + r22) * R1 +
                     (r11*r22 - r12*r21) * I) %*% b_prime 
        X[,k] = GOD[,1]
        k = k+2
      }
    }
  }
  return(Q1 %*% X %*% t(Q2))
}

L_num = function(n){
  t = matrix(0, n*(n-1)/2,2)
  count = 0
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      count = count + 1
      t[count,1]= i; t[count,2]=j
    }
  }
  LL = data.frame(l1 = t[,1],l2 = t[,2])
  LL = LL[order(LL$l1),]
  rownames(LL) = NULL
  return(LL)
}

elk = function(n,p){
  # n,l
  count = 0
  el1 = el2 = matrix(0,n,n*(n-1)/2)
  
  for(i in (n-1):1){
    temp = matrix(0,n,i)
    temp[n-i,] = 1
    el1[,(count+1):(count+i)] = temp
    el2[(n-i+1):n,(count+1):(count+i)] = diag(1,i,i)
    count = count +i
  }  
  
  # p,k
  count = 0
  ek1 = ek2 = matrix(0,p,p*(p-1)/2)
  for(i in (p-1):1){
    temp = matrix(0,p,i)
    temp[p-i,] = 1
    ek1[,(count+1):(count+i)] = temp
    ek2[(p-i+1):p,(count+1):(count+i)] = diag(1,i,i)
    count = count +i
  }
  return(list(el1 = el1, el2 = el2, ek1 = ek1, ek2 = ek2))
}

bi_ADMM_order = function(X, nu1, nu2, 
                         gamma_1, gamma_2, 
                         kk=5, phi=0.5,niter = 1000,tol = 0.1,output = 1){
  
  n=dim(X)[1]; p=dim(X)[2]
  
  n2 = n*(n-1)/2
  p2 = p*(p-1)/2
  
  elks = elk(n,p)
  el1 = elks$el1
  el2 = elks$el2
  ek1 = elks$ek1
  ek2 = elks$ek2
  
  k_row = k_col = kk
  
  p <- nrow(X)
  n <- ncol(X)
  w_row <- kernel_weights(t(X), phi/n)
  w_col <- kernel_weights(X, phi/p)
  w_row <- knn_weights(w_row, k_row, p)
  w_col <- knn_weights(w_col, k_col, n)
  w_row <- w_row/sum(w_row)
  w_col <- w_col/sum(w_col)
  w_row <- w_row/sqrt(n)
  w_col <- w_col/sqrt(p)
  
  w_l = w_row; u_k = w_col
  
  n=dim(X)[1]; p=dim(X)[2]
  
  A = matrix(0,n,p)
  v = matrix(0,p,n2)
  z = matrix(0,n,p2)
  lambda_1 = matrix(0,p,n2)
  lambda_2 = matrix(0,n,p2)
  
  for(iter in 1: niter){
    
    A_old = A; v_old = v; z_old = z; lambda_1_old = lambda_1; lambda_2_old = lambda_2
    
    # update A
    
    En = diag(0:(n - 1)) + diag((n - 1):0) - matrix(1, n, n) + diag(1, n, n)
    Ep = diag(0:(p - 1)) + diag((p - 1):0) - matrix(1, p, p) + diag(1, p, p)
    
    M = diag(1,n,n) + nu1 * En
    
    N = nu2 * Ep
    
    lv = lambda_1+ nu1 * v
    lz = lambda_2 + nu2 * z
    
    C2 = (el1-el2) %*% t(lv)
    C3 = lz %*% t(ek1-ek2)
    C = X +  C2 + C3  
    
    A = sylvester(M,t(N),C)
    
    al1 = t(A) %*% el1; al2 = t(A) %*% el2
    ak1 = A %*% ek1; ak2 = A %*% ek2
    
    
    # update vz
    sigma_1 = gamma_1 * w_l/nu1
    vtemp = al1 - al2 - 1/nu1 * lambda_1
    
    temp1 = ifelse((1 - sigma_1/apply(vtemp^2,2,sum)) < 0, 0,1 - sigma_1/apply(vtemp^2,2,sum))
    temp2 = matrix(temp1,dim(vtemp)[1],dim(vtemp)[2], byrow = TRUE) * vtemp 
    v = temp2
    
    ztemp = ak1 - ak2 - 1/nu2 * lambda_2
    sigma_2 = gamma_2 * u_k/nu2
    
    temp3 = ifelse((1 - sigma_2/apply(ztemp^2,2,sum)) < 0, 0,1 - sigma_2/apply(ztemp^2,2,sum))
    temp4 = matrix(temp3,dim(ztemp)[1],dim(ztemp)[2], byrow = TRUE) * ztemp 
    
    z = temp4
    
    # update lambda
    lambda_1 = lambda_1 + nu1 * (v - al1 + al2)
    
    # update lambda 2
    lambda_2 = lambda_2 + nu2 * (z - ak1 + ak2)
    if(output == 1){
      print('iter')
      print(iter)
      
      print(paste('A',sum(abs(A - A_old))))
      print(paste('v',sum(abs(v - v_old))))
      print(paste('z',sum(abs(z -z_old))))
      print(paste('l',sum(abs(lambda_1 - lambda_1_old))))
      print(paste('2',sum(abs(lambda_2 - lambda_2_old))))
    }
    
    
    # whether coverage
    if(sum(abs(A - A_old)) < tol & 
       sum(abs(v - v_old)) < 0.01& 
       sum(abs(z - z_old)) < 0.01 & 
       sum(abs(lambda_1 - lambda_1_old)) < 0.01 & 
       sum(abs(lambda_2 - lambda_2_old)) <0.01){
      return(list(A = A, v = v, z = z, 
                  lambad_1 = lambda_1, lambad_2 = lambda_2, niter = iter))
      break
    }
  }
  
  if(iter == niter){
    print(paste('not converge within',iter, 'times'))
    return(list(A = A, v = v, z = z, 
                lambad_1 = lambda_1, lambad_2 = lambda_2, niter = iter))
  }
}

dist_weight = function (X, phi, dist.type, p){
  dist_X <- as.numeric(dist(t(X), method = dist.type))
  exp(-phi * dist_X^p)
}

create_adjacency <- function(V,Phi) {
  differences <- apply(V,2,FUN=function(x) {norm(as.matrix(x),'f')})
  connected_ix <- which(differences == 0)
  n <- ncol(Phi)
  m <- length(connected_ix)
  A <- Matrix(0, nrow = n, ncol = n, sparse = TRUE)
  
  if (m > 0) {
    ix <- integer(m)
    jx <- integer(m)
    for (i in 1:m) {
      ix[i] <- which(Phi[connected_ix[i],]==1)
      jx[i] <- which(Phi[connected_ix[i],]==-1)
    }
    A[(jx-1)*n + ix] <- 1
  }
  return(A)
} 

find_clusters <- function(A) {
  G <- graph.adjacency(A, mode = 'upper')
  n <- nrow(A)
  node_seen <- logical(n)
  cluster <- integer(n)
  k <- 1
  for (i in 1:n) {
    if (!node_seen[i]) {
      connected_set <- graph.bfs(G, root=i, unreachable = FALSE)$order
      node_seen[connected_set] <- TRUE
      cluster[connected_set] <- k
      k <- k + 1
    }
  }
  nClusters <- k - 1
  size <- integer(nClusters)
  for (j in 1:nClusters) {
    size[j] <- length(which(cluster == j))
  }
  return(list(cluster=cluster, size=size))
}

## Clusterpath preprocessing
tri2vec <- function(i,j,n) {
  return(n*(i-1) - i*(i-1)/2 + j -i)
}

vec2tri <- function(k,n) {
  i <- ceiling(0.5*(2*n-1 - sqrt((2*n-1)^2 - 8*k)))
  j <- k - n*(i-1) + i*(i-1)/2 + i
  return(as.matrix(cbind(i,j)))
}

knn_weights <- function(w,k,n) {
  i <- 1
  neighbors <- tri2vec(i,(i+1):n,n)
  keep <- neighbors[sort(w[neighbors],decreasing=TRUE,index.return=TRUE)$ix[1:k]]
  for (i in 2:(n-1)) {
    group_A <- tri2vec(i,(i+1):n,n)
    group_B <- tri2vec(1:(i-1),i,n)
    neighbors <- c(group_A,group_B)
    knn <- neighbors[sort(w[neighbors],decreasing=TRUE,index.return=TRUE)$ix[1:k]]
    keep <- union(knn,keep)
  }
  i <- n
  neighbors <- tri2vec(1:(i-1),i,n)
  knn <- neighbors[sort(w[neighbors],decreasing=TRUE,index.return=TRUE)$ix[1:k]]
  keep <- union(knn,keep)
  if (length(keep) > 0)
    w[-keep] <- 0
  return(Matrix(data=w,ncol=1,sparse=TRUE))
}

kernel_weights <- function(X,phi=1) {
  storage.mode(X) <- "double"
  p <- as.integer(nrow(X))
  n <- as.integer(ncol(X))
  phi <- as.double(phi)
  w <- double(n*(n-1)/2)
  sol <- .C('kernel_weights',X=X,p=p,n=n,phi=phi,w=w)
  return(weights=sol$w)
}

cluster_assign = function(X, kk, result, method = 'ADMM'){
  phi = 0.5
  # calculate the weights
  k_row = k_col = kk
  p <- nrow(X)
  n <- ncol(X)
  w_row <- kernel_weights(t(X), phi/n)
  w_col <- kernel_weights(X, phi/p)
  w_row <- knn_weights(w_row, k_row, p)
  w_col <- knn_weights(w_col, k_col, n)
  w_row <- w_row/sum(w_row)
  w_col <- w_col/sum(w_col)
  w_row <- w_row/sqrt(n)
  w_col <- w_col/sqrt(p)
  
  nrow = dim(X)[1]; ncol = dim(X)[2]
  A = result$A
  
  if(method == 'ADMM'){
    V1 = result$v[,which(w_row!=0)]
    Z1 = result$z[,which(w_col!=0)]
  }
  
  if(method == 'AMA'){
    V1 = matrix(0,ncol,nrow*(nrow-1)/2)
    Z1 = matrix(0,nrow,ncol*(ncol-1)/2)
    
    count = 1
    for(i in 1:(nrow-1)){
      for(j in (i+1):nrow){
        V1[,count] = A[i,] - A[j,]
        count = count + 1
      }
    }
    
    count = 1
    for(i in 1:(ncol-1)){
      for(j in (i+1):ncol){
        Z1[,count] = A[,i] - A[,j]
        count = count + 1
      }
    }
    V1 = V1[,which(w_row!=0)]
    Z1 = Z1[,which(w_col!=0)]
  }
  
  nrow = dim(X)[1]; ncol = dim(X)[2]
  w = as.vector(w_row)
  P = which(w!=0)
  P = vec2tri(P,nrow)
  E = matrix(0,nrow,dim(P)[1])
  for(i in 1:dim(P)[1]){
    E[P[i,1],i] = 1
    E[P[i,2],i] = -1
  }
  E1 = E
  
  ncol = dim(X)[2]
  w = as.vector(w_col)
  P = which(w!=0)
  P = vec2tri(P,ncol)
  E = matrix(0,ncol,dim(P)[1])
  for(i in 1:dim(P)[1]){
    E[P[i,1],i] = 1
    E[P[i,2],i] = -1
  }
  E2 = E
  
  create_adjacency3 <- function(V,Phi) {
    differences <- apply(V,2,FUN=function(x) {norm(as.matrix(x),'f')})
    connected_ix <- which(abs(differences) == 0)
    n <- ncol(Phi)
    m <- length(connected_ix)
    A <- Matrix(0, nrow = n, ncol = n, sparse = TRUE)
    
    if (m > 0) {
      ix <- integer(m)
      jx <- integer(m)
      for (i in 1:m) {
        ix[i] <- which(Phi[connected_ix[i],]==1)
        jx[i] <- which(Phi[connected_ix[i],]==-1)
      }
      A[(jx-1)*n + ix] <- 1
    }
    return(A)
  }
  create_adjacency4 <- function(V,Phi) {
    V = ifelse(abs(V) < 1e-04, 0, V)
    differences <- apply(V,2,FUN=function(x) {norm(as.matrix(x),'f')})
    connected_ix <- which(abs(differences) < 0.03)
    n <- ncol(Phi)
    m <- length(connected_ix)
    A <- Matrix(0, nrow = n, ncol = n, sparse = TRUE)
    
    if (m > 0) {
      ix <- integer(m)
      jx <- integer(m)
      for (i in 1:m) {
        ix[i] <- which(Phi[connected_ix[i],]==1)
        jx[i] <- which(Phi[connected_ix[i],]==-1)
      }
      A[(jx-1)*n + ix] <- 1
    }
    return(A)
  }
  
  biclust_smooth <- function(X,clusters_row,clusters_col) {
    p <- nrow(X); n <- ncol(X)
    Y <- matrix(NA,p,n)
    M <- get_subgroup_means_full(X,clusters_row,clusters_col)
    num_clusters_row <- length(clusters_row$size)
    num_clusters_col <- length(clusters_col$size)
    for (i in 1:num_clusters_row) {
      ixi <- which(clusters_row$cluster == i)
      for (j in 1:num_clusters_col) {
        ixj <- which(clusters_col$cluster == j)
        Y[ixi,ixj] <- M[i,j]
      }
    }
    return(Y)
  }
  
  if(method == 'ADMM'){
    clusters_row = find_clusters(create_adjacency3(V1,t(E1)))
    clusters_col = find_clusters(create_adjacency3(Z1,t(E2)))
  }
  if(method == 'AMA'){
    clusters_row = find_clusters(create_adjacency4(V1,t(E1)))
    clusters_col = find_clusters(create_adjacency4(Z1,t(E2)))
  }
  
  M = biclust_smooth(X,clusters_row,clusters_col)
  
  return(list(M = M, clusters_row = clusters_row, clusters_col = clusters_col))
}

