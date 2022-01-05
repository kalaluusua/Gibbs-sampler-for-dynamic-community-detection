rm(list=ls())
## We consider an undirected graph without self-loops

## Function for log-likelihood related to Jth observation
loglike <- function(clusterassign,param,data,J,n) #here J means Jth observation
{
  ################################################################
  
  ## Input: clusterassign = clustering configuration, a n by 1 vector ##
  ##        param = probability matrix, a k by k matrix ##
  ##        data = the adjacency matrix, a n by n matrix ##
  ##        J = observation index ##
  ##        n = number of observations ##
  
  ## Output: log-likelihood related to Jth observation ##
  
  #################################################################
  clustersize = max(clusterassign)
  param = as.matrix(param)
  
  if (J==1) {result2 = 0
  for (ii in c((J+1):n))
  {
    result2 = result2 + data[J,ii]*log(param[clusterassign[J],clusterassign[ii]])+(1-data[J,ii])*log(1-param[clusterassign[J],clusterassign[ii]])
  }
  output = sum(result2)} else if (J==n){
    result = 0
    for (ii in c(1:(J-1)))
    {
      result = result + data[ii,J]*log(param[clusterassign[ii],clusterassign[J]])+(1-data[ii,J])*log(1-param[clusterassign[ii],clusterassign[J]])
    }
    output = sum(result)
  } else {
    result = 0
    for (ii in c(1:(J-1)))
    {
      result = result + data[ii,J]*log(param[clusterassign[ii],clusterassign[J]])+(1-data[ii,J])*log(1-param[clusterassign[ii],clusterassign[J]])
    }
    
    result2 = 0
    for (ii in c((J+1):n))
      
    {
      result2 = result2 + data[J,ii]*log(param[clusterassign[J],clusterassign[ii]])+(1-data[J,ii])*log(1-param[clusterassign[J],clusterassign[ii]])
    }
    output = sum(result)+sum(result2)}
  output
}


## function for sampler for SBM (main algorithm)
SBM_new <- function(data, data1, niterations, beta.a, beta.b, GAMMA, LAMBDA, nClusters, nLayers = 1)
{
  ## Model: A_{ij}|z,Q \sim Ber(Q_{z_i,z_j}) ##
  ##        Q_{rs} \sim Beta(beta.a,beta.b), r,s = 1,...,k ##
  ##        P(z_i = j) = \omega_j, j = 1,...,k ##
  ##        \omega \sim Dirichlet_k(GAMMA,...,GAMMA) ##
  
  
  ################################################################
  
  ## Input: data = a list of nLayers adjacency matrices, each of which is a n by n matrix ##
  ##        data1 = a list of nLayers upper triangles for the adjacency matrices, each of which is a n by n matrix ##
  ##        niterations = the total number of iterations in Gibbs sampler ##
  ##        beta.a, beta.b = hyperparameters for the prior on elements in Q matrix in Beta distribution ##
  ##        GAMMA = the parameter in Dirichlet distribution that controls the relative size of clusters ##
  ##        nClusters = the number of clusters ##
  
  
  ## Output: 
  ##         zout = clustering configuration, a n by 1 vector##
  ##         Qout = probability matrix, a k by k matrix ##
  
  #################################################################
  n = dim(data[[1]])[1]
  gamma <- GAMMA
  N=n ## n is the number of oberservations
  
  # initialization of clustering configuration
  clusterAssign <- c(sample(1:nClusters, size = nClusters, replace = FALSE),
                     sample(1:nClusters, size = n-nClusters, replace = TRUE))
  
  Q<-matrix(0, nClusters,nClusters)
  for (i in 1:nClusters){
    for (j in i:nClusters){
      Q[i,j] = rbeta(1,beta.a,beta.b)
      Q[j,i] = Q[i,j]
    }
  }
  
  # sample omega from Dirichlet distribution indirectly
  omega <- matrix(0, nClusters)
  for (i in 1:nClusters){
    omega[i] = rgamma(1,gamma)
  }
  omega <- omega/sum(omega)
  
  History <- vector("list", niterations)
  
  ##start Gibb's sampling
  for (iter in 1:niterations)
  {
    ## update z ##
    for (i in 1:n){
      #finding the probs for sampling process
      clusterProbs = sapply(1:nClusters, function(x) {
        clusterAssign_temp = clusterAssign
        clusterAssign_temp[i] = x
        pow = 0
        for (t in 1:nLayers){
          pow = pow + loglike(clusterAssign_temp,Q,data[[t]],i,n)
        }
        log(omega[x]) + pow
      })
      clusterProbs <- exp(clusterProbs - max(clusterProbs))
      clusterProbs <- clusterProbs/sum(clusterProbs)
      #choose the cluster number for ith observation
      cluster.i <- sample(1:nClusters, size = 1,
                          prob = clusterProbs)
      clusterAssign[i] <- cluster.i
    }
    ## update Q ##
    Q = matrix(0, nClusters,nClusters)
    AA = matrix(0,nClusters,nClusters)
    NN = matrix(0,nClusters,nClusters)
    for (r in 1:nClusters){
      for (s in r:nClusters)
      {
        for (t in 1:nLayers){
          AA[r,s] <- AA[r,s] + sum(data1[[t]][clusterAssign==r,clusterAssign==s]) + sum(data1[[t]][clusterAssign==s,clusterAssign==r]) - 
            (r==s)*sum(data1[[t]][clusterAssign==s,clusterAssign==r])
        }
        med = matrix(0,n,n)
        med[which(clusterAssign==r),which(clusterAssign==s)] = 1
        med1 = matrix(0,n,n)
        med1[which(clusterAssign==s),which(clusterAssign==r)] = 1
        NN[r,s] = sum(med*lower.tri(med)) + sum(med1*lower.tri(med1)) - (r==s)*sum(med1*lower.tri(med1))
        Q[r,s] = rbeta(1,AA[r,s]+beta.a,nLayers*NN[r,s]-AA[r,s]+beta.b)
        Q[s,r] = Q[r,s]
      }
    }
    ## update omega ##
    omega <- matrix(0, nClusters)
    for (i in 1:nClusters){
      omega[i] = rgamma(1,sum(clusterAssign==i) + gamma)
    }
    omega <- omega/sum(omega)
    
    History[[iter]] <- list(zout = clusterAssign,Qout = Q,omegaout = omega)
    # cat(" iteration:", iter,"\n",clusterAssign,"\n")
  }# for loop over iterations
  
  list(Iterates = History)
}

## Dahl's method to summarize the samples from the MCMC
getDahl <- function(data1, fit, burn, nClusters, nLayers = 1)
{
  ################################################################
  
  ## Input: fit = the result from SBM_new ##
  ##        burn = the number of burn-in iterations ##
  
  ## Output: 
  ##         zout = estimated clustering configuration, a n by 1 vector##
  ##         Qout = estimated probability matrix, a k by k matrix ##
  
  #################################################################
  iters <- fit$Iterates[-(1:burn)]
  n <- length(iters[[1]][[1]])
  niters <- length(iters)
  membershipMatrices <- lapply(iters, function(x){
    clusterAssign <- x[[1]]
    outer(clusterAssign, clusterAssign, FUN = "==")
  })
  membershipAverage <- Reduce("+", membershipMatrices)/niters
  print(membershipAverage)
  SqError <- sapply(membershipMatrices, function(x, av) sum((x - av)^2),
                    av = membershipAverage)
  DahlIndex <- which.min(SqError)
  DahlAns <- iters[[DahlIndex]]
  attr(DahlAns, "iterIndex") <- burn + DahlIndex
  attr(DahlAns, "burnin") <- burn
  ## get Q ##
  AA = matrix(0,nClusters,nClusters)
  NN = matrix(0,nClusters,nClusters)
  for (r in 1:nClusters){
    for (s in r:nClusters)
    {
      for (t in 1:nLayers){
        AA[r,s] <- AA[r,s] + sum(data1[[t]][DahlAns$zout==r,DahlAns$zout==s]) + sum(data1[[t]][DahlAns$zout==s,DahlAns$zout==r]) - 
          (r==s)*sum(data1[[t]][DahlAns$zout==s,DahlAns$zout==r])
      }
      med = matrix(0,n,n)
      med[which(DahlAns$zout==r),which(DahlAns$zout==s)] = 1
      med1 = matrix(0,n,n)
      med1[which(DahlAns$zout==s),which(DahlAns$zout==r)] = 1
      NN[r,s] = sum(med*lower.tri(med)) + sum(med1*lower.tri(med1)) - (r==s)*sum(med1*lower.tri(med1))
      DahlAns$Qout[r,s] = (1/nLayers)*AA[r,s]/NN[r,s]
      DahlAns$Qout[s,r] = DahlAns$Qout[r,s]
    }
  }
  ## get omega ##
  for (i in 1:nClusters){
    DahlAns$omegaout[i] = sum(DahlAns$zout==i)/n
  }
  DahlAns
}

# from package fossil
adj.rand.index <- 
  function(group1, group2) {
    a <- length(table(group1))
    N <- length(group1)
    ctab <- matrix(NA, a, a)
    for (j in 1:a) {
      for (i in 1:a) {
        ctab[j,i] <- length(which(group2[which(group1==i)]==j))
      }
    }
    sumnij <- sum(choose(ctab, 2))
    sumai <- sum(choose(colSums(ctab), 2))
    sumbj <- sum(choose(rowSums(ctab), 2))
    Ntwo <- choose(N, 2)
    ari <- abs((sumnij - (sumai*sumbj)/Ntwo)/(0.5*(sumai+sumbj)-(sumai*sumbj)/Ntwo))
    return(ari)
  }

NN <- 1 # repeat the simulation study for NN identically distributed networks
errors <- c()
# library(fossil) # contains rand.index
adj_rand_indices <- c()

for (s in 1:NN) {
  ###### one example in simulation study
  ## data generation
  # set.seed(33)
  n = 100 ## number of observations
  kk = 2 ## number of clusters
  Z <- c(sample(1:kk, size = kk, replace = FALSE),
         sample(1:kk, size = n-kk, replace = TRUE,prob = c(1,1))) ## clustering configuration
  Z = Z[order(Z)]
  theta <- matrix(0.1,kk,kk) ## off-diagonal value for Q matrix
  diag(theta) = 0.15 ## diagonal value for Q matrix
  NL = 5 # Number of layers
  AT = list() # List of iid adjacency matrices
  AAAT = list() # List of upper triangles of iid adjacency matrices
  
  for (t in 1:NL){
    A = matrix(0,n,n) ##the adjacency matrix
    AAA = matrix(0,n,n) ##the upper triangle for the adjacency matrix
    for (i in 1:n){
      for (j in i:n){
        A[i,j] = rbinom(1,1,prob=theta[Z[i],Z[j]])
        A[j,i] = A[i,j]
        AAA[i,j] = A[i,j]
      }
    }  
    diag(AAA) = 0 ## make it without-selfloop network
    diag(A) = 0 ## make it without-selfloop network
    AT[[t]] = A
    AAAT[[t]] = AAA
  }
  ## taking the data into the SBM algorithm
  # set.seed(1)
  fit1 = SBM_new(data = AT, data1 = AAAT, niterations = 120, beta.a = 1, beta.b = 1, GAMMA=100*kk, LAMBDA = 1, nClusters = 2, nLayers = NL)
  ## fit1$Iterates[[i]] is a list of length two, which denotes the ith sample in MCMC output.
  ## fit1$Iterates[[i]][[1]] denotes the clustering configuration z in ith iteration.
  ## fit1$Iterates[[i]][[2]] denotes the Q matrix in ith iteration.
  
  ## estimated configuration using Dahl's method, choosing first 50 iterations in MCMC as burn-in
  result1 = getDahl(data1 = AAAT, fit1, burn = 100, nClusters = 2, nLayers = NL)
  ## result1[[1]] denotes the estimated clustering configuration.
  ## result1[[2]] denotes the estimated Q matrix.
  errors[s] <- min(length(which(result1[[1]] != Z)), length(which(result1[[1]] != abs(Z - 3))))/n
  adj_rand_indices[s] <- adj.rand.index(result1[[1]], Z)
  cat("epoch:", s, "error:", errors[s], "rand (adj):", adj_rand_indices[s],"\n")
}
sum(errors)/NN
sum(adj_rand_indices)/NN

