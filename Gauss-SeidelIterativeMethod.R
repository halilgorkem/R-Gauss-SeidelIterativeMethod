A <- matrix(c(10,-1,2,0,
              -1,11,-1,3,
              2,-1,10,-1,
              0,3,-1,8), nrow = 4, byrow = T)

b <-c(6,25,-11,15)

GaussSeidel <-  function(A, b, TOL = 1e-5, maxIter = 100)
{
  n <- nrow(A)
  iter <- 1
  x <- rep(0,n)
  x0 <- rep(0,n)
  X <- data.frame(c(x0))
  #step2
  while (iter <= maxIter) {
    
    # Adım 3
    x[1] <- (b[1] - sum(A[1,2:n] * x0[2:n])) / A[1,1]
    
    for(i in 2:(n-1)){
      x[i] <- (b[i] - sum(A[i,1:(i-1)] * x[1:(i-1)]) - sum(A[i,(i+1):n] * x0[(i+1):n])) / A[i,i]
    }
    
    x[n] <- (b[n] - sum(A[n,1:(n-1)] * x[1:(n-1)])) / A[n,n]
    
    X[,iter+1] <- x
    
    # Adım 4
    if (norm((x-x0), type = "2") < TOL){
      print("Algoritma başarılı...")
      return(list(x = as.vector(x), X, NbIter = iter)) 
    }
    
    # Adım 5
    iter <- iter + 1
    
    # Adım 6
    x0 <- x
  }
}
GaussSeidel(A,b)

