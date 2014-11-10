#' Efficient Frontier function
#' 
#' Takes in any portfolio of stocks return and optimizes it using the Efficient Frontier 
#' @param Data A portfolio with its parameters to be optimize 
#' @return The allocation percent for the portfolio
#' @export
eff_frontier <- function (Data){
  # return argument should be a m x n matrix with one column per security
  # short argument is whether short-selling is allowed; default is no (short
  # selling prohibited)max_allocation is the maximum % allowed for any one
  # security (reduces concentration) risk_premium_up is the upper limit of the
  # risk premium modeled (see for loop below) and risk_increment is the
  # increment (by) value used in the for loop
  
  # Depend on the quadprog library to be running
  library("quadprog")
  
  #Convert the json data to numeric matrix
  data <- data.matrix(Data$data,rownames.force=TRUE) # Converts to data matrix
  colNames <- rownames(data) # Gets the rows names
  data <- as.data.frame(sapply(data, as.numeric)) # Convert all value from char to numeric
  names(data)<-colNames # Add the col names 
  data <- data.matrix(data) # Convert the data frame to a numeric matrix for processing
  
  #find and remove all the dimensions
  dimensionCount<-length(which(Data$type == 'dimension')) # Count # of dimensions
    
  returns<-data[,dimensionCount*-1] # remove dimensions from matrix
  
  #extract, convert and save all parameters
  short<-as.numeric(Data$param$short)
  max_allocation<-as.numeric(Data$param$max_allocation)
  risk_premium_up<-as.numeric(Data$param$risk_premium_up)
  risk_increment<-as.numeric(Data$param$risk_increment)
  
  # Let's Optimize!
  covariance <- cov(returns)
  
  n <- ncol(covariance)
  
  # Create initial Amat and bvec assuming only equality constraint
  # (short-selling is allowed, no allocation constraints)
  Amat <- matrix (1, nrow=n)
  bvec <- 1
  meq <- 1
  
  # Then modify the Amat and bvec if short-selling is prohibited
  # 0 means no short-selling
  # 1 means yes short-selling
  if(short==0){
    Amat <- cbind(1, diag(n))
    bvec <- c(bvec, rep(0, n))
  }
  
  # And modify Amat and bvec if a max allocation (concentration) is specified
  if(!is.null(max_allocation)){
    if(max_allocation > 1 | max_allocation <0){
      stop("max_allocation must be greater than 0 and less than 1")
    }
    if(max_allocation * n < 1){
      stop("Need to set max_allocation higher; not enough assets to add to 1")
    }
    Amat <- cbind(Amat, -diag(n))
    bvec <- c(bvec, rep(-max_allocation, n))
  }
  
  # Calculate the number of loops
  loops <- risk_premium_up / risk_increment + 1
  loop <- 1
  
  # Initialize a matrix to contain allocation and statistics
  # This is not necessary, but speeds up processing and uses less memory
  eff <- matrix(nrow=loops, ncol=n+3)
  # Now I need to give the matrix column names
  colnames(eff) <- c(colnames(returns), "Std.Dev", "Exp.Return", "sharpe")
  
  # Loop through the quadratic program solver
  for (i in seq(from=0, to=risk_premium_up, by=risk_increment)){
    dvec <- colMeans(returns) * i # This moves the solution along the EF
    sol <- solve.QP(covariance, dvec=dvec, Amat=Amat, bvec=bvec, meq=meq)
    eff[loop,"Std.Dev"] <- sqrt(sum(sol$solution*colSums((covariance*sol$solution))))
    eff[loop,"Exp.Return"] <- as.numeric(sol$solution %*% colMeans(returns))
    eff[loop,"sharpe"] <- eff[loop,"Exp.Return"] / eff[loop,"Std.Dev"]
    eff[loop,1:n] <- sol$solution
    loop <- loop+1
  }
  
  eff <-as.data.frame(eff)
  
  #Asuumption: check if the more then one result is max for Sharpe then take the first one
  if (nrow(eff[eff$sharpe==max(eff$sharpe),])>1) {
    optimized<-eff[eff$sharpe==max(eff$sharpe),!(names(eff) %in% c("Std.Dev", "Exp.Return", "sharpe"))]
    optimized<-optimized[1,]
  } else {
    optimized<-eff[eff$sharpe==max(eff$sharpe),!(names(eff) %in% c("Std.Dev", "Exp.Return", "sharpe"))]
  }
  
  return(optimized)
}