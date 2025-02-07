
# 1. Pre-Processing  ------------------------------------------------------
pre_processing_function <- function(data){
  #'  Function to restructure Bolding and Franks Data
  #'
  #' @param data (rds) rds file containing raw data
  #' @export  
  #' Outputs a nested list of trial spike times for each neuron (same structure as simulation study)

  h=vector(mode="list", length=dim(data)[1])
  
  trial_names = NULL;
  for(i in 1:dim(data)[1]){
    trial_names[i] = paste0("trial ", i)
  }
  names(h) = trial_names 
  data_output = rep(list(h), dim(data)[2]-1)
  
  neuron_names = NULL;
  for(i in 1:dim(data)[2]-1){
    neuron_names[i] = paste0("Neuron ", i)
  }
  
  names(data_output) = neuron_names
  
  
  for(i in 2:dim(data)[2]){ # remove first column as it is aggregate of the rest
    for(j in 1:dim(data)[1]){
      data_output[[i-1]][[j]] = (data[,i][[j]][1,])
    }
  }
  
  return(data_output)
}



# 2. Firing rate plots ----------------------------------------------------

spike_density = function(data, bin_size, a, b, title, c, d){
  #'  Function to plot firing rate plots from Bolding and Franks Data
  #'
  #' @param data pre-processed data in same format as simulation study
  #' @param bin_size time bins for data
  #' @param a,b xlimits for plots
  #' @param title title of plot
  #' @param c,d ylimits for plots 
  #' @export  
  #' Outputs a nested list of trial spike times for each neuron (same structure as simulation study)
  xs = seq(-5,10, bin_size)
  check = lapply(data, function(x) cut(x, breaks = xs))
  check2 = lapply(check, table)
  counts = NULL;
  for(i in 1:(length(xs)-1)){
    counts[i] = mean(unlist(map(check2, i)))
  }
  plot(seq(-5,(10-bin_size),bin_size), counts/bin_size, type='l', 
       xlab = "Time (s)", ylab = "Firing Rate (Hz)", main = title, xlim = c(a,b), ylim = c(c,d),
       cex.main = 2, cex.axis = 2, cex.lab = 2, family = "serif", las =1)
  rect(xleft=0, xright=1, ybottom = par("usr")[3], ytop = par("usr")[4], border = NA,col=adjustcolor("red", alpha=0.3))
  rect(xleft=1, xright=9, ybottom = par("usr")[3], ytop = par("usr")[4], border = NA,col=adjustcolor("blue", alpha=0.3))
  return(counts)
}


# 3. Trial-Frequency Smoothed Periodogram Estimate ------------------------

Trial_freq_Estimator = function(om, data, Big_T){
  #'  Function to estimate trial-frequency smoothed periodogram estimate
  #'
  #' @param om vector of frequencies 
  #' @param data data of interest in format of nested list of times
  #' @param Big_T length of time window 
  #' @export  
  S_hat = lapply(om, function(x) periodogram(data, x, Big_T)) #Get periodogram estimates for each frequecy of interest
  bands = list(seq(1,4), seq(4,8), seq(8,12), seq(12,30), seq(30,50))
  pooled_s = list();
  for(i in 1:length(bands)){
    p = S_hat[bands[[i]]]
    pooled_s[[i]] = Reduce("+", p)/length(p)
  }
  return(list(pooled_s = pooled_s, f = lengths(bands)))
}



# 4. Glasso estimator -----------------------------------------------------

glasso = function(n.stream,  chosen_lambda, S_hat, Max_iter){
  #'  Function to obtain glasso estimate for a given scenario and lambda value
  #'
  #' @param n.stream no. of data streams/neurons
  #' @param chosen_lambda regularization parameter 
  #' @param S_hat multi-trial periodogram estimate
  #' @param Max_iter Max no. of iterations for the ADMM algorithm
  #' @export  
  
  theta_warm = diag(1, nrow=n.stream, ncol=n.stream)
  res = ADMM(S_hat, Max_iter, chosen_lambda, rho=1, theta_warm)
  theta = res$theta
  r = res$r
  z = res$z
  
  
  
  return(list(theta = theta, z=z, r=r))
}




# 5. Function to evaluate eBIC --------------------------------------------

eBIC = function(theta, S, z, gam, n.stream, f, n.trials){
  #'  Function to evaluate eBIC 
  #'  #'
  #' @param theta estimate of inverse SDM for given lambda 
  #' @param S Peirodogram estimate
  #' @param gam tuning parameter
  #' @param n.stream number of dimensions/neurons
  #' @param f number of frequencies considered in trial-frequency smoothed estimator
  #' @param n.trials number of experimental trials (m)
  #' @export  
 
  lik = -log(Det(theta))+tr(S%*%theta)
  u = z[upper.tri(z, diag = F)]
  edges = length(u[u!=0])
  
  op = 2*n.trials*f*lik +  (edges *log(f*n.trials)) + 4*edges*gam*log(n.stream)
  
  return(op)
} 

# 6. ADMM Function --------------------------------------------------------

ADMM = function(S_hat, Max_iter,lambda,rho,theta_warm){
  #' ADMM algorithm for the estimation of the RSE with the l1 penalty.
  #'
  #' This function implements the Alternating Direction Method of Multipliers algorithm (ADMM) described in the supplementary material for the estimation of the RSE with the l1 penalty.
  #' @param S_hat Multi-taper periodogram estimator.
  #' @param Max_iter Maximum number of iterations for the ADMM algorithm.
  #' @param lambda Regularisation parameter (scalar).
  #' @param rho Augmented Lagrangian parameter
  #' @param theta_warm Initial value for theta matrix
  #' @returns A list with components
  #' \item{theta}{Estimate of Inverse spectral density matrix.}
  #' \item{z}{Z estimate from the ADMM algorithm.}
  #' \item{r}{ADMM stopping criteria.}
  #' @export

  
  
  #Initialization 
  
  theta = z = u = RHS_optim_cond = Q = Q_H = D = eigen_D = W = D_tilde = list();
  
  dim = dim(S_hat)[1] #How many processes are in total?
  
  
  
  r = s = NULL; #r2=s2=NULL;
  
  theta[[1]] = theta_warm; #previous solution to ADMM
  
  
  
  z[[1]] = u[[1]]= matrix(0,nrow=dim,ncol=dim)
  
  
  
  
  e_tol = 1e-4
  
  for(k in 2:Max_iter){
    
    
    ###############################
    ####Theta minimization step####
    ###############################
    
    RHS_optim_cond[[k]] = rho*(z[[k-1]]-u[[k-1]])-S_hat #RHS of optimality condition 
    
    Q[[k]] = eigen(RHS_optim_cond[[k]])$vectors #eigen vectors of above matrix
    
    Q_H[[k]] = H(Q[[k]])
    
    D[[k]] = diag(eigen(RHS_optim_cond[[k]])$values) #D-matrix - eigenvalues of RHS optim cond on diagonal
    
    #populate D-tilde matrix - does same as below loop
    D_tilde[[k]] = matrix(0,nrow=dim,ncol=dim)
    diag(D_tilde[[k]]) = (diag(D[[k]]) + sqrt(diag(D[[k]])^2+4*rho))/(2*rho) #populate D-tilde matrix - does same as below loop
    
    theta[[k]] = Q[[k]] %*% D_tilde[[k]] %*% Q_H[[k]]
    
    ###############################
    ##### Z-minimization step######
    ###############################
    
    kappa = lambda/rho
    
    W[[k]] = theta[[k]] + u[[k-1]] #W[[k]] - symmetric matrix?
    
    z[[k]] = pmax((1-kappa/Mod(W[[k]])),0) #checks if entries in 1-kappa/Mod(W[[k]]) < 0 -> if so replaces with a zero
    
    z[[k]] = z[[k]] * W[[k]] #want entry wise multiplication NOT matrix multiplication
    
    
    
    ###############################
    ###########U-update############
    ###############################
    
    u[[k]] = u[[k-1]]+theta[[k]]-z[[k]]
    
    ###############################
    ####### ADMM Residuals ########
    ###############################
    
    #Need to check
    
    d1 = theta[[k]] - z[[k]]
    d2 = z[[k]] - z[[k-1]]
    
    r[k] = sum(abs(d1))
    s[k] = sum(abs(d2))
    
    
    if(min(r[k], s[k])<e_tol & min(s, na.rm = T)!=0){break}
    
  }
  
  return(list(theta = theta[[k]], z = z[[k]], theta_all = theta, z_all = z, r = r,s=s, u=u))
  
}


# 7. Partial coherence ----------------------------------------------------

partial_co = function(theta){
  #' Estimation of the partial coherence matrix.
  #'
  #' This function is used to obtain the partial coherence matrix from an estimate of or the theoretical inverse spectral density matrix.
  #' @param theta p x p inverse spectral density matrix.
  #' @return p x p partial coherence matrix.
  #' @export
  
  #get indices to use - so it generalizes to higher dimensions
  indexes = which(upper.tri(theta), arr.ind = T) 
  index_list = list();
  
  for(i in 1:(length(indexes)/2)){
    index_list[[i]] = indexes[i,]
  }
  co = diag(1, nrow=dim(theta)[1], ncol=dim(theta)[2]);
  for(i in 1:length(index_list)){
    k = as.numeric(index_list[[i]][1])
    j = as.numeric(index_list[[i]][2])
    co[k,j] = Mod(-theta[k,j]/sqrt(theta[k,k]*theta[j,j]))^2
  }
  co[lower.tri(co)] = t(co)[lower.tri(co)] #populate lower triangle
  
  
  return(co)
}


# 8. Function to extract result from estimation procedure -----------------

calc_result = function(res){
  #'  Function to extract results from estimation procedure 
  #'  #'
  #' @param res rds file calcaulted in Estimation.R file
  #' @export  
  # get ebic estimates
  ebic = unlist(map(res, 1))
  idx = which.min(Re(ebic))
  # get partial coherence associate with this lambda
  z = map(res, 3)
  pc = partial_co(z[[idx]])
  # extract admm criteria for this etsimator
  r = map(res, 4)
  return(list(ebic = ebic, idx = idx, pc=pc, r=r[[idx]]))
  
}

# 9. Function to plot ADMM Stopping Criteria -----------------------------

plot_ADMM = function(tol, title){
  #'  Function to to plot ADMM Stopping Criteria
  #'  
  #' @param tol list of stopping criteria vectors for each scenario 
  #' @param title title of the plot (string)
  #' @export  
  par(mfrow=c(1,1), mar = c(5.1,4.1,4.1,2.1), font.main = 1)
  plot(seq(1, length(tol[[1]])-1), tol[[1]][-1], type='l', xlab = "No. Iterations", ylab = "Stopping Criteria", 
       main = paste0("Laser ", title), family="serif", cex =1 , cex.lab =1, cex.axis =1, las =1, cex.main =1)
  cols = c("black", "red", "blue")
  for(i in 2:3){
    lines(seq(1, length(tol[[i]])-1), tol[[i]][-1], col=cols[i])
  }
  legend("topright", c( expression(0~mW/mm^2),  expression(10~mW/mm^2),  expression(50~mW/mm^2)), col=c("black", "red", "blue"), lty = c(1,1,1))
}


# 10. Function to plot eBIC curves ----------------------------------------

plot_eBIC = function(eBIC, title){
  #'  Function to to plot eBIC curves
  #'
  #' @param eBIC list of ebic for each scenario 
  #' @param title title of the plot (string)
  #' @export 
  lambdas = logspace(-3,1,100)
  par(mfrow=c(1,1), mar = c(5.1,4.1,4.1,2.1), font.main = 1)
  eBIC = lapply(eBIC, Re)
  a = min(unlist(lapply(eBIC, min))); b = max(unlist(lapply(eBIC, max)))
  plot(lambdas, eBIC[[1]], xlab = expression(lambda), ylab = "eBIC", type="l", main = paste0("Laser ",title), 
       ylim = c(a, b), family="serif", cex =1 , cex.lab =1, cex.axis =1, cex.main =1,  las =1)
  cols = c("black", "red", "blue")
  for(i in 2:3){
    lines(lambdas, eBIC[[i]], col=cols[i])
  }
  legend("bottomright", c( expression(0~mW/mm^2),  expression(10~mW/mm^2),  expression(50~mW/mm^2)), col=c("black", "red", "blue"), lty = c(1,1,1))
}



# 11. Plot networks -------------------------------------------------------

plot_UGM_new <- function(pc_list, title_list, points_list){
  #' Network plots for Graphical estimator.
  #'
  #' This function is used to produce network plots obtained via the Graphical RSE.
  #' @param pc_list  list of partial coherence estimates to be compared.
  #' @param titles vector of titles for the considered plots.
  #' @param points vector of points to annotate the network plots.
  #' @return Network plots.
  #' @export
  
  upper = lapply(pc_list, function(x) {t(x)[lower.tri(x, diag=T)]})
  
  g1 = lapply(pc_list, function(x) graph_from_adjacency_matrix(x, mode="undirected", weighted=T))
  
  edges = lapply(upper, function(x){x[x!=0]})
  
  
  df = lapply(edges, function(x) data.frame(x))
  
  df = lapply(upper, function(x) data.frame(x))
  
  c = lapply(upper, function(x) which(x!=0)) #indices of non zero entries
  common = Reduce(intersect, c) #indices of common edges
  
  #check for commonalities among laser settings
  u2 = upper[-1] #get rid of list containing no laser estimates
  
  c2 = lapply(u2, function(x) which(x!=0)) #indices of non zero entries
  common2 = Reduce(intersect, c2)
  
  for(i in 1:length(pc_list)){
    df[[i]]$col = "white"
      
    for(j in 1:length(df[[i]]$x)){
      if(df[[i]]$x[j]!=1 & df[[i]]$x[j]!=0){df[[i]]$col[j] = "darkgrey"}
    }
    
    
    for(k in 1:length(df[[i]]$col[common])){
      if(df[[i]]$col[common][k]!="white"){df[[i]]$col[common][k] = "red"}
    }
    
  }
  
  for(i in 2:length(pc_list)){
    for(k in 1:length(df[[i]]$col[common2])){
      if(df[[i]]$col[common2][k]!="white"&df[[i]]$col[common2][k]!="red"){df[[i]]$col[common2][k] = "blue"}
    }
  }
  library(dplyr)
  df1 = lapply(df, function(x) filter(x, x>0))
  #curves <- autocurve.edges(g1[[1]])
  par(mfrow=c(1,1), mar=c(1,1,1,1), font.main = 1)
  for(i in 1:length(pc_list)){
    plot(g1[[i]], layout=layout.circle, edge.color=df1[[i]]$col, 
         vertex.color="white", vertex.label.color = "black", main ="", cex=2,vertex.size=20, 
         vertex.label.size=10, vertex.label.cex=1.3, family ="serif")
    title(title_list[i], cex.main = 2,  adj =0, line=-1, family="serif")
    title(points_list[i], cex.main = 2, adj =0, line = -24, family="serif")
  }
  
  
}


# 12. Periodogram Function ------------------------------------------------

periodogram = function(data, omega, Big_T){
  #' Multi-trial Periodogram Function.
  #'
  #' This function is used to estimate the multi-taper periodogram for a specified frequency of interest.
  #' @param data Input data of specified structure. 
  #' @param omega Specified frequency of interest.
  #' @param Big_T Time horizon of the data.
  #' @return Periodogram estimate of dimension p x p
  #' @export
  
  img=sqrt(as.complex(-1)) #imaginary number i^2=-1
  H_omega = (Big_T/sqrt(Big_T)) * exp((-img*omega*Big_T)/2) * sinc((omega*Big_T)/2) #FT of Taper
  p = length(data) #dimensionality 
  m = length(data[[1]]) #trials
  
  N_t = NULL;
  N_t = rep(list(N_t), m); bar_d = N_t; I = list();
  
  for(j in 1:m){
    for(i in 1:p){
      N_t[[j]][i] = length(data[[i]][[j]])
      bar_d[[j]][i] = (1/sqrt(Big_T) * sum(exp(-img*omega*data[[i]][[j]])) ) - ((N_t[[j]][i]/Big_T) * H_omega )
    }
    
    I[[j]] = (outer((bar_d[[j]]), Conj(bar_d[[j]])))/(2*pi)
  }
  
  p1 = Reduce("+", I)/length(I)
  
  return(p1)
  
}
