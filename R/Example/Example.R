
# Example with R Package --------------------------------------------------

############ If not already installed ############ 

#library(devtools)
#install_github('pinkney19/RiSE')

library(RiSE)


# Required Libraries ------------------------------------------------------
library(Matrix)
library(phonTools) 
library(QZ) 
library(purrr) 
library(matrixStats) 
library(complexplus) 
library(psych) 
library(pracma) 
library(hawkes) 
library(igraph)



# Estimation --------------------------------------------------------------

# Laser off
library(pracma)
lambdas = logspace(log10(0.001), log10(10), 100) #grid of lambdas to search over
library(doParallel)
library(doRNG)
cl <- makeCluster(4) 
start <- Sys.time()
registerDoParallel(cl)
registerDoRNG(seed = 123)
res_off <- foreach(j=lambdas,  .packages = c("Matrix","phonTools", "QZ", "purrr", "matrixStats", "complexplus", 
                                         "psych", "pracma", "hawkes", "RiSE")) %dopar%{
                                           
                                           
                                           n.trials = 10
                                           n.stream = 26
                                           
                                           S_hat <- readRDS("~/Downloads/RiSE/R/S_hat_off.RDS")
                                           
                                           freqs <- readRDS("~/Downloads/RiSE/R/freqs.RDS")
                                           
                                           S_hat = S_hat[[1]] # 1 indicates delta band
                                           f = freqs[1]
                                           
                                           Max_iter = 1000
                                           
                                           out = glasso(n.stream, j, S_hat, Max_iter)
                                           
                                           theta = out$theta 
                                           zs = out$z
                                           
                                           gam = 0.5
                                           
                                           lik = eBIC(theta, S_hat, zs, gam, n.stream, f, n.trials)
                                           
                                           return(list(eBIC = lik, theta = theta, z = zs, r = out$r ) )
                                           
                                         }
end <-Sys.time() 
stopCluster(cl) 
print(end-start)

# Laser on 
library(pracma)
lambdas = logspace(log10(0.001), log10(10), 100) #grid of lambdas to search over
library(doParallel)
library(doRNG)
cl <- makeCluster(4) 
start <- Sys.time()
registerDoParallel(cl)
registerDoRNG(seed = 123)
res_on <- foreach(j=lambdas,  .packages = c("Matrix","phonTools", "QZ", "purrr", "matrixStats", "complexplus", 
                                             "psych", "pracma", "hawkes", "RiSE")) %dopar%{
                                               
                                               
                                               n.trials = 10
                                               n.stream = 26
                                               
                                               S_hat <- readRDS("~/Downloads/RiSE/R/S_hat_on.RDS")
                                              
                                               freqs <- readRDS("~/Downloads/RiSE/R/freqs.RDS") 
                                               
                                               S_hat = S_hat[[1]] # 1 indicates delta band
                                               f = freqs[1]
                                               
                                               Max_iter = 1000
                                               
                                               out = glasso(n.stream, j, S_hat, Max_iter)
                                               
                                               theta = out$theta 
                                               zs = out$z
                                               
                                               gam = 0.5
                                               
                                               lik = eBIC(theta, S_hat, zs, gam, n.stream, f, n.trials)
                                               
                                               return(list(eBIC = lik, theta = theta, z = zs, r = out$r ) )
                                               
                                             }
end <-Sys.time() 
stopCluster(cl) 
print(end-start)


# Plots -------------------------------------------------------------------

out_off = calc_result(res_off)
out_on = calc_result(res_on)

# ebic curves -------------------------------------------------------------

par(mfrow=c(1,2))
plot_eBIC(list(out_on$ebic), "On")
plot_eBIC(list(out_off$ebic), "Off")



plot_ADMM(list(out_on$r), "On")
plot_ADMM(list(out_ff$r), "Off")


# Partial coherence estimates ---------------------------------------------

# edges 
edges = function(mat){
  u = mat[upper.tri(mat, diag=F)]
  edges = length(u[u!=0])
  return(edges)
}

edges_list_on = c(edges(out_on$pc))
edges_list_off = c(edges(out_off$pc)) 


edges_list_on
edges_list_off


title_list = c("")


library(igraph)
ten = list(out_off$pc, out_on$pc)
ten_points = c(("AvgPoints = 778"), expression("AvgPoints = 712"))
plot_UGM_new(ten, title_list ,ten_points)
















eBIC_off = unlist(map(res_off, 1))
eBIC_on = unlist(map(res_on, 1))

par(mfrow =c (1,2))
plot(lambdas, eBIC_off, ylab ="eBIC", xlab = expression(lambda), main = "Laser Off")
plot(lambdas, eBIC_on, ylab ="eBIC", xlab = expression(lambda), main = "Laser On")

idx_off = which.min(eBIC_off)
idx_on = which.min(eBIC_on)

pc_off = partial_co(map(res_off, 3)[[idx_off]])
pc_on = partial_co(map(res_on, 3)[[idx_on]])

plot_UGM_new(list(pc_off))